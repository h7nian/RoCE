"""Fresh honest nuisance fitting in a finite-stratum causal diagnostic model.

Half-count regularized saturated outcome/propensity fits and positive calibrated
cell weights have closed forms. This is not the production sparse RoCE fitting
program. Exact conditional drifts and feasible product bounds are both retained.
"""
import numpy as np
from scipy.stats import beta
from repairs import DIRECTIONS


def binomial_bounds(successes, trials, failure):
    successes, trials = np.broadcast_arrays(successes, trials)
    lower = np.zeros(successes.shape)
    upper = np.ones(successes.shape)
    positive = successes > 0
    incomplete = successes < trials
    lower[positive] = beta.ppf(failure/2, successes[positive], trials[positive]-successes[positive]+1)
    upper[incomplete] = beta.ppf(1-failure/2, successes[incomplete]+1, trials[incomplete]-successes[incomplete])
    return lower, upper


def fitted_stratum_data(setting, repeats, bias_failure=.005, retain_source_counts=False):
    count = setting["source_count"]
    training, pilot_size, evaluation = 500, 250, 250
    target_evaluation = 500
    cells = np.array([(x1,x2) for x1 in [-1,1] for x2 in [-1,1]])
    cell_count = len(cells)
    states = np.array([(cell,arm,outcome) for cell in range(cell_count) for arm in [1,0] for outcome in [0,1]])
    index, arm_index, outcome = states[:,0], 1-states[:,1], states[:,2]
    target_mass = np.full(cell_count, .25)
    target_propensity = .5 + .1*cells[:,0] - .05*cells[:,1]
    heterogeneity = setting["shared_scale"]
    true_outcomes = np.column_stack((.55+heterogeneity*cells[:,0]/2+.08*cells[:,1],
                                     .45-heterogeneity*cells[:,0]/2+.08*cells[:,1]))
    def probabilities(mass, propensity, shift):
        outcome_probability = true_outcomes[index,arm_index]+shift[arm_index]
        treatment_probability = np.where(arm_index==0,propensity[index],1-propensity[index])
        if np.any((outcome_probability<=0)|(outcome_probability>=1)):
            raise ValueError("Invalid fitted-model outcome probabilities")
        return mass[index]*treatment_probability*np.where(outcome==1,outcome_probability,1-outcome_probability)
    target_probability = probabilities(target_mass,target_propensity,np.zeros(2))
    rng = np.random.RandomState(721038 + setting.get("seed_offset",0))
    target_train = rng.multinomial(training,target_probability,size=repeats).reshape(repeats,cell_count,2,2)
    target_eval = rng.multinomial(target_evaluation,target_probability,size=repeats)
    target_arm_counts = target_train.sum(axis=3)
    target_cell_counts = target_arm_counts.sum(axis=2)
    common = (target_train[:,:,:,1]+.5)/(target_arm_counts+1)
    target_assignment = (target_arm_counts[:,:,0]+.5)/(target_cell_counts+1)
    target_inverse = 1/np.stack((target_assignment,1-target_assignment),axis=2)
    fitted_target_mass = (target_cell_counts+.5)/(training+.5*cell_count)

    failure_per_bound = bias_failure/(4*cell_count+2*count*cell_count)
    mean_lower,mean_upper = binomial_bounds(target_train[:,:,:,1],target_arm_counts,failure_per_bound)
    outcome_error = np.maximum(abs(common-mean_lower),abs(common-mean_upper))
    mass_lower,mass_upper = binomial_bounds(target_cell_counts,training,failure_per_bound)
    mass_error = np.maximum(abs(fitted_target_mass-mass_lower),abs(fitted_target_mass-mass_upper))
    propensity_lower,propensity_upper = binomial_bounds(target_arm_counts[:,:,0],target_cell_counts,failure_per_bound)
    propensity_error = np.maximum(abs(target_assignment-propensity_lower),abs(target_assignment-propensity_upper))
    target_feasible = np.max(propensity_error[:,:,None]*target_inverse*outcome_error,axis=1)
    # Shared target nuisance events are paid for once, independently of K.
    target_failure = (bias_failure/2)/(4*cell_count)
    shared_mean_lower,shared_mean_upper = binomial_bounds(target_train[:,:,:,1],target_arm_counts,target_failure)
    shared_outcome_error = np.maximum(abs(common-shared_mean_lower),abs(common-shared_mean_upper))
    shared_mass_lower,shared_mass_upper = binomial_bounds(target_cell_counts,training,target_failure)
    shared_mass_error = np.maximum(abs(fitted_target_mass-shared_mass_lower),abs(fitted_target_mass-shared_mass_upper))
    shared_e_lower,shared_e_upper = binomial_bounds(target_arm_counts[:,:,0],target_cell_counts,target_failure)
    shared_e_error = np.maximum(abs(target_assignment-shared_e_lower),abs(target_assignment-shared_e_upper))
    shared_target_feasible = np.max(shared_e_error[:,:,None]*target_inverse*shared_outcome_error,axis=1)

    common_state = common[:,index,:]
    target_residual = np.zeros((repeats,len(states),2))
    for arm in [0,1]:
        target_residual[:,:,arm] = (arm_index==arm)[None,:]*target_inverse[:,index,arm]*(outcome[None,:]-common_state[:,:,arm])
    target_support = np.concatenate(((common_state[:,:,0]-common_state[:,:,1])[:,:,None],target_residual),axis=2)
    target = np.einsum("bs,bsi->bi",target_eval,target_support)/target_evaluation
    target_expected = np.einsum("s,bsi->bi",target_probability,target_support)
    target_second = np.einsum("bs,bsi,bsj->bij",target_eval,target_support,target_support)
    target_sample_covariance = (target_second-target_evaluation*np.einsum("bi,bj->bij",target,target))/(target_evaluation*(target_evaluation-1))
    target_population_second = np.einsum("s,bsi,bsj->bij",target_probability,target_support,target_support)
    target_covariance = (target_population_second-np.einsum("bi,bj->bij",target_expected,target_expected))/target_evaluation
    # Full-target two-fold AIPW reference uses the same 1000 target patients.
    # Its usual score-normal interval is an empirical comparator, not a finite
    # sample Bernstein certificate. It is not allowed to tune source inference.
    reverse_train = target_eval.reshape(repeats,cell_count,2,2)
    reverse_arm = reverse_train.sum(axis=3)
    reverse_cell = reverse_arm.sum(axis=2)
    reverse_common = (reverse_train[:,:,:,1]+.5)/(reverse_arm+1)
    reverse_assignment = (reverse_arm[:,:,0]+.5)/(reverse_cell+1)
    reverse_inverse = 1/np.stack((reverse_assignment,1-reverse_assignment),axis=2)
    reverse_score = reverse_common[:,index,0]-reverse_common[:,index,1]
    for arm,sign in [(0,1),(1,-1)]:
        reverse_score += sign*(arm_index==arm)[None,:]*reverse_inverse[:,index,arm]*(outcome[None,:]-reverse_common[:,index,arm])
    first_counts = target_train.reshape(repeats,len(states))
    forward_score = target_support@DIRECTIONS[3]
    total_target = training+target_evaluation
    reference_point = (np.sum(first_counts*reverse_score,axis=1)+np.sum(target_eval*forward_score,axis=1))/total_target
    reference_second = np.sum(first_counts*reverse_score**2,axis=1)+np.sum(target_eval*forward_score**2,axis=1)
    reference_variance = (reference_second-total_target*reference_point**2)/(total_target*(total_target-1))
    # A fair finite-sample target-only comparison. Each half receives .0025
    # nuisance-error probability and .0225 held-out mean-error probability.
    # Dependence between the two fitted folds is handled by a union bound.
    certified_half_radii = []
    reference_noise_log = np.log(4/.0225)
    for train_counts,eval_counts,prediction,assignment,score in [
        (target_train,target_eval,common,target_assignment,forward_score),
        (reverse_train,first_counts,reverse_common,reverse_assignment,reverse_score)]:
        arm_counts = train_counts.sum(axis=3)
        cell_counts = arm_counts.sum(axis=2)
        ml,mu = binomial_bounds(train_counts[:,:,:,1],arm_counts,.0025/(3*cell_count))
        el,eu = binomial_bounds(arm_counts[:,:,0],cell_counts,.0025/(3*cell_count))
        outcome_deviation = np.maximum(abs(prediction-ml),abs(prediction-mu))
        propensity_deviation = np.maximum(abs(assignment-el),abs(assignment-eu))
        inverse = 1/np.stack((assignment,1-assignment),axis=2)
        remainder = np.max(outcome_deviation*propensity_deviation[:,:,None]*inverse,axis=1).sum(axis=1)
        average = np.sum(eval_counts*score,axis=1)/500
        score_variance = (np.sum(eval_counts*score**2,axis=1)-500*average**2)/(500*499)
        radius = np.sqrt(2*np.maximum(score_variance,0)*reference_noise_log) + \
            7*np.ptp(score,axis=1)*reference_noise_log/(3*499) + remainder
        certified_half_radii.append(radius)
    certified_radius = (certified_half_radii[0]+certified_half_radii[1])/2
    certified_reference_interval = np.clip(np.column_stack((reference_point-certified_radius,
                                                            reference_point+certified_radius)),-1,1)
    certified_reference_width = np.diff(certified_reference_interval,axis=1)[:,0]
    common_population = np.mean(common,axis=1)
    residual_truth = np.mean(true_outcomes,axis=0)-common_population
    target_truth = np.column_stack((common_population[:,0]-common_population[:,1],residual_truth))
    target_oracle_budget = abs(target_expected-target_truth)

    source = np.empty((repeats,count,2))
    variance = np.empty_like(source)
    sample_variance = np.empty_like(source)
    pilot_variance = np.empty_like(source)
    source_ranges = np.empty_like(source)
    source_absolute_bounds = np.empty_like(source)
    source_weights = np.empty((repeats, count, cell_count, 2))
    pilot_joint_frequencies = np.empty_like(source_weights)
    evaluation_joint_frequencies = np.empty_like(source_weights)
    pilot_source_mean = np.empty_like(source)
    source_training_mean = np.empty_like(source)
    source_score_lower = np.empty_like(source)
    source_score_upper = np.empty_like(source)
    source_contrast_lower = np.empty((repeats, count))
    source_contrast_upper = np.empty((repeats, count))
    source_sample_cross_covariance = np.empty((repeats, count))
    pilot_sample_cross_covariance = np.empty((repeats, count))
    nuisance_drift = np.empty_like(source)
    source_expected = np.empty_like(source)
    feasible_source = np.empty_like(source)
    shared_feasible_source = np.empty_like(source)
    structural_shift = np.zeros((count,2))
    if setting["pattern"]=="weak_boundary":
        structural_shift[np.arange(count)%4==0,0]=2*count**(-.25)*np.sqrt(.5/evaluation)
    elif setting["pattern"]=="weak_disjoint":
        shift=2*count**(-.25)*np.sqrt(.5/evaluation)
        structural_shift[np.arange(count)%4==0,0]=shift
        structural_shift[np.arange(count)%4==1,1]=-shift
    elif setting["pattern"]=="arm_cancel":
        structural_shift[np.arange(count)%4==0,:]=.2
    elif setting["pattern"]=="strong":
        structural_shift[np.arange(count)%4==0,0]=.2
    elif setting["pattern"]!="all_valid":
        raise ValueError("Unknown fitted-model bias pattern")
    count_data = {}
    if retain_source_counts:
        count_data['source_training_counts'] = np.empty((repeats,count,cell_count,2,2),dtype=np.int16)
        count_data['source_held_counts'] = np.empty((repeats,count,cell_count,2,2),dtype=np.int16)
    source_rng = np.random.RandomState(831017+setting.get("seed_offset",0))
    for site in range(count):
        orientation = 1 if site%2==0 else -1
        mass = np.exp(.3*orientation*cells[:,0]);mass/=mass.sum()
        propensity = .5+.1*orientation*cells[:,0]+.05*cells[:,1]
        probability = probabilities(mass,propensity,structural_shift[site])
        train_counts = source_rng.multinomial(training,probability,size=repeats).reshape(repeats,cell_count,2,2)
        joint_counts = train_counts.sum(axis=3)
        fitted_joint = (joint_counts+.5)/(training+cell_count)
        weights = fitted_target_mass[:,:,None]/fitted_joint
        source_weights[:, site] = weights
        joint_lower,joint_upper = binomial_bounds(joint_counts,training,failure_per_bound)
        joint_error = np.maximum(abs(fitted_joint-joint_lower),abs(fitted_joint-joint_upper))
        feasible_source[:,site] = np.sum((weights*joint_error+mass_error[:,:,None])*outcome_error,axis=1)
        shared_jl,shared_ju = binomial_bounds(joint_counts,training,(bias_failure/2)/(2*count*cell_count))
        shared_joint_error = np.maximum(abs(fitted_joint-shared_jl),abs(fitted_joint-shared_ju))
        shared_feasible_source[:,site] = np.sum((weights*shared_joint_error+shared_mass_error[:,:,None])*shared_outcome_error,axis=1)
        support = np.zeros((repeats,len(states),2))
        for arm in [0,1]:
            support[:,:,arm]=(arm_index==arm)[None,:]*weights[:,index,arm]*(outcome[None,:]-common_state[:,:,arm])
        source_ranges[:,site]=np.ptp(support,axis=1)
        source_score_lower[:,site] = support.min(axis=1)
        source_score_upper[:,site] = support.max(axis=1)
        contrast_support = support[:,:,0] - support[:,:,1]
        source_contrast_lower[:,site] = contrast_support.min(axis=1)
        source_contrast_upper[:,site] = contrast_support.max(axis=1)
        source_training_mean[:,site] = np.einsum('bs,bsa->ba', train_counts.reshape(repeats,-1), support) / training
        source_absolute_bounds[:,site]=np.max(abs(support),axis=1)
        true_joint = mass[:,None]*np.column_stack((propensity,1-propensity))
        null_expected = np.sum(true_joint[None,:,:]*weights*(true_outcomes[None,:,:]-common),axis=1)
        nuisance_drift[:,site] = null_expected-residual_truth
        expected = np.einsum("s,bsi->bi",probability,support)
        source_expected[:,site] = expected
        second = np.einsum("s,bsi->bi",probability,support**2)
        variance[:,site] = (second-expected**2)/evaluation
        if retain_source_counts:
            count_data['source_training_counts'][:,site] = train_counts
            held_counts = np.zeros_like(train_counts)
        for sample_size,kind in [(pilot_size,"pilot"),(evaluation,"evaluation")]:
            counts = source_rng.multinomial(sample_size,probability,size=repeats)
            if retain_source_counts:
                held_counts += counts.reshape(repeats,cell_count,2,2)
            means = np.einsum("bs,bsi->bi",counts,support)/sample_size
            second = np.einsum("bs,bsi->bi",counts,support**2)
            estimates = (second-sample_size*means**2)/((sample_size-1)*evaluation)
            # R^1 R^0 is zero for each patient; its centered covariance need not be.
            cross_covariance = -sample_size * means[:,0] * means[:,1] / ((sample_size-1) * evaluation)
            if kind=="pilot":
                pilot_sample_cross_covariance[:,site] = cross_covariance
                pilot_variance[:,site]=estimates
                pilot_source_mean[:,site]=means
                pilot_joint_frequencies[:,site]=counts.reshape(repeats,cell_count,2,2).sum(axis=3)/sample_size
            else:
                source[:,site]=means
                source_sample_cross_covariance[:,site] = cross_covariance
                evaluation_joint_frequencies[:,site]=counts.reshape(repeats,cell_count,2,2).sum(axis=3)/sample_size
                sample_variance[:,site]=estimates
        if retain_source_counts:
            count_data['source_held_counts'][:,site] = held_counts
    exact_source_budget = np.max(abs(nuisance_drift),axis=1)
    feasible_source_budget = np.max(feasible_source,axis=1)
    return dict(nuisance_structure='saturated_design_weights',
        source=source,variance=variance,sample_variance=sample_variance,pilot_variance=pilot_variance,
        source_training_mean=source_training_mean,
        source_score_lower=source_score_lower,source_score_upper=source_score_upper,
        source_contrast_lower=source_contrast_lower,source_contrast_upper=source_contrast_upper,
        source_sample_cross_covariance=source_sample_cross_covariance,
        pilot_sample_cross_covariance=pilot_sample_cross_covariance,
        target_prediction_bounds=np.column_stack(((common[:,:,0]-common[:,:,1]).min(axis=1),
                                                 (common[:,:,0]-common[:,:,1]).max(axis=1))),
        target_training_prediction_mean=np.sum(target_cell_counts*(common[:,:,0]-common[:,:,1]),axis=1)/training,
        score_absolute_bounds=source_absolute_bounds,
        source_weights=source_weights,pilot_joint_frequencies=pilot_joint_frequencies,
        evaluation_joint_frequencies=evaluation_joint_frequencies,pilot_source_mean=pilot_source_mean,
        target_training_counts=target_train,target_evaluation_counts=reverse_train,
        common_prediction=common,target_inverse=target_inverse,
        target_full_joint_counts=target_arm_counts+reverse_arm,
        score_ranges=source_ranges,sample_sizes=np.full(count,evaluation),
        target=target,target_covariance=target_covariance,target_sample_covariance=target_sample_covariance,
        target_direction_ranges=np.ptp(target_support@DIRECTIONS.T,axis=1),target_sample_size=target_evaluation,
        reference_point=reference_point,reference_width=2*1.959963984540054*np.sqrt(reference_variance),
        certified_reference_width=certified_reference_width,
        certified_reference_interval=certified_reference_interval,
        truth=.1,residual_truth=residual_truth,target_truth=target_truth,
        oracle_source_budget=exact_source_budget,feasible_source_budget=feasible_source_budget,
        oracle_target_budget=target_oracle_budget,
        feasible_target_budget=np.column_stack((np.zeros(repeats),target_feasible)),
        shared_source_budget=np.max(shared_feasible_source,axis=1),
        shared_target_budget=np.column_stack((np.zeros(repeats),shared_target_feasible)),
        nuisance_drift=nuisance_drift,structural_shift=structural_shift,
        source_expected=source_expected,parameter_bounds=(-1.,1.),
        pilot_size=pilot_size,source_evaluation_size=evaluation,
        source_training_size=training,target_training_size=training,**count_data)
