"""Independent patient-level checks for the joint contrast procedure."""
import itertools
import unittest
import numpy as np
from scipy.stats import binom
from joint_contrast import weighted_empirical_mean_radius, infer_joint_contrast, nuisance_radius
from fitted_strata import fitted_stratum_data
from conditional_remainders import reuse_source_design_for_evaluation, conditional_design_coefficients


class JointContrastChecks(unittest.TestCase):
    def test_weighted_noise_matches_scaled_patient_bound(self):
        rng = np.random.RandomState(3189)
        groups = [rng.uniform(-2, 3, (8, n)) for n in [7, 11, 13]]
        weights = np.array([1., .4, -.2])
        centers = np.array([.3, -.1, .2])
        sizes = np.array([g.shape[1] for g in groups])
        means = np.column_stack([g.mean(axis=1) for g in groups])
        variances = np.column_stack([g.var(axis=1,ddof=1)/g.shape[1] for g in groups])
        radius, variance = weighted_empirical_mean_radius(means,variances,sizes,weights,-2.,3.,centers,.05)
        total = sizes.sum()
        scaled = np.concatenate([total*weights[j]*(g-centers[j])/sizes[j] for j,g in enumerate(groups)],axis=1)
        direct_variance = scaled.var(axis=1,ddof=1)/total
        lower = np.minimum(total*weights*(-2-centers)/sizes,total*weights*(3-centers)/sizes).min()
        upper = np.maximum(total*weights*(-2-centers)/sizes,total*weights*(3-centers)/sizes).max()
        direct_radius = np.sqrt(2*direct_variance*np.log(80))+7*(upper-lower)*np.log(80)/(3*(total-1))
        np.testing.assert_allclose(variance,direct_variance,rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(radius,direct_radius,rtol=1e-12,atol=1e-12)

    def test_nonidentical_weighted_binary_exact_enumeration(self):
        observations = np.array(list(itertools.product([0.,1.],repeat=6)))
        probabilities = np.array([.2,.2,.7,.7,.7,.7])
        law = np.prod(np.where(observations==1,probabilities,1-probabilities),axis=1)
        groups = [observations[:,:2],observations[:,2:]]
        means = np.column_stack([g.mean(axis=1) for g in groups])
        variances = np.column_stack([g.var(axis=1,ddof=1)/g.shape[1] for g in groups])
        for alpha in [.05,.2]:
            radius,_ = weighted_empirical_mean_radius(means,variances,[2,4],[1.,.5],0.,1.,[.1,.6],alpha)
            error = abs(means@np.array([1.,.5])-(.2+.5*.7))
            self.assertLessEqual(law@(error>radius),alpha+1e-12)

    def test_source_arm_covariance_and_reuse(self):
        data = fitted_stratum_data(dict(source_count=8,shared_scale=.3,pattern='strong',seed_offset=4141),30)
        for packet in [data,reuse_source_design_for_evaluation(data)]:
            n = packet['source_evaluation_size']
            expected_cross = -np.prod(packet['source'],axis=2)/(n-1)
            np.testing.assert_allclose(packet['source_sample_cross_covariance'],expected_cross,atol=1e-12)
            self.assertGreater(np.max(abs(expected_cross)),0)
            var = packet['sample_variance'].sum(axis=2)-2*expected_cross
            raw_second = packet['sample_variance'].sum(axis=2)*(n-1)/n+np.sum(packet['source']**2,axis=2)/n
            contrast = packet['source'][:,:,0]-packet['source'][:,:,1]
            np.testing.assert_allclose(var,(n*raw_second-contrast**2)/(n-1),atol=1e-12)

    def test_joint_drift_geometry_allows_disjoint_valid_sets(self):
        rng = np.random.RandomState(3032)
        for _ in range(200):
            biases = rng.normal(size=(8,2))
            first = [0,1,2]
            second = [3,4,5]
            fixed_contrast = biases[first,0].mean()-biases[second,1].mean()
            full_contrast = np.mean(biases[:,0]-biases[:,1])
            bound = abs(fixed_contrast)+np.sqrt(5/3*np.var(biases,axis=0)).sum()
            self.assertLessEqual(abs(full_contrast),bound+1e-12)

    def test_joint_nuisance_and_interval_envelope(self):
        data = reuse_source_design_for_evaluation(fitted_stratum_data(
            dict(source_count=8,shared_scale=.3,pattern='weak_boundary',seed_offset=914),40))
        coefficients = conditional_design_coefficients(data,6,.001,.002)
        joint = nuisance_radius(data,coefficients,.002,'joint')
        separate = nuisance_radius(data,coefficients,.002,'separate')
        self.assertTrue(np.all(joint<=separate+1e-12))
        result = infer_joint_contrast(data,6,coefficients=coefficients)
        widths = np.diff(result['interval'],axis=1)[:,0]
        self.assertTrue(np.all(widths<=result['length_envelope']+1e-12))
        self.assertTrue(np.all((result['interval'][:,0]<=.1)&(.1<=result['interval'][:,1])))
        with self.assertRaises(ValueError):
            infer_joint_contrast(data,6,noise_failure=.04)

    def test_contrast_dispersion_requires_valid_intersection(self):
        data=reuse_source_design_for_evaluation(fitted_stratum_data(
            dict(source_count=8,shared_scale=.3,pattern='arm_cancel',seed_offset=113),10))
        with self.assertRaises(ValueError):infer_joint_contrast(data,4,dispersion_mode='contrast')
        result=infer_joint_contrast(data,6,dispersion_mode='contrast')
        np.testing.assert_array_equal(result['certificate_valid_minimum'],[4,4])
        self.assertEqual(result['dispersion_upper'].shape,(10,1))


if __name__=='__main__':
    unittest.main(verbosity=2)
