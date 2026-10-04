import itertools
import unittest
import numpy as np
from scipy.special import logsumexp

from continuous_balance import (linear_balance_design, evaluate_linear_balance,
    leaveout_dispersion_constants, leaveout_dispersion_upper)
from population_protection import (target_coefficient_certificate, composition_radius,
    projection_error_budget, projection_geometry)
from test_continuous_balance import dense_quadratic, stack_summaries


class PopulationProtectionTests(unittest.TestCase):
    def test_joint_region_bounds_range_and_empirical_variance(self):
        rng = np.random.RandomState(24001)
        features = rng.uniform(-1, 1, (300, 4))
        treatment = np.tile([0, 1], 150)
        outcomes = rng.binomial(1, .4+.1*treatment+.04*features[:, 0])
        certificate = target_coefficient_certificate(features, treatment, outcomes)
        factor = np.linalg.cholesky(certificate['proxy'])
        directions = rng.normal(size=(1000, 4))
        directions /= np.linalg.norm(directions, axis=1)[:, None]
        coefficients = certificate['coefficient'] + certificate['ellipsoid_radius']*directions.dot(factor.T)
        # Only model-admissible CATE vectors are relevant to the universal cap 2.
        for coefficient in coefficients:
            if 2*np.abs(coefficient).sum() <= 2:
                self.assertLessEqual(2*np.abs(coefficient).sum(), certificate['range_upper']+1e-12)
                variance = np.var(features.dot(coefficient), ddof=1)
                self.assertLessEqual(variance, certificate['sample_variance_upper']+1e-12)
        self.assertGreater(composition_radius(certificate, 300, .01, 'bernstein'), 0)
        with self.assertRaises(ValueError):
            target_coefficient_certificate(features, treatment, outcomes, support_lower=0.)

    def test_binomial_summary_matches_patient_fit(self):
        from run_projection_stress import fixed_geometry, summaries
        from scipy.stats import binom
        from run_design_stability import chernoff_design_bound
        values, counts, diagonal, constants = fixed_geometry()
        x = np.repeat(values, counts)
        design = linear_balance_design(np.column_stack((np.ones(len(x)), x)), np.array([1., 0.]))
        rng = np.random.RandomState(24201)
        for _ in range(10):
            successes = rng.binomial(counts, [.3, .7, .3])
            y = np.concatenate([np.r_[np.ones(s), np.zeros(n-s)] for s,n in zip(successes,counts)])
            direct = evaluate_linear_balance(design, y)
            grouped = summaries(successes, values, counts, diagonal, constants)
            self.assertAlmostEqual(direct['estimate'], grouped['estimate'], places=12)
            self.assertAlmostEqual(direct['variance'], grouped['variance'], places=12)
        # Exact scalar Gram failure for Bernoulli(.5) treatment counts.
        for count in [1,4,16]:
            exact = 1-(1-2*binom.cdf(25,100,.5))**count
            self.assertLessEqual(exact, chernoff_design_bound(count,0,100,.5)['any_source_gram_failure'])

    def test_nonlinear_dgp_support_before_simulation(self):
        from run_target_protection import plan, validate_nonlinear_support
        settings = plan('misspecification')
        for setting in settings:
            validate_nonlinear_support(setting)
        invalid = dict(settings[-1], nonlinearity=.2)
        with self.assertRaises((ArithmeticError, ValueError)):
            validate_nonlinear_support(invalid)

    def test_exact_binary_coefficient_region(self):
        x = np.tile(np.array([-1., -.3, .3, 1.]), 2)[:, None]
        a = np.repeat([1, 0], 4)
        probability = .3+.2*a+(.05+.05*a)*x[:, 0]
        draws = np.array(list(itertools.product([0., 1.], repeat=8)))
        mass = np.prod(np.where(draws == 1, probability, 1-probability), axis=1)
        errors, failed = [], []
        for y in draws:
            c = target_coefficient_certificate(x, a, y, alpha=.1)
            error = c['coefficient']-np.array([.05])
            errors.append(error)
            failed.append(error.dot(np.linalg.solve(c['proxy'], error)) > c['ellipsoid_radius']**2)
        self.assertLessEqual(mass[np.array(failed)].sum(), .1)
        for direction in [-20., -3., 3., 20.]:
            log_mgf = logsumexp(np.log(mass)+np.array(errors)[:, 0]*direction)
            self.assertLessEqual(log_mgf, .5*direction**2*c['proxy'][0, 0]+1e-12)

    def test_misspecified_variance_identity_and_mgf(self):
        x = np.column_stack((np.ones(4), [-.9, -.3, .2, .8]))
        designs = [linear_balance_design(x, np.array([1., .1])) for _ in range(2)]
        probabilities = np.array([.2, .8, .3, .7, .6, .2, .7, .4])
        q, weights, a = dense_quadratic(designs)
        theta = weights.T.dot(probabilities-.5)
        energy = theta.dot(a).dot(theta)
        rho, bias = [], []
        summaries = []
        for j, design in enumerate(designs):
            m = probabilities[j*4:(j+1)*4]-.5
            residual = m-design['basis'].dot(design['basis'].T.dot(m))
            rho.append(np.linalg.norm(residual))
            bias.append(np.sum(design['diagonal']*m*residual))
            summary = evaluate_linear_balance(design, np.ones(4)*.5)
            summary.update(projection_geometry(design)); summaries.append(summary)
            dense = (np.eye(4)-design['basis'].dot(design['basis'].T)).dot(np.diag(design['diagonal'])).dot(design['basis'])
            self.assertAlmostEqual(np.linalg.norm(dense, 2), summary['projection_commutator_norm'], places=12)
        summary = stack_summaries(summaries)
        budget = projection_error_budget(summary, np.array(rho))
        self.assertLessEqual(abs(sum(bias)), budget['variance_bias_sum']+1e-12)
        draws = np.array(list(itertools.product([0., 1.], repeat=8)))
        mass = np.prod(np.where(draws == 1, probabilities, 1-probabilities), axis=1)
        statistic = np.einsum('bi,ij,bj->b', draws-.5, q, draws-.5)
        expectation = energy-.25*sum(bias)
        self.assertAlmostEqual(mass.dot(statistic), expectation, places=12)
        linear, constant, scale = leaveout_dispersion_constants(summary)
        gsum = summary['diagonal_square_sum'].sum()
        projection = budget['weighted_projection_error']
        constant += .25**2/2*(projection+np.sqrt(projection*gsum))
        for t in [.1, 1., .2/scale]:
            log_mgf = logsumexp(np.log(mass)-t*(statistic-expectation))
            self.assertLessEqual(log_mgf, t*t*(linear*energy+constant)/(2*(1-scale*t))+1e-12)
        uppers = []
        for y in draws:
            values = stack_summaries([evaluate_linear_balance(designs[j], y[j*4:(j+1)*4]) for j in range(2)])
            uppers.append(leaveout_dispersion_upper(values, .05, **budget))
        self.assertLessEqual(mass[np.array(uppers)<energy].sum(), .05+1e-12)
        self.assertEqual(projection_error_budget(summary, np.array(rho), 0)['variance_bias_sum'], 0.)


if __name__ == '__main__':
    unittest.main()
