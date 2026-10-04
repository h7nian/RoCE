"""Independent finite-distribution and paired checks for training averaging."""
import itertools
import unittest
import numpy as np
from scipy.special import logsumexp
from scipy.stats import binom
from training_average import binomial_ratio_mgf, ratio_average_interval, training_average_budgets
from fitted_strata import fitted_stratum_data


class TrainingAverageChecks(unittest.TestCase):
    def test_uniform_probability_envelope_between_grid_points(self):
        envelope = binomial_ratio_mgf(20, 4, .03, grid_size=64, tilts=np.array([.01, .2, .6, 1.5]))
        successes = np.arange(21)
        for probability in np.r_[np.geomspace(.03, .9, 127), 1.]:
            log_probability = binom.logpmf(successes, 20, probability)
            ratio = 24 * probability / (successes + .5)
            for index, tilt in enumerate(envelope['tilts']):
                upper = logsumexp(log_probability + tilt * (ratio - 1))
                lower = logsumexp(log_probability - tilt * (ratio - 1))
                self.assertLessEqual(upper, envelope['log_upper'][index] + 1e-12)
                self.assertLessEqual(lower, envelope['log_lower'][index] + 1e-12)

    def test_exact_heterogeneous_binomial_pair_coverage(self):
        envelope = binomial_ratio_mgf(8, 2, .1, grid_size=128)
        configurations = np.array(list(itertools.product(range(9), repeat=2)))
        for probabilities in [[.1, .3], [.2, .7], [.5, .9]]:
            law = np.prod(binom.pmf(configurations, 8, probabilities), axis=1)
            average = np.mean(10 * np.array(probabilities) / (configurations + .5), axis=1)
            for failure in [.05, .2]:
                lower, upper = ratio_average_interval(envelope, 2, failure)
                self.assertLessEqual(law @ ((average < lower) | (average > upper)), failure + 1e-12)

    def test_refined_probability_grid_is_conservative_and_tighter(self):
        coarse = binomial_ratio_mgf(25, 4, .04, grid_size=32)
        fine = binomial_ratio_mgf(25, 4, .04, grid_size=64)
        self.assertTrue(np.all(fine['log_upper'] <= coarse['log_upper'] + 1e-12))
        self.assertTrue(np.all(fine['log_lower'] <= coarse['log_lower'] + 1e-12))
        first = ratio_average_interval(coarse, 96, .003/8)
        second = ratio_average_interval(fine, 96, .003/8)
        self.assertGreaterEqual(second[0] + 1e-12, first[0])
        self.assertLessEqual(second[1], first[1] + 1e-12)

    def test_fitted_score_fixed_valid_subset_budget(self):
        envelope = binomial_ratio_mgf(500, 4, .04, grid_size=128, tilts=np.geomspace(.01, 1, 30))
        data = fitted_stratum_data(dict(source_count=32, shared_scale=.3, pattern='weak_boundary', seed_offset=190321), 60)
        budget = training_average_budgets(data, 24, envelope)
        for arm in [0, 1]:
            valid = np.flatnonzero(data['structural_shift'][:, arm] == 0)[:24]
            actual = abs(data['nuisance_drift'][:, valid, arm].mean(axis=1))
            self.assertTrue(np.all(actual <= budget['source_mean_budget'][:, arm] + 1e-12))
        self.assertTrue(np.all(data['oracle_target_budget'] <= budget['target_budget'] + 1e-12))


if __name__ == '__main__':
    unittest.main(verbosity=2)
