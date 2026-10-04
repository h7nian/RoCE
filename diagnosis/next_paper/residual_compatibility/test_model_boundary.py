import itertools
import unittest
import numpy as np
from scipy.stats import binom

from continuous_balance import linear_balance_design, evaluate_linear_balance
from run_model_boundary import contamination_overlap, saturated_geometry, saturated_summaries


class ModelBoundaryTests(unittest.TestCase):
    def test_exact_contamination_mixture_and_count_conditioning(self):
        for size in [4, 125, 500, 2000]:
            audit = contamination_overlap(size, size, 2048)
            self.assertLess(audit['mixture_error'], 1e-14)
            for name in ['common_mass', 'q0_mass', 'q1_mass']:
                self.assertAlmostEqual(audit[name], 1., places=12)
            self.assertGreaterEqual(audit['minimum_contamination_mass'], -1e-14)
            self.assertGreater(audit['scaled_lower_bound'], .07)
        # The finite-source-count penalty is essential at small K.
        self.assertLess(contamination_overlap(500, 500, 8)['interval_length_lower_bound'],
                        contamination_overlap(500, 500, 2048)['interval_length_lower_bound'])
        draws = np.array(list(itertools.product([0, 1], repeat=4)))
        delta = .1/2
        masses = [np.prod(np.where(draws, theta, 1-theta), axis=1) for theta in [.5, .5+delta]]
        vector_tv = np.abs(masses[0]-masses[1]).sum()/2
        count_tv = np.abs(binom.pmf(np.arange(5), 4, .5)-binom.pmf(np.arange(5), 4, .5+delta)).sum()/2
        self.assertAlmostEqual(vector_tv, count_tv, places=13)

    def test_saturated_group_formula_matches_patient_leaveout(self):
        values, counts, diagonal, constants = saturated_geometry(40)
        x = np.repeat(values, counts)
        matrix = np.column_stack((np.ones(len(x)), x, x*x))
        design = linear_balance_design(matrix, matrix.mean(axis=0))
        rng = np.random.RandomState(32401)
        for _ in range(12):
            successes = rng.binomial(counts, [.7, .3, .8])
            outcomes = np.concatenate([np.r_[np.ones(s), np.zeros(n-s)] for s, n in zip(successes, counts)])
            raw = evaluate_linear_balance(design, outcomes)
            grouped = saturated_summaries(successes, counts, diagonal, constants)
            for name in ['estimate', 'variance', 'weight_square_sum', 'diagonal_square_sum']:
                self.assertAlmostEqual(raw[name], grouped[name], places=13)


if __name__ == '__main__':
    unittest.main()
