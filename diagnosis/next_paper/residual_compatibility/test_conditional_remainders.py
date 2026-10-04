import itertools
import unittest
import numpy as np
from scipy.stats import binom
from conditional_remainders import outcome_linear_radius, conditional_remainder_budgets, conditional_target_budget
from fitted_strata import fitted_stratum_data


class ConditionalRemainderChecks(unittest.TestCase):
    def test_scalar_outcome_bound_by_exact_enumeration(self):
        counts = np.array([2, 3, 4])
        coefficients = np.array([.3, -.2, .5])
        probabilities = np.array([.1, .5, .8])
        outcomes = np.array(list(itertools.product(range(3), range(4), range(5))))
        law = np.prod(binom.pmf(outcomes, counts, probabilities), axis=1)
        drift = (probabilities - (outcomes + .5) / (counts + 1)) @ coefficients
        for failure in [.01, .05, .2]:
            radius = outcome_linear_radius(abs(coefficients), counts, failure)
            self.assertLessEqual(law @ (abs(drift) > radius), failure)

    def test_fitted_fixed_subset_and_target_budgets(self):
        data = fitted_stratum_data(dict(source_count=8, shared_scale=.3,
                                       pattern='weak_boundary', seed_offset=103000), 100)
        budgets = conditional_remainder_budgets(data, 6)
        for arm in [0, 1]:
            valid = np.flatnonzero(data['structural_shift'][:, arm] == 0)[:6]
            drift = abs(data['nuisance_drift'][:, valid, arm].mean(axis=1))
            self.assertTrue(np.all(drift <= budgets['source_mean_budget'][:, arm] + 1e-12))
        self.assertTrue(np.all(data['oracle_target_budget'] <= budgets['target_budget'] + 1e-12))
        target = conditional_target_budget(data['target_training_counts'], data['target_full_joint_counts'], .0025)
        self.assertTrue(np.all(data['oracle_target_budget'][:, 1:] <= target + 1e-12))

    def test_reused_blocks_match_pooled_patient_moments(self):
        from conditional_remainders import reuse_source_design_for_evaluation
        rng = np.random.RandomState(1032)
        pilot = rng.uniform(-1, 1, (5, 4, 7, 2))
        evaluation = rng.uniform(-1, 1, (5, 4, 11, 2))
        packet = dict(pilot_size=7, source_evaluation_size=11,
            pilot_source_mean=pilot.mean(axis=2), source=evaluation.mean(axis=2),
            pilot_variance=pilot.var(axis=2, ddof=1)/11,
            sample_variance=evaluation.var(axis=2, ddof=1)/11,
            pilot_joint_frequencies=np.full((5, 4, 4, 2), .125),
            evaluation_joint_frequencies=np.full((5, 4, 4, 2), .125), variance=np.ones((5,4,2))/11)
        reused = reuse_source_design_for_evaluation(packet)
        pooled = np.concatenate((pilot, evaluation), axis=2)
        np.testing.assert_allclose(reused['source'], pooled.mean(axis=2))
        np.testing.assert_allclose(reused['sample_variance'], pooled.var(axis=2, ddof=1)/18)
        self.assertEqual(reused['source_design_sample_size'], 18)
        self.assertEqual(reused['pilot_size'], 0)

    def test_reuse_rejects_unproved_fourier_calibration(self):
        from repairs import infer_repair
        with self.assertRaises(ValueError):
            infer_repair(dict(source_evaluation_reuses_design=True), source_mode='hybrid')

    def test_unsplit_target_handles_empty_arms(self):
        from conditional_remainders import full_target_stratified_interval
        states = np.array(list(itertools.product([0, 1], repeat=2)))
        configurations = np.array(list(itertools.product(range(4), repeat=4)))
        halves = []
        for block in [configurations[:, :2], configurations[:, 2:]]:
            counts = np.zeros((len(block), 1, 2, 2), dtype=int)
            for row, patients in enumerate(block):
                for index in patients:
                    arm, outcome = states[index]
                    counts[row, 0, arm, outcome] += 1
            halves.append(counts)
        result = full_target_stratified_interval(dict(target_training_counts=halves[0], target_evaluation_counts=halves[1]))
        self.assertTrue(np.all(np.isfinite(result['interval'])))
        self.assertTrue(np.all(result['interval'][:, 0] <= result['interval'][:, 1]))
        probability = np.array([.15, .35, .4, .1])
        law = np.prod(probability[configurations], axis=1)
        covered = (result['interval'][:, 0] <= .5) & (.5 <= result['interval'][:, 1])
        self.assertGreaterEqual(law @ covered, .95 - 1e-12)

    def test_outcome_dependent_weight_contract_is_rejected(self):
        with self.assertRaises(ValueError):
            conditional_remainder_budgets(dict(nuisance_structure='outcome_calibrated_weights'), 2)


if __name__ == '__main__':
    unittest.main(verbosity=2)
