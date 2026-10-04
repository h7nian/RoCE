import unittest
import numpy as np
from scipy.optimize import linprog
from experiments import gaussian_data, settings
from inference import TargetRegion
from repairs import (DIRECTIONS, DirectionalTargetRegion, directional_radii,
                     protected_projection, infer_repair, variance_certificate)


class RepairChecks(unittest.TestCase):
    def test_invalid_inference_contracts_are_rejected(self):
        data=gaussian_data(settings()[0],2)
        for kwargs in [dict(source_bias_budget=-.1),dict(target_bias_budget=np.nan),
                       dict(valid_minimum=2.5),dict(variance_failure=-.01),dict(target_alpha=.06)]:
            with self.assertRaises(ValueError):infer_repair(data,**kwargs)

    def test_fitted_stratum_identities_and_budgets(self):
        from fitted_strata import fitted_stratum_data
        setting=dict(source_count=8,shared_scale=.3,pattern="weak_boundary",seed_offset=9000000)
        data=fitted_stratum_data(setting,100)
        np.testing.assert_allclose(data["target_truth"]@DIRECTIONS[3],.1,atol=1e-12)
        valid=data["structural_shift"]==0
        deviations=data["source_expected"]-data["residual_truth"][:,None,:]
        np.testing.assert_allclose(deviations[:,valid],data["nuisance_drift"][:,valid],atol=1e-12)
        self.assertTrue(np.all(data["oracle_source_budget"]<=data["feasible_source_budget"]+1e-12))
        self.assertTrue(np.all(data["oracle_target_budget"]<=data["feasible_target_budget"]+1e-12))
        self.assertTrue(np.all(data["oracle_source_budget"]<=data["shared_source_budget"]+1e-12))
        self.assertTrue(np.all(data["oracle_target_budget"]<=data["shared_target_budget"]+1e-12))
        self.assertEqual(data["source_training_size"]+data["pilot_size"]+data["source_evaluation_size"],1000)
        self.assertEqual(data["target_training_size"]+data["target_sample_size"],1000)

    def test_directional_projection_matches_linear_programming(self):
        rng = np.random.RandomState(172)
        for _ in range(40):
            center = rng.normal(size=3)
            radii = rng.uniform(.2, 2, 4)
            region = DirectionalTargetRegion(radii)
            bounds = [tuple(sorted(rng.uniform(-1, 1, 2))) for _ in range(2)]
            result = region.project(center, [[bounds[0]], [bounds[1]]])
            matrix = np.vstack((DIRECTIONS, -DIRECTIONS))
            upper = np.r_[DIRECTIONS @ center + radii, -DIRECTIONS @ center + radii]
            values = []
            for sign in [-1, 1]:
                fit = linprog(sign * DIRECTIONS[3], A_ub=matrix, b_ub=upper,
                              bounds=[(None, None)] + bounds, method="revised simplex")
                if fit.success:
                    values.append(float(DIRECTIONS[3] @ fit.x))
            if values:
                np.testing.assert_allclose(result, [min(values), max(values)], atol=1e-8)
            else:
                self.assertIsNone(result)

    def test_fallback_length_cap_and_nonempty_parity(self):
        region = TargetRegion(np.diag([0., 1., 1.]))
        for center in [np.array([.1, 0, 0]), np.array([.1, 10, 10])]:
            source = [[(-.1, .1)], [(-.1, .1)]]
            moments = [(-.2, .2), (-.2, .2)]
            reference = (float(DIRECTIONS[3] @ center)-3, float(DIRECTIONS[3] @ center)+3)
            result = protected_projection(center, region, source, moments, 0, reference, "shorter")
            self.assertLessEqual(result[0][1]-result[0][0], .8 + 1e-8)
            self.assertTrue(result[0][0] <= result[1] <= result[0][1])
            if center[1] == 0:
                np.testing.assert_allclose(result[0], region.project(center, source))
            else:
                self.assertEqual(result[2], 2)

    def test_bias_budgets_expand_region(self):
        covariance = np.array([[.1, .02, -.01], [.02, .4, -.05], [-.01, -.05, .3]])
        ranges = np.array([1., 2., 2., 3.])
        plain = directional_radii(covariance, ranges, 100, .025, "empirical")
        budget = np.array([0, .1, .2])
        protected = directional_radii(covariance, ranges, 100, .025, "empirical", budget)
        np.testing.assert_allclose(protected-plain, [0, .1, .2, .3])

    def test_variance_certificate_in_exact_binary_law(self):
        rng = np.random.RandomState(991)
        repeats, count, pilot, evaluation = 3000, 8, 250, 750
        probability = np.linspace(.1, .8, count)
        successes = rng.binomial(pilot, probability[None, :, None], (repeats, count, 2))
        estimate = successes*(pilot-successes)/(pilot*(pilot-1)*evaluation)
        lower, upper = variance_certificate(estimate, np.ones((count,2)), pilot, evaluation, .01)
        truth = probability*(1-probability)/evaluation
        failed = np.any((lower > truth[None, :, None]) | (upper < truth[None, :, None]), axis=(1,2))
        self.assertLessEqual(failed.mean(), .01)

    def test_shorter_and_target_rules_only_differ_on_empty_sets(self):
        setting = next(row for row in settings() if row["source_count"]==128 and row["shared_scale"]==0)
        setting["seed_offset"] = 5100000
        data = gaussian_data(setting, 100)
        first = infer_repair(data, fallback_rule="target")
        second = infer_repair(data, fallback_rule="shorter")
        nonempty = first["fallback"] == 0
        np.testing.assert_allclose(first["interval"][nonempty], second["interval"][nonempty], atol=1e-10)
        self.assertTrue(np.all(np.diff(second["interval"], axis=1)[:,0] <= second["length_envelope"] + 1e-8))


if __name__ == "__main__":
    unittest.main(verbosity=2)
