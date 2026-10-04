"""Independent geometry and calibration checks; run before freezing experiments."""
import unittest
import numpy as np
from scipy.optimize import minimize
from scipy.stats import ncx2
from inference import (TargetRegion, fourier_interval_set, fourier_summary,
                       count_interval_set, noncentrality_upper, intersect_sets)


class InferenceChecks(unittest.TestCase):
    def test_all_settings_respect_total_site_budget(self):
        from experiments import settings
        plan = settings()
        self.assertEqual(len(plan), 90)
        for row in plan:
            self.assertEqual(row["source_evaluation_size"] + row.get("pilot_size", 0), 1000)
            self.assertEqual(row["target_evaluation_size"], 1000)

    def test_minority_centers_are_not_uniquely_identified_by_sources(self):
        frequencies = np.linspace(.01, 100, 10000)
        for fraction in [.25, .5]:
            population = fraction + (1 - fraction) * np.exp(1j * frequencies)
            for candidate in [0., 1.]:
                self.assertTrue(np.all(abs(population - fraction * np.exp(1j * frequencies * candidate))
                                       <= 1 - fraction + 1e-12))

    def test_analytic_arcs_match_disk_constraints(self):
        rng = np.random.RandomState(191)
        grid = np.linspace(-4, 4, 20001)
        for _ in range(30):
            frequencies = np.sqrt(np.arange(1, 9)) / 2
            transform = rng.normal(size=8) + 1j * rng.normal(size=8)
            slack = rng.uniform(.1, 1, 8)
            intervals = fourier_interval_set(transform, frequencies, slack, .75, (-4, 4))
            direct = np.all(abs(transform[:, None] - .75 * np.exp(1j * frequencies[:, None] * grid))
                            <= .25 + slack[:, None], axis=0)
            included = np.zeros(len(grid), dtype=bool)
            for lower, upper in intervals:
                included |= (grid >= lower) & (grid <= upper)
            np.testing.assert_array_equal(included, direct)

    def test_periodic_branches_and_count_sweep(self):
        intervals = fourier_interval_set([1], [4], [.01], .75, (-5, 5))
        self.assertGreater(len(intervals), 3)
        rng = np.random.RandomState(92)
        grid = np.linspace(-4, 4, 20001)
        for count in [4, 8, 32]:
            estimates = rng.normal(size=count)
            variances = rng.uniform(.01, .3, count)
            accepted = count_interval_set(estimates, variances, count // 2, .05, (-4, 4))
            from scipy.stats import norm
            radius = norm.isf(.05 / (2 * count)) * np.sqrt(variances)
            direct = np.sum(abs(grid[:, None] - estimates) <= radius, axis=1) >= count // 2
            included = np.zeros(len(grid), dtype=bool)
            for lower, upper in accepted:
                included |= (lower <= grid) & (grid <= upper)
            np.testing.assert_array_equal(included, direct)
        self.assertEqual(intersect_sets([(0, 1)], [(2, 3)]), [])

    def test_projection_matches_independent_constrained_solver(self):
        rng = np.random.RandomState(741)
        for _ in range(60):
            draw = rng.normal(size=(3, 3))
            covariance = draw @ draw.T + np.eye(3) * .5
            center = rng.normal(size=3)
            region = TargetRegion(covariance)
            bounds = [tuple(center[index] + np.sort(rng.uniform(-1.5, 1.5, 2))) for index in [1, 2]]
            result = region.project(center, [[bounds[0]], [bounds[1]]])
            inverse = np.linalg.inv(covariance)
            start = center.copy()
            for index in [1, 2]:
                start[index] = np.clip(start[index], *bounds[index - 1])
            values = []
            for sign in [-1, 1]:
                fit = minimize(lambda value: sign * region.contrast @ value, start,
                    jac=lambda value: sign * region.contrast,
                    method="SLSQP", bounds=[(None, None)] + bounds,
                    constraints=[{"type": "ineq", "fun": lambda value:
                        region.critical_squared - (value - center) @ inverse @ (value - center),
                        "jac": lambda value: -2 * inverse @ (value - center)}],
                    options={"ftol": 1e-11, "maxiter": 400})
                if not fit.success and fit.fun is None:
                    self.fail("Independent optimizer did not return a solution")
                self.assertLess(abs(region.critical_squared - (fit.x - center) @ inverse @ (fit.x - center)), 1e-6)
                values.append(float(region.contrast @ fit.x))
            np.testing.assert_allclose(result, [min(values), max(values)], atol=2e-6, rtol=1e-7)

    def test_singular_prediction_and_empty_projection(self):
        covariance = np.array([[0, 0, 0], [0, 1, -.3], [0, -.3, 1.]])
        region = TargetRegion(covariance)
        center = np.array([.1, 0, 0])
        result = region.project(center, [[(-100, 100)], [(-100, 100)]])
        radius = np.sqrt(region.critical_squared * region.effect_variance)
        np.testing.assert_allclose(result, [.1 - radius, .1 + radius], atol=1e-9)
        self.assertIsNone(region.project(center, [[(100, 101)], [(0, 1)]]))
        fixed = region.project(center, [[(.2, .2)], [(-.1, -.1)]])
        np.testing.assert_allclose(fixed, [.4, .4], atol=1e-9)

    def test_noncentrality_conservative_endpoint(self):
        for degrees in [3, 31, 511]:
            statistics = np.array([.1, degrees, degrees + 10 * np.sqrt(degrees)])
            upper = noncentrality_upper(statistics, degrees, .01)
            active = upper > 0
            self.assertTrue(np.all(ncx2.cdf(statistics[active], degrees, upper[active]) <= .01 + 1e-10))

    def test_gaussian_concentration_and_variance_error(self):
        rng = np.random.RandomState(635)
        count, repeats = 128, 4000
        variances = np.linspace(.5, 1.5, count)
        means = np.zeros(count)
        means[::4] = .5
        data = means + rng.normal(size=(repeats, count)) * np.sqrt(variances)
        for calibration in ["hoeffding", "bernstein"]:
            transform, frequencies, slack = fourier_summary(data, variances, .05, calibration=calibration)
            population = np.mean(np.exp(1j * frequencies[0, :, None] * means), axis=1)
            failure = np.mean(np.any(abs(transform - population) > slack, axis=1))
            self.assertLess(failure, .05)
        lower, upper = .8 * variances, 1.2 * variances
        transform, frequencies, slack = fourier_summary(data, lower, .05, variance_upper=np.broadcast_to(upper, data.shape))
        truth_failure = np.mean(np.any(abs(transform - .75) > .25 + slack, axis=1))
        self.assertLess(truth_failure, .05)

    def test_lindeberg_remainder_contains_exact_bernoulli_error(self):
        for probability in [.05, .3, .7]:
            for size in [20, 1000]:
                variance = probability * (1 - probability) / size
                frequencies = np.sqrt(np.arange(1, 9)) / (2 * np.sqrt(variance))
                exact = ((1 - probability) + probability * np.exp(1j * frequencies / size))**size
                corrected = exact * np.exp(frequencies**2 * variance / 2)
                normal_limit = np.exp(1j * frequencies * probability)
                patient_variance = variance * size
                bound = np.exp(frequencies**2 * variance / 2) * frequencies**3 / (6 * size**2) * (
                    patient_variance + 2 * np.sqrt(2 / np.pi) * patient_variance**1.5)
                self.assertTrue(np.all(abs(corrected - normal_limit) <= bound + 1e-10))


if __name__ == "__main__":
    unittest.main(verbosity=2)
