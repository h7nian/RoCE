"""Independent checks of the patient quadratic bound and probability regions."""
import itertools
import unittest
import numpy as np
from scipy.optimize import linprog
from concentrated_bounds import (exponential_dispersion_constants,
    exponential_dispersion_upper, pooled_mean_interval, simplex_linear_max,
    grouped_probability_max)


def patient_matrix(sizes):
    count = len(sizes)
    contrast = np.eye(count) / count - np.ones((count, count)) / count**2
    indices = np.repeat(np.arange(count), sizes)
    matrix = contrast[indices[:, None], indices] / (sizes[indices[:, None]] * sizes[indices])
    for site, size in enumerate(sizes):
        block = np.ix_(indices == site, indices == site)
        matrix[block] *= size / (size - 1)
    np.fill_diagonal(matrix, 0)
    return matrix, contrast, indices


class ConcentratedBoundChecks(unittest.TestCase):
    def test_block_spectrum_and_mgf_constants(self):
        rng = np.random.RandomState(5321)
        for _ in range(40):
            count = rng.randint(2, 8)
            sizes = rng.randint(2, 8, count)
            ranges = rng.uniform(.2, 3, count)
            matrix, contrast, indices = patient_matrix(sizes)
            proxy_sd = ranges[indices] / 2
            scaled = proxy_sd[:, None] * matrix * proxy_sd
            linear, quadratic, scale = exponential_dispersion_constants(ranges, sizes, (1, count))
            self.assertAlmostEqual(np.trace(scaled), 0)
            self.assertAlmostEqual(quadratic[0], 2 * np.sum(scaled**2), places=12)
            self.assertAlmostEqual(scale[0], -2 * np.linalg.eigvalsh(scaled)[0], places=12)
            means = rng.normal(size=count)
            dispersion = np.mean((means - means.mean())**2)
            vector = proxy_sd * (contrast @ means)[indices] / sizes[indices]
            self.assertLessEqual(4 * vector @ vector, linear[0] * dispersion + 1e-12)
            for fraction in [.01, .2, .75]:
                time = fraction / scale[0]
                shifted = np.eye(len(indices)) + 2 * time * scaled
                sign, logdet = np.linalg.slogdet(shifted)
                self.assertEqual(sign, 1)
                gaussian_mgf = -.5 * logdet + 2 * time**2 * vector @ np.linalg.solve(shifted, vector)
                bound = time**2 * (quadratic[0] + linear[0] * dispersion) / (2 * (1 - scale[0] * time))
                self.assertLessEqual(gaussian_mgf, bound + 1e-8)

    def test_exact_non_gaussian_enumeration(self):
        # Six independent patients with four possible outcomes: 4096 atoms.
        support = np.array([-.8, -.2, .4, 1.1])
        probabilities = np.array([[.1, .2, .3, .4], [.4, .3, .2, .1], [.2, .4, .1, .3]])
        configurations = np.array(list(itertools.product(range(4), repeat=6)))
        observations = support[configurations].reshape(-1, 3, 2)
        probability = np.prod(probabilities[np.repeat(np.arange(3), 2), configurations], axis=1)
        means = observations.mean(axis=2)
        variances = observations.var(axis=2, ddof=1) / 2
        centers = probabilities @ support
        dispersion = np.var(centers)
        estimate = np.var(means, axis=1) - 2 / 9 * variances.sum(axis=1)
        self.assertAlmostEqual(probability @ estimate, dispersion, places=12)
        ranges = np.full(3, np.ptp(support))
        linear, quadratic, scale = exponential_dispersion_constants(ranges, 2, (1, 3))
        matrix, contrast, indices = patient_matrix(np.full(3, 2))
        sd = ranges[indices] / 2
        scaled = sd[:, None] * matrix * sd
        vector = sd * (contrast @ centers)[indices] / 2
        for fraction in [.02, .1, .5]:
            time = fraction / scale[0]
            exact = np.log(probability @ np.exp(time * (dispersion - estimate)))
            shifted = np.eye(6) + 2 * time * scaled
            gaussian = -.5 * np.linalg.slogdet(shifted)[1] + 2 * time**2 * vector @ np.linalg.solve(shifted, vector)
            self.assertLessEqual(exact, gaussian + 1e-12)
        for alpha in [.01, .05, .2]:
            upper = exponential_dispersion_upper(means, variances, 2, ranges, alpha)
            self.assertLessEqual(probability @ (dispersion > upper), alpha)

    def test_pooled_mean_nonidentical_patient_variance(self):
        rng = np.random.RandomState(205)
        observations = rng.uniform(size=(5, 7, 10)) + np.arange(7)[None, :, None] / 4
        means = observations.mean(axis=2)
        mean_variances = observations.var(axis=2, ddof=1) / 10
        result = pooled_mean_interval(means, mean_variances, 10, 3., .05)
        pooled = observations.reshape(5, -1)
        radius = np.sqrt(2 * pooled.var(axis=1, ddof=1) / 70 * np.log(80)) + 7 * 6 * np.log(80) / (3 * 69)
        np.testing.assert_allclose(result, np.column_stack((pooled.mean(axis=1) - radius, pooled.mean(axis=1) + radius)))

    def test_probability_regions_match_linear_programming(self):
        rng = np.random.RandomState(629)
        for _ in range(100):
            probability = rng.dirichlet(np.ones(8)).reshape(4, 2)
            lower = probability * rng.uniform(size=(4, 2))
            upper = probability + rng.uniform(0, .2, (4, 2))
            group_lower = probability.sum(axis=1) * rng.uniform(size=4)
            group_upper = probability.sum(axis=1) + rng.uniform(0, .2, 4)
            coefficients = rng.normal(size=(4, 2))
            group_matrix = np.repeat(np.eye(4), 2, axis=1)
            fit = linprog(-coefficients.ravel(), A_ub=np.vstack((group_matrix, -group_matrix)),
                b_ub=np.r_[group_upper, -group_lower], A_eq=np.ones((1, 8)), b_eq=[1],
                bounds=list(zip(lower.ravel(), upper.ravel())), method='interior-point', options={'tol': 1e-11})
            self.assertTrue(fit.success)
            result = grouped_probability_max(coefficients, lower, upper, group_lower, group_upper)
            self.assertAlmostEqual(result, -fit.fun, places=8)
            simplex = simplex_linear_max(coefficients.ravel(), lower.ravel(), upper.ravel())
            fit = linprog(-coefficients.ravel(), A_eq=np.ones((1, 8)), b_eq=[1],
                bounds=list(zip(lower.ravel(), upper.ravel())), method='interior-point', options={'tol': 1e-11})
            self.assertAlmostEqual(simplex, -fit.fun, places=8)

    def test_random_variance_never_uses_exact_chi_square(self):
        from unittest.mock import patch
        from inference import dispersion_intervals
        means = np.array([[.1, .3, -.2], [.2, .4, .0]])
        with patch('inference.noncentrality_upper', side_effect=AssertionError('Exact branch used')):
            dispersion_intervals(means, .1, 2, .05, sample_variances=np.full_like(means, .1),
                                 sample_sizes=np.full(3, 50), variance_status='certified')
            dispersion_intervals(means, .1, 2, .05, variance_status='gaussian_plugin')
        with self.assertRaises(ValueError):
            dispersion_intervals(means, .1, 2, .05, variance_status='certified')

    def test_subset_geometry_and_fourier_drift(self):
        rng = np.random.RandomState(372)
        for _ in range(100):
            count, valid = 8, 6
            drift = rng.normal(size=count)
            dispersion = np.var(drift)
            for subset in itertools.combinations(range(count), valid):
                selected = drift[list(subset)]
                self.assertLessEqual(abs(drift.mean() - selected.mean()),
                                     np.sqrt((count-valid)/valid*dispersion) + 1e-12)
                frequency = rng.uniform(0, 10)
                exact = abs(np.sum(np.exp(1j*frequency*selected)-1)/count)
                bound = min(2*valid/count, valid/count*frequency*abs(selected.mean())
                            + frequency**2*np.sum(selected**2)/(2*count))
                self.assertLessEqual(exact, bound + 1e-12)

    def test_fitted_counterfactual_certificate_and_target_region(self):
        from fitted_strata import fitted_stratum_data
        from calibration_remainders import calibration_remainder_budgets, target_remainder_budget
        data = fitted_stratum_data(dict(source_count=8, shared_scale=.3,
                                       pattern='weak_boundary', seed_offset=173000), 40)
        certificate = calibration_remainder_budgets(data, 6)
        drift = data['nuisance_drift']
        for subset in itertools.combinations(range(8), 6):
            selected = drift[:, subset, :]
            self.assertTrue(np.all(abs(selected.mean(axis=1)) <= certificate['source_mean_budget'] + 1e-12))
            self.assertTrue(np.all(np.sum(selected**2, axis=1)/8 <= certificate['source_squared_budget'] + 1e-12))
        self.assertTrue(np.all(data['oracle_target_budget'] <= certificate['target_budget'] + 1e-12))
        target = target_remainder_budget(data['target_training_counts'], data['target_full_joint_counts'], .0025)
        self.assertTrue(np.all(data['oracle_target_budget'][:, 1:] <= target + 1e-12))


if __name__ == '__main__':
    unittest.main(verbosity=2)
