import itertools
import unittest

import numpy as np
from scipy.special import logsumexp

from continuous_balance import (UnsupportedDesign, linear_balance_design, evaluate_linear_balance,
    leaveout_dispersion_constants, leaveout_dispersion_upper, continuous_balance_interval,
    contrast_summaries)


def stack_summaries(values):
    return {name: np.array([entry[name] for entry in values]) for name in values[0]}


def dense_quadratic(designs):
    size = len(designs)
    lengths = [len(d['weights']) for d in designs]
    offsets = np.cumsum([0] + lengths)
    weights = np.zeros((offsets[-1], size))
    correction = np.zeros((offsets[-1], offsets[-1]))
    centering = np.eye(size) / size - np.ones((size, size)) / size**2
    for j, design in enumerate(designs):
        block = slice(offsets[j], offsets[j + 1])
        weights[block, j] = design['weights']
        hat = design['basis'].dot(design['basis'].T)
        diagonal = np.diag(design['diagonal'])
        c = diagonal - (diagonal.dot(hat) + hat.dot(diagonal)) / 2
        correction[block, block] = centering[j, j] * c
    return weights.dot(centering).dot(weights.T) - correction, weights, centering


class ContinuousBalanceTests(unittest.TestCase):
    def test_balance_and_outcome_predictor_cancellation(self):
        rng = np.random.RandomState(18001)
        features = np.column_stack((np.ones(30), rng.uniform(-1, 1, (30, 4))))
        target = np.array([1., .1, -.1, .2, 0.])
        design = linear_balance_design(features, target)
        y = rng.binomial(1, .4, len(features))
        beta = rng.normal(size=5)
        weighted = design['weights'].dot(y)
        augmented = target.dot(beta) + design['weights'].dot(y - features.dot(beta))
        self.assertAlmostEqual(weighted, augmented, places=12)
        np.testing.assert_allclose(features.T.dot(design['weights']), target, atol=1e-12)

    def test_dense_geometry_and_linear_remainder_bound(self):
        rng = np.random.RandomState(19001)
        for _ in range(40):
            designs, values, means = [], [], []
            for size in [8, 11, 13]:
                x = np.column_stack((np.ones(size), rng.uniform(-1, 1, (size, 2))))
                d = linear_balance_design(x, np.array([1., .1, -.15]))
                mu = .35 + x[:, 1:].dot(rng.uniform(-.12, .12, 2))
                designs.append(d); means.append(mu)
                values.append(evaluate_linear_balance(d, rng.binomial(1, mu)))
            summary = stack_summaries(values)
            q, weights, a = dense_quadratic(designs)
            mu = np.concatenate(means) - .5
            theta = weights.T.dot(mu)
            energy = theta.dot(a).dot(theta)
            linear, constant, scale = leaveout_dispersion_constants(summary)
            correction = ((len(designs) - 1) / len(designs)**2)**2 * sum(
                d['diagonal_square_sum'] for d in designs) / 8
            np.testing.assert_allclose(np.diag(q), 0, atol=2e-16)
            self.assertAlmostEqual(constant, correction + np.sum(q*q)/8, places=12)
            self.assertGreaterEqual(np.linalg.eigvalsh(q/4).min(), -scale/2 - 1e-12)
            self.assertLessEqual(np.linalg.norm(q.dot(mu))**2, linear*energy + correction + 1e-12)
            self.assertAlmostEqual(mu.dot(q).dot(mu), energy, places=12)

    def test_non_gaussian_exact_enumeration(self):
        x = np.column_stack((np.ones(4), [-.9, -.3, .2, .8]))
        designs = [linear_balance_design(x, np.array([1., .1])) for _ in range(2)]
        probabilities = np.concatenate((.3 + .1*x[:, 1], .5 - .15*x[:, 1]))
        q, weights, a = dense_quadratic(designs)
        centers = weights.T.dot(probabilities - .5)
        truth = centers.dot(a).dot(centers)
        draws = np.array(list(itertools.product([0., 1.], repeat=8)))
        mass = np.prod(np.where(draws == 1, probabilities, 1-probabilities), axis=1)
        statistics, variances, uppers = [], [], []
        for draw in draws:
            values = [evaluate_linear_balance(designs[j], draw[j*4:(j+1)*4]) for j in range(2)]
            summary = stack_summaries(values)
            statistic = np.var(summary['estimate']) - .25*np.sum(summary['variance'])
            statistics.append(statistic); variances.append(summary['variance'])
            uppers.append(leaveout_dispersion_upper(summary, .05))
            self.assertAlmostEqual(statistic, (draw-.5).dot(q).dot(draw-.5), places=12)
        statistics = np.array(statistics)
        linear, constant, scale = leaveout_dispersion_constants(summary)
        self.assertAlmostEqual(mass.dot(statistics), truth, places=12)
        self.assertTrue(np.any(np.asarray(variances) < 0))
        for j in range(2):
            analytic = np.sum(designs[j]['weights']**2 * probabilities[j*4:(j+1)*4]
                              * (1-probabilities[j*4:(j+1)*4]))
            self.assertAlmostEqual(mass.dot(np.array(variances)[:, j]), analytic, places=12)
        for t in [.1, 1., .2/scale]:
            log_mgf = logsumexp(np.log(mass) - t*(statistics-truth))
            bound = t*t*(linear*truth+constant)/(2*(1-scale*t))
            self.assertLessEqual(log_mgf, bound+1e-12)
        self.assertLessEqual(mass[np.array(uppers) < truth].sum(), .05+1e-12)

    def test_fast_variance_matches_actual_leave_one_out_refits(self):
        rng = np.random.RandomState(21001)
        x = np.column_stack((np.ones(15), rng.uniform(-1, 1, (15, 3))))
        design = linear_balance_design(x, np.array([1., .1, 0., -.1]))
        y = rng.binomial(1, .4, len(x))
        direct = 0.
        for index in range(len(x)):
            keep = np.arange(len(x)) != index
            coefficient = np.linalg.lstsq(x[keep], y[keep]-.5, rcond=None)[0]
            direct += design['weights'][index]**2 * (y[index]-.5) * (y[index]-.5-x[index].dot(coefficient))
        self.assertAlmostEqual(evaluate_linear_balance(design, y)['variance'], direct, places=12)

    def test_target_rank_failure_keeps_the_repetition(self):
        from run_continuous_balance import run_repeat
        setting = dict(n_per_site=5, dimension=4, source_count=4,
                       heterogeneity=.15, source_shift=.4, pattern='all_valid',
                       seed_offset=23001, cell_id=1)
        result = run_repeat(setting, 0)
        self.assertTrue(result['diagnostics']['target_failed'])
        self.assertEqual(len(result['methods']), 7)
        for value in result['methods'].values():
            np.testing.assert_array_equal(value['interval'], [-1., 1.])

    def test_guards_and_summary_contrast(self):
        with self.assertRaises(UnsupportedDesign):
            linear_balance_design(np.ones((3, 3)), np.ones(3))
        rng = np.random.RandomState(22001)
        values = []
        for _ in range(4):
            arms = []
            for arm in range(2):
                x = np.column_stack((np.ones(50), rng.uniform(-1, 1, (50, 2))))
                design = linear_balance_design(x, np.array([1., .1, 0.]))
                arms.append(evaluate_linear_balance(design, rng.binomial(1, .4+.1*arm, 50)))
            values.append(stack_summaries(arms))
        source = {name: np.stack([v[name] for v in values]) for name in values[0]}
        target = values[0]
        contrast = contrast_summaries(source)
        np.testing.assert_allclose(contrast['variance'], source['variance'].sum(axis=1))
        for mode in ['arms', 'contrast', 'adaptive']:
            result = continuous_balance_interval(source, target, [3, 3], 100, mode)
            self.assertLessEqual(result['interval'][0], result['target_projection'])
            self.assertGreaterEqual(result['interval'][1], result['target_projection'])
        with self.assertRaisesRegex(ValueError, 'intersection'):
            continuous_balance_interval(source, target, [2, 2], 100, 'contrast')


if __name__ == '__main__':
    unittest.main()
