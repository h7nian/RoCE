import itertools
import unittest
import numpy as np

from continuous_balance import conditional_linear_bands, evaluate_linear_balance
from design_selection import (select_source_designs, evaluate_retained_sources,
    conditional_selected_bands, retained_projection_budgets)
from population_protection import projection_geometry, projection_error_budget


class DesignSelectionTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.RandomState(32201)
        self.designs = [(np.column_stack((np.ones(80), rng.uniform(-1, 1, (80, 2)))),
                         np.tile([1, 0], 40)) for _ in range(8)]
        self.outcomes = [rng.binomial(1, .45, 80) for _ in self.designs]
        self.target_mean = np.array([1., .03, -.02])

    def selection(self, failed=(), valid=6, policy='remove'):
        designs = [(matrix.copy(), treatment.copy()) for matrix, treatment in self.designs]
        for index in failed:
            designs[index][0][:, 2] = designs[index][0][:, 1]
        return select_source_designs(designs, self.target_mean, valid, policy)

    def target(self, selection):
        values = [evaluate_linear_balance(arm['design'], self.outcomes[0][arm['rows']])
                  for arm in selection['designs'][0]]
        return {key: np.array([v[key] for v in values]) for key in values[0]}

    def test_all_good_and_manual_subset_parity(self):
        for failed in [(), (2,), (2, 5)]:
            selected = self.selection(failed)
            source = evaluate_retained_sources(selected, self.outcomes)
            target = self.target(selected)
            direct = conditional_linear_bands(source, target, 6-len(failed), dispersion_failure=.015)
            actual = conditional_selected_bands(selected, source, target)
            for key in direct:
                np.testing.assert_array_equal(actual[key], direct[key])
            for i, original in enumerate(selected['retained_indices']):
                manual = evaluate_linear_balance(selected['designs'][i][0]['design'], self.outcomes[original][::2])
                self.assertEqual(source['estimate'][i, 0], manual['estimate'])

    def test_selection_is_outcome_independent_and_counts_are_worst_case(self):
        selected = self.selection((1, 4), valid=[5, 6])
        first = evaluate_retained_sources(selected, self.outcomes)
        second = evaluate_retained_sources(selected, [1-y for y in self.outcomes])
        np.testing.assert_allclose(first['estimate']+second['estimate'], 1.)
        np.testing.assert_array_equal(selected['retained_indices'], [0, 2, 3, 5, 6, 7])
        for valid_size, adjusted in zip([5, 6], selected['adjusted_valid']):
            for valid_set in itertools.combinations(range(8), valid_size):
                actual_valid = len(set(valid_set).intersection(selected['retained_indices']))
                self.assertGreaterEqual(actual_valid, adjusted)
        self.assertEqual(selected['adjusted_valid'].tolist(), [3, 4])

    def test_no_valid_guarantee_and_all_required_fall_back(self):
        target = self.target(self.selection())
        for selected in [self.selection(range(8)), self.selection((1, 2), valid=2),
                         self.selection((1,), policy='all_required')]:
            self.assertFalse(selected['eligible'])
            summary = evaluate_retained_sources(selected, self.outcomes)
            result = conditional_selected_bands(selected, summary, target)
            self.assertTrue(result['design_fallback'])
            np.testing.assert_array_equal(result['conditional_interval'], result['target_band'])
        with self.assertRaises(ValueError):
            self.selection(valid=6.5)
        with self.assertRaises(ValueError):
            self.selection(policy='outcome_screening')

    def test_projection_budget_recomputed_on_retained_set(self):
        selected = self.selection((1, 6))
        source = evaluate_retained_sources(selected, self.outcomes)
        budgets = retained_projection_budgets(selected, source, 'bounded_invalid')
        # K'=6, g'=4; the maximum invalid count is still 2, not floor(K'/4).
        for arm, name in enumerate(['mu1', 'mu0']):
            summaries = {key: value[:, arm] for key, value in source.items()}
            geometry = [projection_geometry(arms[arm]['design']) for arms in selected['designs']]
            for key in ['sample_size', 'projection_commutator_norm']:
                summaries[key] = np.array([g[key] for g in geometry])
            direct = projection_error_budget(summaries, .5*np.sqrt(summaries['sample_size']), 2)
            self.assertEqual(budgets[name], direct)
        declared = retained_projection_budgets(selected, source, 'declared', [.5, .5])
        self.assertEqual(declared, budgets)
        with self.assertRaises(ValueError):
            retained_projection_budgets(selected, source, 'declared')
        with self.assertRaises(ValueError):
            conditional_selected_bands(selected, None, self.target(selected))


if __name__ == '__main__':
    unittest.main()
