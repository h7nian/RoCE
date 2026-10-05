import unittest
import numpy as np

from design_selection import restrict_weight_inflation, evaluate_retained_sources, conditional_selected_bands
from design_budget import inference_ledger, apply_design_precision_gate, population_selected_interval
from population_protection import target_coefficient_certificate, composition_radius, expand_conditional_interval
import test_design_selection


class DesignBudgetTests(unittest.TestCase):
    def setUp(self):
        self.fixture = test_design_selection.DesignSelectionTests()
        self.fixture.setUp()

    def test_common_event_ledger_and_exact_target_fallback(self):
        for policy in ['legacy', 'conditional_reallocation', 'shared_composition']:
            borrow = inference_ledger(True, policy)
            fallback = inference_ledger(False, policy)
            self.assertEqual(borrow['composition'], fallback['composition'])
            for row in [borrow, fallback]:
                self.assertLessEqual(sum(row[k] for k in ['source','target','dispersion','coefficient','composition']), .05+1e-14)
        selected = self.fixture.selection((1, 2), valid=2)
        target = self.fixture.target(self.fixture.selection())
        matrix, treatment = self.fixture.designs[0]
        certificate = target_coefficient_certificate(matrix[:, 1:], treatment, self.fixture.outcomes[0])
        result = population_selected_interval(selected, None, target, certificate, len(matrix))
        point = target['estimate'][0]-target['estimate'][1]
        noise = np.sqrt(target['weight_square_sum'].sum()/2*np.log(2/.0225))
        expected = expand_conditional_interval(point+noise*np.array([-1.,1.]),
                                               composition_radius(certificate, len(matrix), .0225))
        np.testing.assert_array_equal(result['interval'], expected)

    def test_legacy_parity_and_weight_guard_counts(self):
        selected = self.fixture.selection((1,))
        target = self.fixture.target(selected)
        source = evaluate_retained_sources(selected, self.fixture.outcomes)
        ledger = inference_ledger(True, 'legacy')
        old = conditional_selected_bands(selected, source, target)
        controlled = conditional_selected_bands(selected, source, target,
            source_failure=ledger['source'], target_failure=ledger['target'], dispersion_failure=ledger['dispersion'])
        for name in old:
            np.testing.assert_array_equal(old[name], controlled[name])
        kept = restrict_weight_inflation(selected, 100.)
        np.testing.assert_array_equal(kept['retained_indices'], selected['retained_indices'])
        dropped = restrict_weight_inflation(selected, 1.)
        self.assertFalse(dropped['eligible'])
        np.testing.assert_array_equal(dropped['adjusted_valid'], np.maximum(0, selected['original_valid']-len(dropped['rejected'])))
        self.assertEqual(len(selected['rejected']), 1)  # No mutation of the original policy.

    def test_precision_gate_uses_only_design_and_can_abstain(self):
        selected = self.fixture.selection()
        target = self.fixture.target(selected)
        gate = apply_design_precision_gate(selected, target['weight_square_sum'])
        np.testing.assert_array_equal(gate['retained_indices'], selected['retained_indices'])
        self.assertEqual(gate['precision_gate_rejected'],
                         gate['design_source_radius'] >= gate['design_target_radius'])
        very_precise_target = apply_design_precision_gate(selected, np.array([1e-12,1e-12]))
        self.assertFalse(very_precise_target['eligible'])


if __name__ == '__main__':
    unittest.main()
