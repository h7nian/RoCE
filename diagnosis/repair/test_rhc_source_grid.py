"""Guard paired source-rule experiments and the unchanged default manifest order."""
import unittest

from prepare_rhc_current import prepare, source_rule_grid


class SourceRuleGridTests(unittest.TestCase):
    def test_original_design_retains_profile_order(self):
        grid = source_rule_grid()
        self.assertEqual(len(grid), 24)
        self.assertEqual([(seed, name) for seed, name, _ in grid[:6]], [
            (0, 'min'), (0, 'final_weight_1se'), (0, 'final_outcome_1se'),
            (0, 'initial_weight_1se'), (0, 'source_all_1se'), (0, 'global_1se')])
        self.assertEqual([grid[index][0] for index in (0, 6, 12, 18)], [0, 101, 202, 303])

    def test_confirmation_pairs_every_rule_with_every_seed(self):
        variants = ['min', 'source_all_1se', 'global_1se']
        grid = source_rule_grid(range(401, 421), variants)
        self.assertEqual(len(grid), 60)
        self.assertEqual({(seed, name) for seed, name, _ in grid},
                         {(seed, name) for seed in range(401, 421) for name in variants})
        for _, name, rules in grid:
            if name == 'source_all_1se':
                self.assertEqual(rules, dict(initial_weight='1se', weight='1se', outcome='1se'))
            else:
                self.assertIsNone(rules)

    def test_invalid_design_fails_before_preparation(self):
        for seeds in ([], [1, 1], [-1], [True], [1.5], [2147483648]):
            with self.subTest(seeds=seeds), self.assertRaises(ValueError):
                source_rule_grid(seeds)
        for variants in ([], ['min', 'min'], ['unknown']):
            with self.subTest(variants=variants), self.assertRaises(ValueError):
                source_rule_grid(variants=variants)

    def test_unknown_target_is_rejected_before_creating_files(self):
        with self.assertRaisesRegex(ValueError, 'four retained RHC insurance sites'):
            prepare('/unused', '/unused', '/unused', target_site='No insurance')

    def test_truncation_settings_cannot_silently_change_other_profiles(self):
        with self.assertRaisesRegex(ValueError, 'Source radii require'):
            prepare('/unused', '/unused', '/unused', source_radii=[3])
        with self.assertRaisesRegex(ValueError, 'keeps the min rule'):
            prepare('/unused', '/unused', '/unused', study_profile='source_truncation', source_variants=['global_1se'])
        for radii in ([], [3, 3], [0], [12]):
            with self.subTest(radii=radii), self.assertRaisesRegex(ValueError, 'prespecified source radii'):
                prepare('/unused', '/unused', '/unused', study_profile='source_truncation', source_radii=radii)


if __name__ == '__main__':
    unittest.main()
