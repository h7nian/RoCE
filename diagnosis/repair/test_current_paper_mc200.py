#!/usr/bin/env python3
"""Check the expansion's scientific identities, missing seeds and boundaries."""
from collections import Counter
import sys
import unittest

sys.dont_write_bytecode = True
import prepare_current_paper_mc200 as plan
import submit_repeat_pilot as pilot
from prepare_highdim_validation import experiment_identity


class CurrentPaperPlanTests(unittest.TestCase):
    def test_full_grid_has_every_prescribed_pair_and_separate_boundaries(self):
        requests, reused, missing = plan.resolve_requests({})
        self.assertFalse(reused)
        self.assertEqual(Counter(row["role"] for row in requests),
                         dict(method=46400, baseline=46400, standard_source=1200, ordinary_target=1200))
        counts = Counter((row["role"], row["task"]["panel"]) for row in requests)
        self.assertEqual(counts[("method", "main")], 24000)
        self.assertEqual(counts[("method", "c4")], 8000)
        self.assertEqual(counts[("method", "half_invalid")], 14400)
        identities = set()
        for request in requests:
            task = request["task"]
            configuration = plan.scientific_configuration(request["role"], task["n_deviated_sites"])
            key = experiment_identity(task, configuration)
            self.assertNotIn(key, identities)
            identities.add(key)
            if task["panel"] == "half_invalid":
                self.assertEqual(task["n_deviated_sites"], task["K"]//2)
                self.assertNotEqual(task["rho"], 0)
        main = [row["task"] for row in requests if row["role"] == "method" and row["task"]["panel"] == "main"]
        cells = Counter((row["p"], row["K"], row["config"], row["rho"]) for row in main)
        self.assertEqual(len(cells), 120)
        self.assertEqual(set(cells.values()), {200})
        self.assertEqual({key[3] for key in cells}, {0, .5, 1, 1.5, 2})
        self.assertEqual(sum(len(group["configs"])*len(group["rhos"])*len(group["seeds"])
            for group in plan.rectangular_groups(missing["method"])), 46400)

    def test_existing_baselines_remove_only_exact_data_and_seeds(self):
        rows = pilot.task_rows(list(plan.SCENARIOS), [4], [0, .5, 1, 2], range(1, 51),
                               "one_round", "treated_arm", [100])
        existing = {experiment_identity(row, plan.scientific_configuration("baseline")):
                    dict(root="frozen", task_id=row["task_id"], complete=True) for row in rows}
        _, reused, missing = plan.resolve_requests(existing)
        self.assertEqual(len(reused), 600)
        self.assertEqual({row["role"] for row in reused}, {"baseline"})
        self.assertEqual(len(missing["method"]), 46400)
        groups = [group for group in plan.rectangular_groups(missing["baseline"])
                  if group["panel"] == "main" and group["dimension"] == 100 and group["sources"] == 4]
        self.assertEqual(len(groups), 2)
        by_rhos = {tuple(group["rhos"]): group["seeds"] for group in groups}
        self.assertEqual(by_rhos[(1.5,)], list(range(1, 201)))
        self.assertEqual(by_rhos[(0, .5, 1, 2)], list(range(51, 201)))

    def test_three_layer_and_half_invalid_results_cannot_replace_main(self):
        row = pilot.task_rows(["C1"], [4], [1], [1], "one_round", "treated_arm", [100])[0]
        current = plan.scientific_configuration()
        three = dict(current, crossfit_layers=3, source_validation_method="calibrated")
        half = dict(current, n_deviated_sites=2)
        self.assertNotEqual(experiment_identity(row, current), experiment_identity(row, three))
        self.assertNotEqual(experiment_identity(row, current), experiment_identity(row, half))
        key = experiment_identity(row, three)
        _, reused, _ = plan.resolve_requests({key: dict(root="old_three", task_id=1, complete=True)})
        self.assertFalse(reused)

    def test_unfinished_existing_work_is_not_duplicated(self):
        row = pilot.task_rows(["C1"], [2], [0], [1], "one_round", "treated_arm", [100])[0]
        key = experiment_identity(row, plan.scientific_configuration("baseline"))
        with self.assertRaisesRegex(ValueError, "needs recovery"):
            plan.resolve_requests({key: dict(root="pending", task_id=1, complete=False)})

    def test_grouping_rejects_missing_cartesian_cells_and_duplicate_seeds(self):
        rows = pilot.task_rows(["C1", "C2"], [2], [0, 1], [1], "one_round", "treated_arm", [100])
        for row in rows:
            row.update(panel="main", n_deviated_sites=1)
        with self.assertRaisesRegex(ValueError, "rectangular"):
            plan.rectangular_groups(rows[:-1])
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            plan.rectangular_groups(rows + [rows[0]])


if __name__ == "__main__":
    unittest.main()
