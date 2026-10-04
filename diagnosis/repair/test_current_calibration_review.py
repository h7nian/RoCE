#!/usr/bin/env python3
"""Check pairing guards and independent numerical references for the MC200 review."""
import copy
import csv
import math
from pathlib import Path
import statistics
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
import review_current_calibration as review


class CalibrationReviewTests(unittest.TestCase):
    def configuration(self):
        return dict(crossfit_layers=2, source_validation_method="outer_fit",
                    n_per_site=1000, aggregation_lambda=.5, dgp_type="bounded",
                    dgp_control=dict(covariate_shift=1, source_treatment_scale=1))

    def test_configuration_guards_keep_the_other_component_fixed(self):
        full = self.configuration()
        source = dict(full, source_nuisance_method="standard")
        target = dict(full, target_nuisance_method="lasso")
        review.validate_configurations(full, source, "source")
        review.validate_configurations(full, target, "target")
        for changed in (dict(source, crossfit_layers=3), dict(source, aggregation_lambda=1),
                        dict(source, target_nuisance_method="lasso")):
            with self.assertRaises(ValueError):
                review.validate_configurations(full, changed, "source")

    def test_prescribed_seeds_and_all_scenarios_are_required(self):
        tasks = [dict(config=config, K="2", p="100", rho="0", sim_id=str(seed),
                      protocol="one_round", deviation_mechanism="treated_arm")
                 for config in ("C1", "C2", "C3") for seed in range(1, 201)]
        self.assertEqual(len(review.index_tasks(tasks, self.configuration(), 100)), 600)
        for changed in (tasks[:-1], tasks + [tasks[0]],
                        tasks[:-1] + [dict(tasks[-1], sim_id="401")]):
            with self.assertRaises(ValueError):
                review.index_tasks(changed, self.configuration(), 100)

    def result_rows(self):
        value = dict(estimate=.2, se=.03, truth=.24)
        for arm in ("mu1", "mu0"):
            for source in ("s1", "s2"):
                value[arm + "_source_" + source + "_source_estimate"] = .4
        return {method: dict(value) for method in
                list(review.MODES.values()) + ["target_only_ate", "target_anchor_ate"]}

    def test_reference_guards_catch_cross_dataset_and_cross_component_changes(self):
        full, ablation = self.result_rows(), self.result_rows()
        review.validate_references(full, ablation, "abc", "abc", "source")
        with self.assertRaisesRegex(ValueError, "hashes"):
            review.validate_references(full, ablation, "abc", "def", "source")
        ablation["target_anchor_ate"]["se"] = .04
        with self.assertRaisesRegex(ValueError, "calibrated-target"):
            review.validate_references(full, ablation, "abc", "abc", "source")
        ablation = self.result_rows()
        del ablation["target_anchor_ate"]
        review.validate_references(full, ablation, "abc", "abc", "target")
        ablation[review.MODES["joint_tate"]]["mu1_source_s1_source_estimate"] = .5
        with self.assertRaisesRegex(ValueError, "source_estimate"):
            review.validate_references(full, ablation, "abc", "abc", "target")

    def test_arm_covariance_must_reconstruct_saved_tate_variance(self):
        results = {method: dict(estimate=2, truth=1, se=2) for method in review.MODES.values()}
        rows = [dict(recipe="score_derivative", aggregation_mode=mode, crossfit_layers=2,
                     quantity=quantity, estimate=estimate, truth=truth, variance=variance,
                     variance_fixed_weights=variance, covariance_mu1_mu0=.5,
                     covariance_mu1_mu0_fixed_weights=.5)
                for mode in review.MODES for quantity, estimate, truth, variance in
                (("mu1", 2.5, 1.5, 1), ("mu0", .5, .5, 4), ("tate", 2, 1, 4))]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            def write(values):
                with (path / "layer_estimates.csv").open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(values[0]))
                    writer.writeheader()
                    writer.writerows(values)
            write(rows)
            self.assertEqual(len(review.read_layers(path, results)), 9)
            changed = copy.deepcopy(rows)
            changed[0]["covariance_mu1_mu0"] = .25
            write(changed)
            with self.assertRaisesRegex(ValueError, "covariance"):
                review.read_layers(path, results)
            write(rows[:-1])
            with self.assertRaisesRegex(ValueError, "Missing"):
                review.read_layers(path, results)

    def test_paired_statistics_against_direct_delete_one_calculation(self):
        first, second = [1, -1, 2, -2], [.5, -.5, 3, -3]
        first_variance, second_variance = [2, 3, 4, 5], [1, 2, 6, 7]
        full = [dict(estimate=x, truth=0, variance=v) for x, v in zip(first, first_variance)]
        other = [dict(estimate=x, truth=0, variance=v) for x, v in zip(second, second_variance)]
        value = review.paired_summary(full, other)
        self.assertAlmostEqual(value["mse_difference"], -2.125)
        self.assertAlmostEqual(value["mse_difference_mcse"], statistics.stdev([.75, .75, -5, -5])/2)
        self.assertAlmostEqual(value["relative_rmse_change"], math.sqrt(2.5/4.625)-1)
        deleted = []
        for index in range(4):
            keep = [j for j in range(4) if j != index]
            deleted.append(statistics.mean(first_variance[j]-second_variance[j] for j in keep) -
                           (statistics.variance(first[j] for j in keep) -
                            statistics.variance(second[j] for j in keep)))
        expected = math.sqrt(.75*sum((x-statistics.mean(deleted))**2 for x in deleted))
        self.assertAlmostEqual(value["variance_discrepancy_jackknife_se"], expected)
        self.assertEqual(review.paired_summary(full, full)["mse_difference"], 0)
        other[0]["truth"] = .01
        with self.assertRaisesRegex(ValueError, "truth"):
            review.paired_summary(full, other)


if __name__ == "__main__":
    unittest.main()
