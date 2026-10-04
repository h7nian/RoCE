#!/usr/bin/env python3
"""Exercise duplicate and provenance guards without contacting Slurm."""
import contextlib
import io
import json
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
import supervise_repeat_pilot as supervisor
import summarize_repeat_pilot as review
import advance_validation_plan as expansion
import prepare_highdim_validation as highdim
import prepare_layer_validation as layers
import review_layer_comparison as layer_review


class RepeatPilotTests(unittest.TestCase):
    def test_certificate_control_requires_checked_gate_and_enters_provenance(self):
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            scratch = Path(directory)
            check = scratch / "check"
            (check / "Rlib/RoCE").mkdir(parents=True)
            (check / "Rlib/RoCE/DESCRIPTION").write_text("test package")
            (check / "TESTS_PASSED").touch()
            (check / "source_manifest.json").write_text("[]")
            args = SimpleNamespace(output=str(scratch / "method"), check_root=str(check),
                configs=["C2"], sources=[2], rhos=[0], repeats=[1], protocol="one_round",
                deviation_mechanism="treated_arm", recipes=None, cpus=2, memory="8G",
                time_limit="01:00:00", partitions="preempt", account="hou00123",
                shared_cache=False, nuisance_cv_certificate=True)
            with mock.patch.object(pilot, "SCRATCH", scratch):
                with self.assertRaisesRegex(ValueError, "full CV certificate checks"):
                    pilot.prepare(args)
                self.assertFalse((scratch / "method").exists())
                gate = check / "CV_CERTIFICATE_CHECKS_PASSED"
                gate.write_text("checked")
                pilot.prepare(args)
                configuration = pilot.verify(scratch / "method")
                self.assertIs(configuration["nuisance_cv_certificate"], True)
                provenance = json.loads((scratch / "method/provenance.json").read_text())
                self.assertIn(str(gate), provenance)
                args.nuisance_cv_certificate = "TRUE"
                with self.assertRaisesRegex(ValueError, "must be boolean"):
                    pilot.prepare(args)
                args.nuisance_cv_certificate = True
                args.output = str(scratch / "baseline")
                args.baseline_methods = list(pilot.BASELINE_METHODS)
                pilot.prepare(args)
                baseline = pilot.verify(scratch / "baseline")
                self.assertTrue(baseline["nuisance_cv_certificate"])
                self.assertEqual(baseline["worker_script"], "run_baseline_repeat.R")

    def test_variance_component_uncertainty_matches_brute_force_deletions(self):
        first, second, reported = [1,3,-1,4], [2,-2,0,5], [.1,.3,.2,.5]
        def discrepancy(x,y,values):
            import statistics
            mx,my = statistics.mean(x),statistics.mean(y)
            empirical = sum((a-mx)*(b-my) for a,b in zip(x,y))/(len(x)-1)
            return statistics.mean(values)-empirical
        result = layer_review.covariance_accuracy(first,second,reported)
        deleted = [discrepancy(first[:i]+first[i+1:],second[:i]+second[i+1:],reported[:i]+reported[i+1:]) for i in range(4)]
        center = sum(deleted)/4
        self.assertAlmostEqual(result["reported_minus_empirical"],discrepancy(first,second,reported))
        self.assertAlmostEqual(result["discrepancy_jackknife_se"],(.75*sum((x-center)**2 for x in deleted))**.5)
        diagonal = layer_review.covariance_accuracy(first,first,reported)
        self.assertAlmostEqual(diagonal["empirical_component"],layer_review.statistics.variance(first))
        self.assertIsNone(layer_review.covariance_accuracy(first[:2],second[:2],reported[:2])["discrepancy_jackknife_se"])
        with self.assertRaises(ValueError):
            layer_review.covariance_accuracy([1],[1],[1])

    def test_paired_variance_discrepancy_has_independent_jackknife_check(self):
        result = layer_review.paired_variance_comparison([-1,0,1],[-.5,0,.5],[1,1,1],[.25,.25,.25])
        self.assertEqual(result["empirical_variance_difference_two_minus_three"],.75)
        self.assertEqual(result["variance_difference_discrepancy"],0)
        self.assertAlmostEqual(result["variance_discrepancy_jackknife_se"],.75)
        self.assertIsNone(layer_review.paired_variance_comparison([0],[0],[1],[1])["empirical_variance_difference_two_minus_three"])

    def test_completed_local_import_keeps_its_actual_execution_record(self):
        with tempfile.TemporaryDirectory() as directory:
            scratch = Path(directory)
            old, new = scratch / "local", scratch / "campaign"
            for root in (old, new):
                (root / "tasks").mkdir(parents=True)
                (root / "submissions").mkdir()
                (root / "configuration.json").write_text("{}")
                (root / "manifest.csv").write_text("task_id,config\n1,C1\n")
            task = old / "tasks/1"
            task.mkdir()
            for name in ("COMPLETE", "results.csv", "data_sha256.txt"):
                (task / name).write_text("committed fixture")
            (task / "status.txt").write_text("COMPLETE")
            with mock.patch.object(pilot, "SCRATCH", scratch), mock.patch.object(pilot, "verify", return_value={}):
                with self.assertRaisesRegex(ValueError, "local validation record"):
                    pilot.import_completed(new, old, [dict(task_id=1)])
                (old / "local_validation.json").write_text(json.dumps(dict(backend="local", exit_code=0, completed_tasks=[1])))
                pilot.import_completed(new, old, [dict(task_id=1)])
                self.assertTrue((new / "tasks/1").is_symlink())
                self.assertFalse((new / "submissions/1.json").exists())
                receipt = json.loads((new / "imported_results.json").read_text())["tasks"][0]
                self.assertEqual(receipt["backend"], "local")
                self.assertEqual(receipt["receipt"], str(old / "local_validation.json"))

    def test_controller_finishes_when_all_results_were_imported_from_local_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "tasks/1").mkdir(parents=True)
            (root / "submissions").mkdir()
            (root / "tasks/1/COMPLETE").touch()
            (root / "configuration.json").write_text("{}")
            (root / "manifest.csv").write_text("task_id,sim_id\n1,401\n")
            rows = [dict(task_id="1", sim_id="401")]
            self.assertEqual(supervisor.submission_capacity(rows, {}, {}, 32, 401, {"1"}), 0)
            with mock.patch.object(pilot, "verify", return_value=dict(first_seed=401)), \
                 mock.patch.object(supervisor, "scheduler_states", return_value={}), \
                 mock.patch.object(review, "summarize", return_value=True), \
                 mock.patch.object(pilot, "submit") as submit, contextlib.redirect_stdout(io.StringIO()):
                self.assertTrue(supervisor.supervise(root, 32, .1, 30))
                submit.assert_not_called()

    def test_layer_metrics_use_empirical_variance_and_keep_unavailable_repeats(self):
        rows = [dict(estimate=value, truth=0, variance=.01, variance_fixed_weights=.0025)
                for value in (-.1, 0, .1)]
        result = layer_review.metrics(rows, 4)
        self.assertEqual(result["n_unavailable"], 1)
        self.assertAlmostEqual(result["empirical_variance"], .01)
        self.assertAlmostEqual(result["variance_ratio"], 1)
        self.assertEqual(result["coverage"], 1)
        self.assertAlmostEqual(result["coverage_fixed_weights"], 1/3)
        self.assertIsNone(layer_review.metrics(rows[:1], 4)["empirical_variance"])
        with self.assertRaisesRegex(ValueError, "Population truth"):
            layer_review.metrics([dict(rows[0], truth=1), rows[1]], 2)

    def test_layer_report_checks_pairing_and_reconstructs_arm_covariance(self):
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            scratch = Path(directory)
            plan = scratch / "plan"
            plan.mkdir()
            campaigns, configurations = [], {}
            for count in (2, 3):
                root = plan / ("layers" + str(count))
                task = root / "tasks/1"
                task.mkdir(parents=True)
                record = dict(task_id="1", p="10", K="1", config="C1", rho="0",
                              protocol="one_round", deviation_mechanism="treated_arm", sim_id="401")
                review.write_csv(root / "manifest.csv", [record, dict(record,task_id="2",sim_id="402")])
                (task / "COMPLETE").touch()
                (task / "status.txt").write_text("COMPLETE")
                (task / "data_sha256.txt").write_text("a"*64)
                (task / "score_derivative_outer_fit_sha256.txt").write_text("b"*64)
                rows = []
                for mode in layer_review.MODES:
                    for quantity, estimate, variance, fixed in (("mu1", .5, .04, .03),
                                                               ("mu0", .2, .01, .02), ("tate", .3, .04, .04)):
                        rows.append(dict(recipe="score_derivative", aggregation_mode=mode,
                            quantity=quantity, crossfit_layers=count, estimate=estimate, truth=estimate,
                            variance=variance, variance_fixed_weights=fixed,
                            covariance_mu1_mu0=.005, covariance_mu1_mu0_fixed_weights=.005))
                review.write_csv(task / "layer_estimates.csv", rows)
                review.write_csv(task / "layer_weights.csv", [dict(recipe="score_derivative",
                    aggregation_mode="common_tate", coordinate="s1", weight=.5)])
                review.write_csv(task / "results.csv", [dict(simulation_elapsed_seconds=10)])
                if count == 2:
                    shutil.copytree(str(task),str(root / "tasks/2"))
                configurations[root] = dict(crossfit_layers=count, recipes=["score_derivative"])
                campaigns.append(dict(root=str(root)))
            (plan / "expansion_plan.json").write_text(json.dumps(dict(campaigns=campaigns)))
            with mock.patch.object(pilot, "SCRATCH", scratch), \
                 mock.patch.object(pilot, "verify", side_effect=lambda root: configurations[root]):
                layer_review.review(plan, scratch / "report")
                covariance = layer_review.read_csv(scratch / "report/arm_covariance.csv")
                self.assertEqual(len(covariance), 6)
                self.assertTrue(all(float(row["mean_covariance"]) == .005 for row in covariance))
                pairs = layer_review.read_csv(scratch / "report/paired_differences.csv")
                self.assertEqual(len(pairs), 9)
                self.assertTrue(all(row["n_pairs"] == "1" for row in pairs))
                paired_metrics = layer_review.read_csv(scratch / "report/paired_metrics.csv")
                self.assertTrue(all(row["n_success"] == "1" for row in paired_metrics))
                self.assertTrue(all(row["n_complete_in_layer"] == "2" for row in paired_metrics if row["crossfit_layers"] == "2"))
                paired_covariance = layer_review.read_csv(scratch / "report/paired_arm_covariance.csv")
                self.assertEqual(len(paired_covariance),6)
                self.assertTrue(all(row["n_success"] == "1" for row in paired_covariance))
                (plan / "layers3/tasks/1/score_derivative_outer_fit_sha256.txt").write_text("c"*64)
                with self.assertRaisesRegex(ValueError, "outer fits differ"):
                    layer_review.review(plan, scratch / "mismatched_report")

    def test_layer_profiles_separate_site_count_and_covariate_dimension(self):
        profiles = layers.profiles()
        self.assertEqual(len(profiles), 6)
        self.assertEqual(sum(profile["max_in_flight"] for profile in profiles), 256)
        for dimension in (10, 100, 200):
            pair = [profile for profile in profiles if profile["dimension"] == dimension]
            self.assertEqual([profile["crossfit_layers"] for profile in pair], [2, 3])
            for profile in pair:
                self.assertEqual(profile["sources"], [2, 4, 6])
                self.assertEqual(profile["configs"], ["C1", "C2", "C3"])
                self.assertEqual(profile["seed_start"], 401)
                self.assertEqual(profile["n_repeats"], 50)
        with self.assertRaises(ValueError):
            layers.profiles(n_repeats=0)

    def test_layer_pair_gate_rejects_different_outer_fits(self):
        import csv
        with tempfile.TemporaryDirectory() as directory:
            scratch = Path(directory)
            roots = [scratch / "two", scratch / "three"]
            configs = {}
            for root, count in zip(roots, (2, 3)):
                task = root / "tasks/1"
                task.mkdir(parents=True)
                (root / "manifest.csv").write_text("task_id,protocol\n1,one_round\n")
                (task / "COMPLETE").touch()
                (task / "status.txt").write_text("COMPLETE\n")
                (task / "data_sha256.txt").write_text("a" * 64)
                (task / "score_derivative_outer_fit_sha256.txt").write_text("b" * 64)
                rows = [dict(estimand_scope="tate", pilot_recipe="score_derivative", method=method,
                             estimate=.9, truth=.1, se=.01) for method in
                        ("one_round_crossfit_ate", "one_round_crossfit_ate_separate_arms",
                         "one_round_crossfit_ate_joint_tate", "target_only_ate")]
                with (task / "results.csv").open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                    writer.writeheader()
                    writer.writerows(rows)
                configs[root] = dict(recipes=["score_derivative"], crossfit_layers=count)
            gate = dict(root=str(roots[0]), paired_root=str(roots[1]), tasks=[1])
            with mock.patch.object(pilot, "SCRATCH", scratch), \
                 mock.patch.object(pilot, "verify", side_effect=lambda root: configs[root]):
                self.assertTrue(expansion.check_gate(gate)[0])
                (roots[1] / "tasks/1/score_derivative_outer_fit_sha256.txt").write_text("c" * 64)
                with self.assertRaisesRegex(ValueError, "outer nuisance fits differ"):
                    expansion.check_gate(gate)
                (roots[1] / "tasks/1/COMPLETE").unlink()
                self.assertFalse(expansion.check_gate(gate)[0])

    def test_layer_cli_and_scientific_identity(self):
        command = ["pilot", "prepare", "unused", "--check-root", "unused", "--crossfit-layers", "2"]
        with mock.patch.object(sys, "argv", command), mock.patch.object(pilot, "prepare") as prepare:
            pilot.main()
            self.assertEqual(prepare.call_args[0][0].crossfit_layers, 2)
        row = pilot.task_rows(["C1"], [2], [0], [401], "one_round", "treated_arm", [10])[0]
        old = highdim.experiment_identity(row, {})
        self.assertEqual(old, highdim.experiment_identity(row, dict(crossfit_layers=3)))
        two = highdim.experiment_identity(row, dict(crossfit_layers=2))
        self.assertNotEqual(two, old)
        self.assertEqual(two, highdim.experiment_identity(row, dict(source_validation_method="outer_fit")))
        self.assertEqual(highdim.data_identity(row, dict(crossfit_layers=2)),
                         highdim.data_identity(row, dict(crossfit_layers=3)))
        with self.assertRaisesRegex(ValueError, "requires crossfit_layers=2"):
            highdim.experiment_identity(row, dict(crossfit_layers=3, source_validation_method="outer_fit"))

    def test_layer_preparation_requires_checked_capability_and_records_configuration(self):
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            scratch = Path(directory)
            check = scratch / "check"
            (check / "Rlib/RoCE").mkdir(parents=True)
            (check / "Rlib/RoCE/DESCRIPTION").write_text("test package")
            (check / "TESTS_PASSED").touch()
            (check / "source_manifest.json").write_text("[]")
            args = SimpleNamespace(output=str(scratch / "two"), check_root=str(check),
                configs=["C1"], sources=[2], rhos=[0], repeats=[401], protocol="one_round",
                deviation_mechanism="treated_arm", recipes=None, cpus=2, memory="4G",
                time_limit="00:30:00", partitions="preempt", account="hou00123",
                shared_cache=False, crossfit_layers=2, source_validation_method="calibrated")
            with mock.patch.object(pilot, "SCRATCH", scratch):
                with self.assertRaisesRegex(ValueError, "two/three-layer regression"):
                    pilot.prepare(args)
                (check / "CROSSFIT_LAYER_CHECKS_PASSED").write_text("passed")
                args.source_validation_method = "initial"
                with self.assertRaisesRegex(ValueError, "legacy initial"):
                    pilot.prepare(args)
                args.source_validation_method = "calibrated"
                pilot.prepare(args)
                root = scratch / "two"
                configuration = pilot.verify(root)
                self.assertEqual(configuration["crossfit_layers"], 2)
                self.assertEqual(configuration["source_validation_method"], "outer_fit")
                self.assertTrue(configuration["layer_diagnostics"])
                self.assertTrue((root / "workflow/layer_comparison_summary.R").is_file())
                self.assertEqual(configuration["n_folds"], 10)
                self.assertEqual(configuration["nlambda"], 100)
                args.baseline_methods = list(pilot.BASELINE_METHODS)
                with self.assertRaisesRegex(ValueError, "baseline worker"):
                    pilot.prepare(args)

    def test_summary_keeps_layer_and_validation_schemes_separate(self):
        import csv
        row = pilot.task_rows(["C1"], [2], [0], [401], "one_round", "treated_arm", [10])[0]
        manifest = [dict(row, task_id=str(index), crossfit_layers=str(layers), source_validation_method=method)
                    for index, (layers, method) in enumerate(
                        ((2, "outer_fit"), (3, "calibrated"), (3, "initial")), start=1)]
        estimates = [dict(task, state="COMPLETE", recipe="score_derivative", method="target_only_ate",
                          estimate=.2, truth=.2, se=.1) for task in manifest]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            review.summarize_metrics(manifest, ["score_derivative"], estimates, root, ["target_only_ate"])
            with (root / "metrics.csv").open() as stream:
                metrics = list(csv.DictReader(stream))
            self.assertEqual(len(metrics), 3)
            self.assertTrue(all(row["n_planned"] == "1" and row["n_success"] == "1" for row in metrics))

    def test_calibration_profiles_fill_only_missing_factorial_cells(self):
        profiles = highdim.profiles("calibration")
        self.assertEqual(len(profiles), 6)
        self.assertEqual(sum(len(highdim.rows_for_profile(profile)) for profile in profiles), 606)
        self.assertEqual(sum(profile["max_in_flight"] for profile in profiles), 134)
        for profile in profiles:
            self.assertEqual(profile["source_nuisance_method"], "standard")
            self.assertEqual(profile["source_validation_method"], "complete")
            self.assertEqual(profile["rhos"], [0])
            self.assertEqual(profile["baseline_max_in_flight"], 0)
            if profile["dimensions"] != [10]:
                self.assertIn(profile["parent"], {entry["name"] for entry in profiles})
        self.assertEqual({profile["target_nuisance_method"] for profile in profiles},
                         {"hou_calibrated", "lasso"})

    def test_source_program_cli_and_duplicate_identity(self):
        command = ["pilot", "prepare", "unused", "--check-root", "unused",
                   "--source-nuisance-method", "standard", "--source-validation-method", "complete"]
        with mock.patch.object(sys, "argv", command), mock.patch.object(pilot, "prepare") as prepare:
            pilot.main()
            self.assertEqual(prepare.call_args[0][0].source_nuisance_method, "standard")
            self.assertEqual(prepare.call_args[0][0].source_validation_method, "complete")
        row = pilot.task_rows(["C1"], [2], [0], [1], "one_round", "treated_arm", [10])[0]
        self.assertEqual(highdim.experiment_identity(row, {}),
                         highdim.experiment_identity(row, dict(source_validation_method="complete")))
        self.assertNotEqual(highdim.experiment_identity(row, {}),
                            highdim.experiment_identity(row, dict(source_nuisance_method="standard")))
        self.assertEqual(highdim.data_identity(row, {}),
                         highdim.data_identity(row, dict(source_nuisance_method="standard")))

    def test_standard_source_preparation_requires_checked_capability(self):
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            scratch = Path(directory)
            check = scratch / "check"
            (check / "Rlib/RoCE").mkdir(parents=True)
            (check / "Rlib/RoCE/DESCRIPTION").write_text("test package")
            (check / "TESTS_PASSED").touch()
            (check / "source_manifest.json").write_text("[]")
            args = SimpleNamespace(output=str(scratch / "standard"), check_root=str(check),
                configs=["C1"], sources=[2], rhos=[0], repeats=[1], protocol="one_round",
                deviation_mechanism="treated_arm", recipes=None, cpus=2, memory="4G",
                time_limit="00:30:00", partitions="preempt", account="hou00123",
                shared_cache=False, source_nuisance_method="standard")
            with mock.patch.object(pilot, "SCRATCH", scratch):
                with self.assertRaisesRegex(ValueError, "standard-source regression"):
                    pilot.prepare(args)
                (check / "SOURCE_STANDARD_CHECKS_PASSED").write_text("passed")
                args.protocol = "two_round"
                with self.assertRaisesRegex(ValueError, "one_round"):
                    pilot.prepare(args)
                args.protocol = "one_round"
                args.source_validation_method = "initial"
                with self.assertRaisesRegex(ValueError, "complete inner validation"):
                    pilot.prepare(args)
                args.source_validation_method = "calibrated"
                pilot.prepare(args)
                root = scratch / "standard"
                configuration = pilot.verify(root)
                self.assertEqual(configuration["source_nuisance_method"], "standard")
                self.assertEqual(configuration["source_validation_method"], "complete")
                self.assertEqual(configuration["nlambda"], 100)
                provenance = json.loads((root / "provenance.json").read_text())
                self.assertIn(str(check / "SOURCE_STANDARD_CHECKS_PASSED"), provenance)
                args.baseline_methods = list(pilot.BASELINE_METHODS)
                with self.assertRaisesRegex(ValueError, "do not apply"):
                    pilot.prepare(args)

    def test_summary_keeps_source_and_target_ablation_cells_separate(self):
        import csv
        row = pilot.task_rows(["C1"], [2], [0], [1], "one_round", "treated_arm", [10])[0]
        manifest = [dict(row, task_id=str(index), source_nuisance_method=source,
                         target_nuisance_method=target)
                    for index, (source, target) in enumerate(
                        (("calibrated", "hou_calibrated"), ("standard", "hou_calibrated"),
                         ("calibrated", "lasso"), ("standard", "lasso")), start=1)]
        estimates = [dict(task, state="COMPLETE", recipe="score_derivative", method="target_only_ate",
                          estimate=.2, truth=.2, se=.1) for task in manifest]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            review.summarize_metrics(manifest, ["score_derivative"], estimates, root, ["target_only_ate"])
            with (root / "metrics.csv").open() as stream:
                metrics = list(csv.DictReader(stream))
            self.assertEqual(len(metrics), 4)
            self.assertTrue(all(row["n_planned"] == "1" and row["n_success"] == "1" for row in metrics))

    def test_source_count_profiles_distinguish_invalid_count_from_fraction(self):
        profiles = highdim.profiles("source_count")
        self.assertEqual(sum(len(highdim.rows_for_profile(p)) for p in profiles), 3000)
        self.assertEqual(sum(p["max_in_flight"] + p["baseline_max_in_flight"] for p in profiles), 432)
        half = [p for p in profiles if p["role"] == "fixed_invalid_fraction"]
        self.assertEqual({(p["sources"][0], p["n_deviated_sites"]) for p in half}, {(4, 2), (6, 3)})
        for profile in profiles:
            self.assertLessEqual(profile["n_deviated_sites"], min(profile["sources"]))
        row = pilot.task_rows(["C1"], [6], [1], [1], "one_round", "treated_arm", [100])[0]
        self.assertNotEqual(highdim.data_identity(row, {}),
                            highdim.data_identity(row, dict(n_deviated_sites=3)))
        self.assertEqual(highdim.data_identity(dict(row, rho=0), {}),
                         highdim.data_identity(dict(row, rho=0), dict(n_deviated_sites=3)))

    def test_source_count_and_deviation_count_are_explicit_cli_arguments(self):
        command = ["pilot", "prepare", "unused", "--check-root", "unused",
                   "--sources", "6", "--n-deviated-sites", "3"]
        with mock.patch.object(sys, "argv", command), mock.patch.object(pilot, "prepare") as prepare:
            pilot.main()
            self.assertEqual(prepare.call_args[0][0].sources, [6])
            self.assertEqual(prepare.call_args[0][0].n_deviated_sites, 3)

    def test_weak_fraction_followup_preserves_half_invalid_counts(self):
        profiles = highdim.profiles("weak_fraction")
        self.assertEqual(sum(len(highdim.rows_for_profile(p)) for p in profiles), 300)
        self.assertEqual(sum(p["max_in_flight"] + p["baseline_max_in_flight"] for p in profiles), 144)
        for profile in profiles:
            self.assertEqual(profile["rhos"], [.5])
            self.assertEqual(profile["dimensions"], [100])
            self.assertEqual(2 * profile["n_deviated_sites"], profile["sources"][0])
            self.assertTrue(profile["existing_parent"].startswith("source_count_validation_v1/"))

    def test_summary_does_not_pool_different_deviation_counts(self):
        row = pilot.task_rows(["C1"], [6], [1], [1], "one_round", "treated_arm", [100])[0]
        manifest = [dict(row, task_id=str(index), n_deviated_sites=str(count))
                    for index, count in enumerate((1, 3), start=1)]
        estimates = [dict(task, state="COMPLETE", recipe="score_derivative", method="target_only_ate",
                          estimate=.2, truth=.2, se=.1) for task in manifest]
        with tempfile.TemporaryDirectory() as directory:
            import csv
            root = Path(directory)
            review.summarize_metrics(manifest, ["score_derivative"], estimates, root,
                                     ["target_only_ate"])
            with (root / "metrics.csv").open() as stream:
                metrics = list(csv.DictReader(stream))
            self.assertEqual({row["n_deviated_sites"] for row in metrics}, {"1", "3"})
            self.assertTrue(all(row["n_planned"] == "1" and row["n_success"] == "1" for row in metrics))

    def test_c4_requires_explicit_selection_and_preserves_default_scenarios(self):
        for arguments, expected in (([], ["C1", "C2", "C3"]), (["--configs", "C4"], ["C4"])):
            command = ["pilot", "prepare", "unused", "--check-root", "unused"] + arguments
            with mock.patch.object(sys, "argv", command), mock.patch.object(pilot, "prepare") as prepare:
                pilot.main()
                self.assertEqual(prepare.call_args[0][0].configs, expected)

    def test_primary_highdim_seed_blocks_complete_200_without_repeating_old_seeds(self):
        cells = {(scenario, p, rho): set() for scenario in ("C1", "C2", "C3")
                 for p in (100, 200) for rho in (1, 2)}
        # Previously authorized p100 deviation block remains owned by its old jobs.
        for (scenario, p, rho), seeds in cells.items():
            if p == 100:
                seeds.update(range(1, 51))
        for profile in highdim.profiles():
            if profile["role"] not in ("primary", "deviation_validation"):
                continue
            for row in highdim.rows_for_profile(profile):
                if row["rho"] == .5:
                    continue
                seeds = cells[(row["config"], row["p"], row["rho"])]
                self.assertNotIn(row["sim_id"], seeds)
                seeds.add(row["sim_id"])
        self.assertTrue(all(seeds == set(range(1, 201)) for seeds in cells.values()))

    def test_baselines_are_shared_across_protocols_but_methods_are_distinct(self):
        row = pilot.task_rows(["C1"], [2], [1], [1], "one_round", "treated_arm", [200])[0]
        two_round = dict(row, protocol="two_round")
        baseline = dict(baseline_methods=list(pilot.BASELINE_METHODS))
        self.assertEqual(highdim.experiment_identity(row, baseline),
                         highdim.experiment_identity(two_round, baseline))
        self.assertNotEqual(highdim.experiment_identity(row, {}),
                            highdim.experiment_identity(two_round, {}))
        self.assertNotEqual(highdim.data_identity(row, {}),
                            highdim.data_identity(row, dict(dgp_control=dict(covariate_shift=2))))

    def test_new_submissions_do_not_inherit_controller_restart_counts(self):
        with mock.patch.dict(pilot.os.environ, SLURM_RESTART_COUNT="7", MODULEPATH="kept_modules"):
            environment = pilot.submission_environment()
            self.assertNotIn("SLURM_RESTART_COUNT", environment)
            self.assertEqual(environment["MODULEPATH"], "kept_modules")
            self.assertEqual(pilot.os.environ["SLURM_RESTART_COUNT"], "7")

    def test_expansion_can_use_a_repaired_controller_without_changing_workers(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "campaign"
            root.mkdir()
            workflow = root / "repaired_controller"
            campaign = dict(name="check", root=str(root), max_in_flight=4, gates=[],
                            controller_workflow=str(workflow))

            def submit(command, **kwargs):
                self.assertEqual(command[2], str(workflow / "submit_repeat_pilot.py"))
                self.assertEqual(command[3:], ["supervise", str(root), "--max-in-flight", "4"])
                (root / "controller_submission.json").write_text(
                    json.dumps(dict(state="submitted", stdout="123\n")))
                return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

            with mock.patch.object(pilot, "SCRATCH", Path(directory)), \
                 mock.patch.object(expansion.subprocess, "run", side_effect=submit) as submit_mock:
                self.assertTrue(expansion.advance(dict(campaigns=[campaign]), root / "status.json"))
                self.assertTrue(expansion.advance(dict(campaigns=[campaign]), root / "status.json"))
                self.assertEqual(submit_mock.call_count, 1)

    def test_expansion_gate_checks_execution_without_selecting_good_estimates(self):
        import csv
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "campaign"
            task = root / "tasks/1"
            task.mkdir(parents=True)
            (root / "manifest.csv").write_text("task_id,protocol\n1,one_round\n")
            gate = dict(root=str(root), tasks=[1])
            config = dict(recipes=["score_derivative"])
            with mock.patch.object(pilot, "SCRATCH", Path(directory)), \
                 mock.patch.object(pilot, "verify", return_value=config):
                self.assertFalse(expansion.check_gate(gate)[0])
                (task / "COMPLETE").touch()
                (task / "status.txt").write_text("COMPLETE\n")
                (task / "data_sha256.txt").write_text("a" * 64)
                rows = [dict(estimand_scope="tate", pilot_recipe="score_derivative", method=method,
                             estimate=.9, truth=.1, se=.01) for method in
                        ("one_round_crossfit_ate", "one_round_crossfit_ate_separate_arms",
                         "one_round_crossfit_ate_joint_tate", "target_only_ate")]
                def save_rows():
                    with (task / "results.csv").open("w", newline="") as stream:
                        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                        writer.writeheader()
                        writer.writerows(rows)
                save_rows()
                self.assertTrue(expansion.check_gate(gate)[0])
                with self.assertRaisesRegex(ValueError, "expected design"):
                    expansion.check_gate(dict(gate, expected_data_hashes={"1": "b" * 64}))
                self.assertTrue(expansion.check_gate(dict(gate, expected_data_hashes={"1": "a" * 64}))[0])
                rows[0]["se"] = 0
                save_rows()
                with self.assertRaisesRegex(ValueError, "Invalid completed"):
                    expansion.check_gate(gate)

    def test_staged_checkpoints_preserve_source_and_restart_progress(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "checkpoint_seeds/1/recipe/namespace"
            source.mkdir(parents=True)
            (source / "completed.rds").write_bytes(b"immutable fit")
            (source / "unfinished.tmp").write_bytes(b"partial fit")
            task = root / "tasks/1"
            task.mkdir(parents=True)
            pilot.restore_nuisance_checkpoints(root, 1, task)
            restored = task / "checkpoints/recipe/namespace"
            self.assertEqual((restored / "completed.rds").read_bytes(), b"immutable fit")
            self.assertFalse((restored / "unfinished.tmp").exists())
            (restored / "new_fit.rds").write_bytes(b"new progress")
            pilot.restore_nuisance_checkpoints(root, 1, task)
            self.assertEqual((restored / "new_fit.rds").read_bytes(), b"new progress")
            (restored / "completed.rds").unlink()
            self.assertEqual((source / "completed.rds").read_bytes(), b"immutable fit")

    def test_controller_wrapper_finds_workflow_after_slurm_copies_script(self):
        for explicit_workflow in (False, True):
            with self.subTest(explicit_workflow=explicit_workflow), tempfile.TemporaryDirectory() as directory:
                root = Path(directory) / "campaign with spaces"
                workflow = root / ("controller_workflow" if explicit_workflow else "workflow")
                workflow.mkdir(parents=True)
                (workflow / "supervise_repeat_pilot.py").write_text(
                    "import json, sys\nprint(json.dumps(sys.argv[1:]))\n")
                spool_script = Path(directory) / "slurm_script"
                spool_script.write_bytes(Path(__file__).with_name("supervise_repeat_pilot.sh").read_bytes())
                command = ["bash", str(spool_script), str(root), "6", "0.9"]
                if explicit_workflow:
                    command.append(str(workflow))
                result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True, timeout=10, check=True)
                self.assertEqual(json.loads(result.stdout),
                                 [str(root), "--max-in-flight", "6", "--hours", "0.9"])

    def test_first_wave_covers_every_dimension_before_new_seeds(self):
        rows = pilot.task_rows(["C1", "C2", "C3"], [2], [0], [1, 2],
            "one_round", "treated_arm", [10, 20, 50])
        self.assertEqual(len(rows), 18)
        self.assertEqual({row["sim_id"] for row in rows[:9]}, {1})
        rows = [{key: str(value) for key, value in row.items()} for row in rows]
        self.assertEqual(supervisor.submission_capacity(rows, {}, {}, 64, 1, set()), 9)
        jobs = {str(i): str(1000+i) for i in range(1, 10)}
        self.assertEqual(supervisor.submission_capacity(rows, jobs, {}, 64, 1, set()), 0)
        terminal = {job: {} for job in jobs.values()}
        self.assertEqual(supervisor.submission_capacity(rows, jobs, terminal, 4, 1, set(jobs)), 4)

    def test_scheduler_terminal_states_are_not_lost(self):
        response = subprocess.CompletedProcess([], 0,
            stdout="1001|COMPLETED|00:01:00|0:0\n1002|CANCELLED by 1|00:00:12|0:15\n", stderr="")
        with mock.patch.object(supervisor.subprocess, "run", return_value=response):
            states = supervisor.scheduler_states(["1001", "1002"], "2026-09-20")
        self.assertEqual(states["1002"]["scheduler_state"], "CANCELLED")

    def test_controller_waits_for_committed_output_after_accounting_completion(self):
        for publish_result in [True, False]:
            with self.subTest(publish_result=publish_result), tempfile.TemporaryDirectory() as directory, \
                 contextlib.redirect_stdout(io.StringIO()):
                root = Path(directory)
                (root / "configuration.json").write_text("{}")
                (root / "manifest.csv").write_text("task_id,sim_id\n1,1\n")
                (root / "tasks/1").mkdir(parents=True)
                clock = [0.0]
                states = {"12345": dict(job_id="12345", scheduler_state="COMPLETED", elapsed="00:01:00", exit_code="0:0")}
                def wait_for_visibility(interval):
                    clock[0] += 130
                    if publish_result: (root / "tasks/1/COMPLETE").write_text("done")
                def summarize(*args):
                    self.assertEqual((root / "tasks/1/COMPLETE").exists(), publish_result)
                    self.assertGreater(clock[0], 0)
                    return publish_result
                with mock.patch.object(pilot, "verify", return_value=dict(first_seed=1, checkpoint=True)), \
                     mock.patch.object(supervisor, "submitted_jobs", return_value={"1":"12345"}), \
                     mock.patch.object(supervisor, "scheduler_states", return_value=states), \
                     mock.patch.object(supervisor.time, "monotonic", side_effect=lambda:clock[0]), \
                     mock.patch.object(supervisor.time, "sleep", side_effect=wait_for_visibility), \
                     mock.patch.object(review, "summarize", side_effect=summarize):
                    self.assertEqual(supervisor.supervise(root,512,1,30), publish_result)

    def test_controller_recovers_original_job_id_without_resubmitting(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "configuration.json").write_text("{}")
            receipt = dict(state="submitting")
            path = root / "3.json"
            token = pilot.submission_comment(root, "3")
            live = subprocess.CompletedProcess([], 0, stdout="12345|" + token + "\n", stderr="")
            with mock.patch.object(supervisor.subprocess, "run", return_value=live):
                recovered = supervisor.recover_submission(root, path, receipt)
            self.assertEqual(recovered["stdout"], "12345\n")
            self.assertEqual(json.loads(path.read_text())["state"], "submitted")
            retry_token = token + ":startup_retry1"
            retry_live = subprocess.CompletedProcess([], 0,
                stdout="12345|" + token + "\n12346|" + retry_token + "\n", stderr="")
            with mock.patch.object(supervisor.subprocess, "run", return_value=retry_live):
                retried = supervisor.recover_submission(root, path,
                    dict(state="submitting", command=["sbatch", "--comment=" + retry_token]))
            self.assertEqual(retried["stdout"], "12346\n")
            empty = subprocess.CompletedProcess([], 0, stdout="", stderr="")
            with mock.patch.object(supervisor.subprocess, "run", return_value=empty):
                with self.assertRaisesRegex(ValueError, "Ambiguous"):
                    supervisor.recover_submission(root, path, dict(state="submitting"))

    def test_live_queue_overrides_stale_preempted_accounting(self):
        accounting = subprocess.CompletedProcess([], 0, stdout="1001|PREEMPTED|00:01:00|0:0\n", stderr="")
        live = subprocess.CompletedProcess([], 0, stdout="1001|PENDING\n", stderr="")
        with mock.patch.object(supervisor.subprocess, "run", side_effect=[accounting, live]):
            states = supervisor.scheduler_states(["1001"], "2026-09-20")
        self.assertEqual(states["1001"]["scheduler_state"], "PENDING")

    def test_scheduler_read_failures_do_not_resubmit_workers(self):
        timeout = subprocess.TimeoutExpired(["squeue"], 30)
        command_failure = subprocess.CalledProcessError(1, ["sacct"])
        for failures in ([timeout], [command_failure], [timeout] * 3):
            with self.subTest(failures=len(failures)), tempfile.TemporaryDirectory() as directory, \
                 contextlib.redirect_stdout(io.StringIO()):
                root = Path(directory)
                (root / "configuration.json").write_text("{}")
                (root / "manifest.csv").write_text("task_id,sim_id\n1,1\n")
                (root / "tasks/1").mkdir(parents=True)
                (root / "tasks/1/COMPLETE").touch()
                complete = {"1001": dict(job_id="1001", scheduler_state="COMPLETED",
                                        elapsed="00:01:00", exit_code="0:0")}
                def check_wait(interval):
                    status = json.loads((root / "supervisor_status.json").read_text())
                    self.assertEqual(status["state"], "WAITING_FOR_SCHEDULER")
                    self.assertEqual(status["next_batch"], 0)
                with mock.patch.object(pilot, "verify", return_value=dict(first_seed=1)), \
                     mock.patch.object(supervisor, "STOP_SIGNAL", None), \
                     mock.patch.object(supervisor, "submitted_jobs", return_value={"1": "1001"}), \
                     mock.patch.object(supervisor, "scheduler_states", side_effect=failures+[complete]), \
                     mock.patch.object(supervisor.time, "sleep", side_effect=check_wait), \
                     mock.patch.object(review, "summarize", return_value=True) as summarize, \
                     mock.patch.object(pilot, "submit") as submit:
                    result = supervisor.supervise(root, 32, 1, 30)
                submit.assert_not_called()
                status = json.loads((root / "supervisor_status.json").read_text())
                if len(failures) == 3:
                    self.assertIsNone(result)
                    self.assertEqual(status["state"], "CHECKPOINTED_FOR_REQUEUE")
                    self.assertEqual(status["scheduler_poll_failures"], 3)
                    summarize.assert_not_called()
                else:
                    self.assertTrue(result)
                    self.assertEqual(status["state"], "COMPLETE")

    def test_failed_live_read_discards_partial_accounting_snapshot(self):
        accounting = subprocess.CompletedProcess([], 0,
            stdout="1001|PREEMPTED|00:01:00|0:0\n", stderr="")
        with mock.patch.object(supervisor.subprocess, "run", side_effect=[
                accounting, subprocess.TimeoutExpired(["squeue"], 30)]):
            with self.assertRaises(subprocess.TimeoutExpired):
                supervisor.scheduler_states(["1001"], "2026-09-20")

    def test_requeue_recovers_only_the_same_job_and_skips_committed_results(self):
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            root = Path(directory) / "campaign"
            (root / "tasks").mkdir(parents=True)
            (root / "provenance.json").write_text("{}")
            (root / "manifest.csv").write_text("task_id\n1\n")
            args = SimpleNamespace(output=str(root), task_id=1)
            with mock.patch.object(pilot, "SCRATCH", Path(directory)), \
                 mock.patch.object(pilot, "verify", return_value=dict(checkpoint=True, max_restarts=5)), \
                 mock.patch.dict(pilot.os.environ, SLURM_JOB_ID="12345", SLURM_RESTART_COUNT="0"):
                pilot.validate_run(args)
                with self.assertRaises(FileExistsError):
                    pilot.validate_run(args)
                pilot.os.environ["SLURM_RESTART_COUNT"] = "1"
                pilot.validate_run(args)
                owner = json.loads((root / "tasks/1/owner.json").read_text())
                self.assertEqual(owner["restart"], 1)
                pilot.os.environ["SLURM_JOB_ID"] = "99999"
                pilot.os.environ["SLURM_RESTART_COUNT"] = "2"
                with self.assertRaisesRegex(ValueError, "identity mismatch"):
                    pilot.validate_run(args)
                pilot.os.environ["SLURM_JOB_ID"] = "12345"
                (root / "tasks/1/COMPLETE").write_text("done")
                pilot.validate_run(args)

    def test_large_concurrency_preserves_first_wave_and_counts_pending_jobs(self):
        rows = pilot.task_rows(["C1", "C2", "C3"], [2], [0], range(1, 201),
            "one_round", "treated_arm", [10, 20, 50])
        rows = [{key: str(value) for key, value in row.items()} for row in rows]
        first_jobs = {str(i): str(10000 + i) for i in range(1, 10)}
        completed = set(first_jobs)
        terminal = {job: {} for job in first_jobs.values()}
        self.assertEqual(supervisor.submission_capacity(rows, {}, {}, 512, 1, set()), 9)
        self.assertEqual(supervisor.submission_capacity(rows, first_jobs, terminal,
            512, 1, completed), 512)
        jobs = {str(i): str(10000 + i) for i in range(1, 522)}
        self.assertEqual(supervisor.submission_capacity(rows, jobs, terminal,
            512, 1, completed), 0)
        self.assertEqual(supervisor.submission_capacity(rows, jobs, terminal,
            2000, 1, completed), 1279)

    def test_monte_carlo_metrics_retain_the_planned_denominator(self):
        rows = pilot.task_rows(["C1"], [2], [0], [1, 2, 3], "one_round", "treated_arm", [10])
        rows = [{key: str(value) for key, value in row.items()} for row in rows]
        estimates = []
        for task, error in zip(rows[:2], [-.1, .1]):
            for method, shift in [("target_only_ate", 0), ("one_round_crossfit_ate", .05)]:
                estimates.append(dict(task, recipe="score_derivative", state="COMPLETE",
                    method=method, estimate=.2+error+shift, se=.1, truth=.2))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            review.summarize_metrics(rows, ["score_derivative"], estimates, root)
            import csv
            with (root / "metrics.csv").open() as stream:
                metrics = list(csv.DictReader(stream))
            reference = next(row for row in metrics if row["method"] == "target_only_ate")
            self.assertEqual(reference["n_planned"], "3")
            self.assertEqual(reference["n_success"], "2")
            self.assertAlmostEqual(float(reference["rmse"]), .1)
            self.assertGreater(float(reference["coverage_lower"]), 0)
            self.assertLess(float(reference["coverage_lower"]), 1)

    def test_baseline_metrics_use_requested_methods_and_preserve_missing_repeats(self):
        rows = pilot.task_rows(["C1"], [2], [0], [1, 2], "one_round", "treated_arm", [100])
        rows = [{key: str(value) for key, value in row.items()} for row in rows]
        estimates = [dict(rows[0], recipe="baselines", state="COMPLETE", method=method,
                          estimate=.21, se=.1, truth=.2)
                     for method in ("target_only_ate", "pooled_dr_ate")]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            review.summarize_metrics(rows, ["baselines"], estimates, root,
                                     ["target_only_ate", "pooled_dr_ate"])
            import csv
            with (root / "metrics.csv").open() as stream:
                metrics = list(csv.DictReader(stream))
            self.assertEqual({row["method"] for row in metrics},
                             {"target_only_ate", "pooled_dr_ate"})
            self.assertTrue(all(row["n_success"] == "1" and row["n_unavailable"] == "1" for row in metrics))

    def test_partition_override_changes_only_submission_routing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            configuration = dict(partitions="preempt", checkpoint=True)
            self.assertEqual(pilot.submission_partitions(root, configuration), "preempt")
            pilot.replace_json(root / "routing.json", dict(partitions=pilot.DEFAULT_PARTITIONS))
            self.assertEqual(pilot.submission_partitions(root, configuration), pilot.DEFAULT_PARTITIONS)
            self.assertEqual(configuration["partitions"], "preempt")
            with self.assertRaisesRegex(ValueError, "checkpoint/requeue"):
                pilot.submission_partitions(root, dict(partitions="agsmall", checkpoint=False))
            pilot.replace_json(root / "routing.json", dict(partitions="preempt", nuisance_tol=1))
            with self.assertRaisesRegex(ValueError, "only the partitions"):
                pilot.submission_partitions(root, configuration)
            pilot.replace_json(root / "routing.json", dict(partitions="preempt", exclude_nodes="acl22,acl74"))
            self.assertEqual(pilot.submission_partitions(root, configuration), "preempt")
            self.assertEqual(pilot.submission_excluded_nodes(root, configuration), "acl22,acl74")

    def test_baseline_preparation_dispatches_the_baseline_worker(self):
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            scratch = Path(directory)
            check = scratch / "check"
            (check / "Rlib/RoCE").mkdir(parents=True)
            (check / "Rlib/RoCE/DESCRIPTION").write_text("test package")
            (check / "TESTS_PASSED").touch()
            (check / "source_manifest.json").write_text("[]")
            root = scratch / "baselines"
            args = SimpleNamespace(output=str(root), check_root=str(check), configs=["C1"],
                sources=[6], n_deviated_sites=3, rhos=[1], repeats=[1], protocol="one_round", deviation_mechanism="treated_arm",
                recipes=None, cpus=2, memory="4G", time_limit="00:30:00", partitions="preempt",
                account="hou00123", shared_cache=False, baseline_methods=list(pilot.BASELINE_METHODS))
            with mock.patch.object(pilot, "SCRATCH", scratch):
                pilot.prepare(args)
                config = pilot.verify(root)
                self.assertEqual(config["recipes"], ["baselines"])
                self.assertEqual(config["worker_script"], "run_baseline_repeat.R")
                self.assertEqual(config["n_deviated_sites"], 3)
                self.assertTrue((root / "workflow/run_baseline_repeat.R").is_file())
                response = subprocess.CompletedProcess([], 0, stdout="12345\n", stderr="")
                with mock.patch.object(pilot.subprocess, "run", return_value=response) as sbatch:
                    pilot.submit(SimpleNamespace(output=str(root), max_jobs=1, test_only=False))
                    self.assertEqual(sbatch.call_args[0][0][-1], "run_baseline_repeat.R")
                args.target_nuisance_method = "lasso"
                with self.assertRaisesRegex(ValueError, "do not apply"):
                    pilot.prepare(args)
                args.target_nuisance_method = "hou_calibrated"
                args.n_deviated_sites = 7
                with self.assertRaisesRegex(ValueError, "every requested source count"):
                    pilot.prepare(args)

    def test_manifest_keeps_each_scenario_repeat_unique(self):
        rows = pilot.task_rows(["C1", "C3"], [2], [0, 2.5], [1, 2], "one_round", "treated_arm")
        self.assertEqual(len(rows), 8)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            pilot.task_rows(["C1", "C1"], [2], [0], [1], "one_round", "treated_arm")
        with self.assertRaisesRegex(ValueError, "C1 only"):
            pilot.task_rows(["C3"], [2], [0], [1], "one_round", "both_arms")

    def test_submissions_and_workers_cannot_duplicate_a_repeat(self):
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            scratch = Path(directory)
            check = scratch / "check"
            (check / "Rlib" / "RoCE").mkdir(parents=True)
            (check / "Rlib" / "RoCE" / "DESCRIPTION").write_text("test package")
            (check / "TESTS_PASSED").touch()
            (check / "source_manifest.json").write_text("[]")
            root = scratch / "pilot"
            configuration = SimpleNamespace(output=str(root), check_root=str(check),
                configs=["C1"], sources=[2], rhos=[0], repeats=[1, 2], protocol="one_round",
                deviation_mechanism="treated_arm", recipes=["legacy", "score_derivative"],
                cpus=4, memory="16G", time_limit="06:00:00",
                partitions=pilot.DEFAULT_PARTITIONS, exclude_nodes="acl03,acl74",
                aggregation_cutoff=2.0, account="hou00123", shared_cache=False)
            with mock.patch.object(pilot, "SCRATCH", scratch):
                pilot.prepare(configuration)
                self.assertEqual(pilot.verify(root)["aggregation_lambda"], 0.5)
                request = SimpleNamespace(output=str(root), max_jobs=2, test_only=False)
                response = subprocess.CompletedProcess([], 0, stdout="12345\n", stderr="")
                with mock.patch.object(pilot.subprocess, "run", return_value=response) as sbatch:
                    pilot.submit(request)
                    self.assertEqual(sbatch.call_count, 2)
                    for call in sbatch.call_args_list:
                        command = call[0][0]
                        self.assertIn("--partition=" + pilot.DEFAULT_PARTITIONS, command)
                        self.assertIn("--exclude=acl03,acl74", command)
                        self.assertFalse(any(option.startswith("--array") for option in command))
                    pilot.submit(request)
                    self.assertEqual(sbatch.call_count, 2)
                    pilot.submit_supervisor(SimpleNamespace(output=str(root), max_in_flight=512))
                    self.assertIn("--exclude=acl03,acl74", sbatch.call_args[0][0])
                    self.assertEqual(sbatch.call_args[0][0][-1], str(Path(pilot.__file__).resolve().parent))
                task = SimpleNamespace(output=str(root), task_id=1)
                pilot.validate_run(task)
                with self.assertRaises(FileExistsError):
                    pilot.validate_run(task)
                (check / "Rlib" / "RoCE" / "DESCRIPTION").write_text("changed")
                with self.assertRaisesRegex(ValueError, "Changed input"):
                    pilot.verify(root)


if __name__ == "__main__":
    unittest.main()
