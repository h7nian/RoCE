#!/usr/bin/env python3
"""Population checks of calibration against the derivative of a clipped score.

No RoCE implementation is imported. These examples diagnose orthogonality,
not finite-sample coverage or the FACE DGP's performance.
"""

import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.optimize import root
from scipy.special import expit


def solve_equation(score, jacobian, initial):
    fitted = root(score, np.asarray(initial, dtype=float), jac=jacobian, tol=1e-12)
    residual = np.max(np.abs(score(fitted.x)))
    if not np.all(np.isfinite(fitted.x)) or residual > 1e-10:
        raise RuntimeError("Population equation failed: " + fitted.message)
    return fitted.x


def fit_outcome(design, mass, outcome, weights):
    effective_mass = mass * weights
    return solve_equation(
        lambda coefficient: design.T @ (effective_mass * (expit(design @ coefficient) - outcome)),
        lambda coefficient: design.T @ ((effective_mass * expit(design @ coefficient) *
            (1 - expit(design @ coefficient)))[:, None] * design), np.zeros(design.shape[1]))


def fit_tilt(design, source_mass, linear_moment, radius, initial=None):
    def weight(coefficient):
        return np.exp(-np.clip(design @ coefficient, -radius, radius))
    return solve_equation(
        lambda coefficient: linear_moment - design.T @ (source_mass * weight(coefficient)),
        lambda coefficient: design.T @ ((source_mass * weight(coefficient) *
            (np.abs(design @ coefficient) < radius))[:, None] * design),
        np.zeros(design.shape[1]) if initial is None else initial)


def mean_derivative(predictor, radius=np.inf):
    mean = expit(np.clip(predictor, -radius, radius))
    return mean * (1 - mean)


def check_case(kind, correct, design, target_mass, source_mass, probability,
               truth, outcome_initial, weight_initial, radius):
    rows = []
    initial_predictor = design @ weight_initial
    for recipe in ("legacy_clipped_plugins", "score_derivatives"):
        matched = recipe == "score_derivatives"
        derivative = mean_derivative(design @ outcome_initial, np.inf if matched else radius)
        if kind == "target":
            moment = design.T @ (target_mass * (1 - probability) * derivative)
            weight_mass = target_mass * probability * derivative
        else:
            moment = design.T @ (target_mass * derivative)
            weight_mass = source_mass * probability * derivative
        weight_final = fit_tilt(design, weight_mass, moment, radius, weight_initial)
        outcome_weight = probability * np.exp(-np.clip(initial_predictor, -radius, radius))
        if matched:
            outcome_weight *= np.abs(initial_predictor) < radius
        outcome_final = fit_outcome(design, source_mass, truth, outcome_weight)
        prediction = expit(design @ outcome_final)
        tilt = np.exp(-np.clip(design @ weight_final, -radius, radius))
        tilt_derivative = tilt * (np.abs(design @ weight_final) < radius)
        if kind == "target":
            inverse_probability = 1 + tilt
            estimate = np.sum(target_mass * (prediction + probability * inverse_probability * (truth - prediction)))
            sensitivity_outcome = design.T @ (target_mass * (1 - probability * inverse_probability) *
                                               mean_derivative(design @ outcome_final))
        else:
            estimate = np.sum(target_mass * prediction) + np.sum(source_mass * probability * tilt * (truth - prediction))
            sensitivity_outcome = design.T @ ((target_mass - source_mass * probability * tilt) *
                                               mean_derivative(design @ outcome_final))
        sensitivity_weight = -design.T @ (source_mass * probability * tilt_derivative * (truth - prediction))
        expected = np.sum(target_mass * truth)
        row = dict(site=kind, correct_branch=correct, recipe=recipe,
                   estimate_error=float(estimate - expected),
                   outcome_sensitivity=float(np.max(np.abs(sensitivity_outcome))),
                   weight_sensitivity=float(np.max(np.abs(sensitivity_weight))),
                   initial_final_weight_distance=float(np.max(np.abs(weight_final - weight_initial))),
                   initial_final_outcome_distance=float(np.max(np.abs(outcome_final - outcome_initial))))
        if abs(row["estimate_error"]) > 1e-10:
            raise AssertionError("The population point estimate should be correct in this example")
        maximum = max(row["outcome_sensitivity"], row["weight_sensitivity"])
        if matched and maximum > 1e-10:
            raise AssertionError("Score-derivative calibration failed its population sensitivity check")
        if not matched and maximum < 1e-7:
            raise AssertionError("Example did not expose the legacy sensitivity mismatch")
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    if Path("/scratch.global/zhan9381/FACE-HD") not in output.parents or output.exists():
        parser.error("Choose a new directory below the designated FACE-HD scratch root")
    x = np.array([-3., -1., 0., 1., 3.])
    design = np.column_stack((np.ones(len(x)), x))
    mass = np.array([.10, .15, .15, .20, .40])
    outcome_coefficient = np.array([1.2, 2.])
    outcome_truth = expit(design @ outcome_coefficient)
    probability = np.clip(expit(-.2 + .25 * x + .35 * x**2), .1, .9)
    rows = []
    weight_initial = fit_tilt(design, mass * probability, design.T @ (mass * (1 - probability)), 5)
    rows += check_case("target", "outcome", design, mass, mass, probability,
                       outcome_truth, outcome_coefficient, weight_initial, 5)
    target_mass = np.array([.20, .20, .30, .20, .10])
    weight_initial = fit_tilt(design, mass * probability, design.T @ target_mass, 5)
    rows += check_case("source", "outcome", design, target_mass, mass, probability,
                       outcome_truth, outcome_coefficient, weight_initial, 5)
    # Correct bounded propensity. Supplying its true initial coefficient makes
    # the legacy OR-loss mismatch visible even under ideal initialization.
    radius = np.log(9.)
    weight_truth = np.array([.2, 1.3])
    probability = expit(np.clip(design @ weight_truth, -radius, radius))
    outcome_truth = expit(.1 + .25 * x + .3 * x**2)
    outcome_initial = fit_outcome(design, mass, outcome_truth, probability)
    rows += check_case("target", "propensity", design, mass, mass, probability,
                       outcome_truth, outcome_initial, weight_truth, radius)
    # Correct clipped merged transport weight. Choose a constant treatment
    # probability so that the target empirical measure is normalized exactly.
    radius = np.log(4.)
    weight_truth = np.array([-.4, .9])
    tilt = np.exp(-np.clip(design @ weight_truth, -radius, radius))
    probability = np.repeat(1 / np.sum(mass * tilt), len(x))
    target_mass = mass * probability * tilt
    outcome_initial = fit_outcome(design, mass, outcome_truth, probability)
    rows += check_case("source", "merged_weight", design, target_mass, mass, probability,
                       outcome_truth, outcome_initial, weight_truth, radius)
    output.mkdir(parents=True)
    with (output / "population_sensitivities.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / "scope.json").write_text(json.dumps(dict(
        scope="Independent finite-support population equations; not a FACE Monte Carlo experiment",
        matched_recipe="Unclipped initial OR derivative in weight loss; derivative of clipped initial weight in OR loss",
        caveats=["Rates, regularization, selection inference and fold dependence are not proved here",
                 "The correct clipped-PS example uses the matching clipping radius and an ideal consistent initial PS"]), indent=2) + "\n")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
