#!/usr/bin/env python3
"""Independently check the arm-specific calibration mechanism used by Hou et al.

This checks unpenalized population equations on finite support, where truncation
is inactive. It is not an implementation of SMMAL, a high-dimensional proof, or
a FACE performance study. theta models logit P(A=arm), beta models that arm's OR.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from population_checks import central_gradient, json_value, logistic_projection, solve_score


def calibrated_propensity(design, mass, arm_probability, initial_outcome):
    outcome_mean = expit(design @ initial_outcome)
    derivative = outcome_mean * (1 - outcome_mean)

    def score(coefficient):
        inverse_odds = np.exp(-design @ coefficient)
        return design.T @ (mass * derivative *
                           (1 - arm_probability - arm_probability * inverse_odds))

    def jacobian(coefficient):
        curvature = mass * derivative * arm_probability * np.exp(-design @ coefficient)
        return design.T @ (curvature[:, None] * design)

    return solve_score(score, jacobian, np.zeros(design.shape[1]))


def check_case(case, arm):
    covariate = np.arange(-2., 3.)
    design = np.column_stack([np.ones(5), covariate])
    mass = np.ones(5) / 5
    true_treated_probability = (np.array([.15, .8, .25, .9, .35]) if case == "outcome_correct"
                                else expit(design @ np.array([.2, .45])))
    arm_probability = true_treated_probability if arm == 1 else 1 - true_treated_probability
    true_outcome_coefficient = np.array([-.4, .7]) if arm == 1 else np.array([-.7, -.3])
    outcome_predictor = design @ true_outcome_coefficient
    if case == "propensity_correct":
        outcome_predictor += .5 * covariate ** 2
    true_outcome_mean = expit(outcome_predictor)
    theta_initial, _ = logistic_projection(design, mass, arm_probability)
    beta_initial, _ = logistic_projection(design, mass * arm_probability, true_outcome_mean)
    theta_final, propensity_score_residual = calibrated_propensity(
        design, mass, arm_probability, beta_initial)
    beta_final, outcome_score_residual = logistic_projection(
        design, mass * arm_probability * np.exp(-design @ theta_initial), true_outcome_mean)
    fitted_probability = expit(design @ theta_final)
    fitted_mean = expit(design @ beta_final)
    outcome_gradient = design.T @ (mass * (1 - arm_probability / fitted_probability) *
                                   fitted_mean * (1 - fitted_mean))
    propensity_gradient = -design.T @ (mass * arm_probability * np.exp(-design @ theta_final) *
                                       (true_outcome_mean - fitted_mean))

    def functional(parameter):
        probability = expit(design @ parameter[:2])
        mean = expit(design @ parameter[2:])
        return float(mass @ (mean + arm_probability / probability * (true_outcome_mean - mean)))

    parameter = np.r_[theta_final, beta_final]
    analytic_gradient = np.r_[propensity_gradient, outcome_gradient]
    finite_gradient = central_gradient(functional, parameter)
    truth = float(mass @ true_outcome_mean)
    estimate = functional(parameter)
    assert abs(estimate - truth) < 1e-10
    assert np.max(np.abs(analytic_gradient)) < 1e-10
    assert np.max(np.abs(finite_gradient - analytic_gradient)) < 1e-9
    if case == "outcome_correct":
        assert np.linalg.norm(theta_initial - theta_final) > 1e-3
        assert np.max(np.abs(beta_initial - beta_final)) < 1e-9
    elif case == "propensity_correct":
        assert np.linalg.norm(beta_initial - beta_final) > 1e-3
        assert np.max(np.abs(theta_initial - theta_final)) < 1e-9
    else:
        assert np.max(np.abs(theta_initial - theta_final)) < 1e-9
        assert np.max(np.abs(beta_initial - beta_final)) < 1e-9
    return dict(case=case, arm=arm, estimate=estimate, truth=truth,
                max_gradient=float(np.max(np.abs(analytic_gradient))),
                max_finite_difference_error=float(np.max(np.abs(finite_gradient - analytic_gradient))),
                theta_initial_final_difference=float(np.linalg.norm(theta_initial - theta_final)),
                beta_initial_final_difference=float(np.linalg.norm(beta_initial - beta_final)),
                theta_initial=theta_initial, theta_final=theta_final,
                beta_initial=beta_initial, beta_final=beta_final,
                propensity_score_residual=propensity_score_residual,
                outcome_score_residual=outcome_score_residual,
                covariate=covariate, mass=mass, arm_probability=arm_probability,
                true_outcome_mean=true_outcome_mean)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    root = Path("/scratch.global/zhan9381/FACE-HD")
    if root not in output.parents or output.exists():
        parser.error("Choose a new output directory under /scratch.global/zhan9381/FACE-HD")
    results = [check_case(case, arm) for case in
               ["both_correct", "propensity_correct", "outcome_correct"] for arm in [0, 1]]
    output.mkdir(parents=True)
    source = Path(__file__).resolve()
    payload = dict(scope=__doc__, results=results, source_sha256={
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in [source, source.with_name("population_checks.py")]})
    (output / "population_checks.json").write_text(
        json.dumps(payload, indent=2, default=json_value, allow_nan=False) + "\n")
    fields = ["case", "arm", "estimate", "truth", "max_gradient", "max_finite_difference_error",
              "theta_initial_final_difference", "beta_initial_final_difference"]
    with (output / "summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({field: result[field] for field in fields} for result in results)
    print("Verified six population calibration cases; results: " + str(output))


if __name__ == "__main__":
    main()
