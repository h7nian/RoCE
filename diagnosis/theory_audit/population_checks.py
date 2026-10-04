#!/usr/bin/env python3
"""Independent population checks for the RoCE theory audit.

These finite-support calculations do not import RoCE or estimate performance of
the FACE simulation. A successful check can reproduce a counterexample; it does
not certify the production estimator. All output must go to scratch.global.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import root
from scipy.special import expit
from scipy.stats import norm


def solve_score(score, jacobian, initial):
    """Solve a small population score, checking its residual independently."""
    fit = root(score, np.asarray(initial, dtype=float), jac=jacobian, tol=1e-11)
    residual = float(np.max(np.abs(score(fit.x))))
    if not np.all(np.isfinite(fit.x)) or residual > 1e-10:
        raise RuntimeError("Population score did not solve: " + str(fit.message))
    return fit.x, residual


def logistic_projection(design, mass, response_mean):
    def score(coefficient):
        return design.T @ (mass * (expit(design @ coefficient) - response_mean))

    def jacobian(coefficient):
        mean = expit(design @ coefficient)
        return design.T @ ((mass * mean * (1 - mean))[:, None] * design)

    return solve_score(score, jacobian, np.zeros(design.shape[1]))


def tilt_projection(design, target_mass, source_arm_mass, derivative):
    target_moment = design.T @ (target_mass * derivative)

    def score(coefficient):
        source_weight = source_arm_mass * derivative * np.exp(-design @ coefficient)
        return target_moment - design.T @ source_weight

    def jacobian(coefficient):
        source_weight = source_arm_mass * derivative * np.exp(-design @ coefficient)
        return design.T @ (source_weight[:, None] * design)

    return solve_score(score, jacobian, np.zeros(design.shape[1]))


def central_gradient(function, parameter, step=1e-5):
    directions = np.eye(len(parameter)) * step
    return np.array([(function(parameter + d) - function(parameter - d)) / (2 * step)
                     for d in directions])


def finite_balance_check():
    """All linear features balance, but the transported mean is inconsistent."""
    x = np.array([-1., 0., 1.])
    design = np.column_stack([np.ones(3), x])
    mass = np.ones(3) / 3
    propensity = np.array([.25, .75, .25])
    outcome_mean = np.array([.8, .2, .8])
    alpha, alpha_residual = logistic_projection(design, mass * propensity, outcome_mean)
    gamma, gamma_residual = tilt_projection(design, mass, mass * propensity, np.ones(3))
    prediction = expit(design @ alpha)
    weighted_mass = mass * propensity * np.exp(-design @ gamma)
    estimate = float(mass @ prediction + weighted_mass @ (outcome_mean - prediction))
    truth = float(mass @ outcome_mean)
    balance = design.T @ (weighted_mass - mass)
    calibrated_gamma, _ = tilt_projection(
        design, mass, mass * propensity, prediction * (1 - prediction))
    calibrated_alpha, _ = logistic_projection(design, weighted_mass, outcome_mean)
    assert np.max(np.abs(balance)) < 1e-10
    assert np.max(np.abs(calibrated_alpha - alpha)) < 1e-9
    assert np.max(np.abs(calibrated_gamma - gamma)) < 1e-9
    assert abs(estimate - .44) < 1e-10 and abs(truth - .6) < 1e-10
    return dict(x=x, mass=mass, propensity=propensity, outcome_mean=outcome_mean,
                alpha=alpha, gamma=gamma, balance=balance, estimate=estimate,
                truth=truth, bias=estimate - truth,
                score_residual=max(alpha_residual, gamma_residual))


def alignment_and_remainder_check():
    """Correct logistic OR: tilt limits differ, yet both final scores vanish.

    This refutes automatic initialization alignment and a pure product bound,
    not double robustness of the actual update on the correct-OR branch.
    """
    x = np.arange(-2., 3.)
    design = np.column_stack([np.ones(5), x])
    target_mass = np.array([.1, .2, .4, .2, .1])
    source_mass = np.array([.25, .15, .2, .1, .3])
    propensity = np.array([.2, .6, .3, .8, .4])
    alpha = np.array([-.4, .7])
    outcome_mean = expit(design @ alpha)
    derivative = outcome_mean * (1 - outcome_mean)
    gamma_initial, initial_residual = tilt_projection(
        design, target_mass, source_mass * propensity, np.ones(5))
    gamma_final, final_residual = tilt_projection(
        design, target_mass, source_mass * propensity, derivative)
    weighted_mass = source_mass * propensity * np.exp(-design @ gamma_final)
    alpha_final, _ = logistic_projection(
        design, source_mass * propensity * np.exp(-design @ gamma_initial), outcome_mean)

    def functional(coefficient):
        prediction = expit(design @ coefficient)
        return float(target_mass @ prediction + weighted_mass @ (outcome_mean - prediction))

    gradient = design.T @ ((target_mass - weighted_mass) * derivative)
    finite_gradient = central_gradient(functional, alpha)
    second_derivative = float(np.sum((target_mass - weighted_mass) * derivative *
                                     (1 - 2 * outcome_mean)))
    remainders = []
    for step in [.1, .01, .001]:
        remainder = functional(alpha + np.array([step, 0])) - functional(alpha)
        remainders.append(dict(step=step, remainder=remainder,
                               remainder_over_step_squared=remainder / step ** 2,
                               product_of_coefficient_errors=0.0))
    assert np.linalg.norm(gamma_initial - gamma_final) > .05
    assert np.max(np.abs(alpha_final - alpha)) < 1e-9
    assert np.max(np.abs(gradient)) < 1e-10
    assert np.max(np.abs(finite_gradient - gradient)) < 1e-9
    assert abs(remainders[-1]["remainder"]) > 1e-10
    assert abs(remainders[-1]["remainder_over_step_squared"] - second_derivative / 2) < 1e-5
    return dict(x=x, target_mass=target_mass, source_mass=source_mass,
                propensity=propensity, alpha=alpha, gamma_initial=gamma_initial,
                gamma_final=gamma_final, gradient=gradient,
                finite_difference_gradient=finite_gradient,
                intercept_second_derivative=second_derivative,
                score_residual=max(initial_residual, final_residual), remainders=remainders)


def target_anchor_check():
    """Exact MLE limit for a correct logistic OR and misspecified logistic PS.

    Refit both nuisance models under case-mass perturbations. This validates the
    actual low-dimensional functional derivative, not a lasso correction.
    """
    x = np.arange(-2., 3.)
    design = np.column_stack([np.ones(5), x])
    mass = np.ones(5) / 5
    propensity = np.array([.15, .8, .25, .9, .35])
    outcome_mean = {0: expit(design @ np.array([-.8, .4])),
                    1: expit(design @ np.array([-.3, .7]))}
    gamma, _ = logistic_projection(design, mass, propensity)
    fitted_propensity = expit(design @ gamma)
    outcome_sensitivity = {}
    outcome_correction = {}
    for arm in [0, 1]:
        p_arm = propensity if arm == 1 else 1 - propensity
        q_arm = fitted_propensity if arm == 1 else 1 - fitted_propensity
        derivative = outcome_mean[arm] * (1 - outcome_mean[arm])
        sensitivity = design.T @ (mass * (1 - p_arm / q_arm) * derivative)
        information = design.T @ ((mass * p_arm * derivative)[:, None] * design)
        outcome_sensitivity[arm] = sensitivity
        outcome_correction[arm] = np.linalg.solve(information, sensitivity)
    cells = [(i, arm, y) for i in range(5) for arm in [0, 1] for y in [0, 1]]
    cell_mass = np.array([mass[i] * (propensity[i] if arm else 1 - propensity[i]) *
                          (outcome_mean[arm][i] if y else 1 - outcome_mean[arm][i])
                          for i, arm, y in cells])
    cell_design = np.array([design[i] for i, _, _ in cells])
    treatment = np.array([arm for _, arm, _ in cells])
    outcome = np.array([y for _, _, y in cells])
    truth = float(mass @ (outcome_mean[1] - outcome_mean[0]))
    oracle_if = []
    corrected_if = []
    for i, arm, y in cells:
        sign = 1 if arm else -1
        q_arm = fitted_propensity[i] if arm else 1 - fitted_propensity[i]
        residual = y - outcome_mean[arm][i]
        direct = outcome_mean[1][i] - outcome_mean[0][i] - truth + sign * residual / q_arm
        oracle_if.append(direct)
        corrected_if.append(direct + sign * (design[i] @ outcome_correction[arm]) * residual)
    oracle_if, corrected_if = np.array(oracle_if), np.array(corrected_if)

    def refitted_functional(weights):
        weights = weights / np.sum(weights)
        ps, _ = logistic_projection(cell_design, weights, treatment)
        q = expit(cell_design @ ps)
        value = 0.
        for arm in [0, 1]:
            mask = treatment == arm
            coefficient, _ = logistic_projection(cell_design, weights * mask, outcome)
            prediction = expit(cell_design @ coefficient)
            q_arm = q if arm else 1 - q
            value += (1 if arm else -1) * (weights @ (
                prediction + mask * (outcome - prediction) / q_arm))
        return float(value)

    step = 1e-5
    directions = np.eye(len(cell_mass)) - cell_mass[None, :]
    finite_if = np.array([(refitted_functional(cell_mass + step * direction) -
                           refitted_functional(cell_mass - step * direction)) / (2 * step)
                          for direction in directions])
    correction_error = float(np.max(np.abs(finite_if - corrected_if)))
    oracle_error = float(np.max(np.abs(finite_if - oracle_if)))
    assert abs(refitted_functional(cell_mass) - truth) < 1e-10
    assert correction_error < 1e-7 and oracle_error > .01
    assert abs(cell_mass @ corrected_if) < 1e-10
    return dict(x=x, mass=mass, propensity=propensity, fitted_propensity=fitted_propensity,
                outcome_mean=outcome_mean, outcome_sensitivity=outcome_sensitivity,
                truth=truth, oracle_if_variance=float(cell_mass @ oracle_if ** 2),
                refitted_if_variance=float(cell_mass @ corrected_if ** 2),
                oracle_derivative_error=oracle_error, corrected_derivative_error=correction_error,
                oracle_if=oracle_if, corrected_if=corrected_if,
                finite_difference_if=finite_if, cell_mass=cell_mass)


def clipping_check():
    """Compare the convex continuation, displayed exponential and clipped score."""
    radius = 5.

    def continued_loss(value):
        clipped = np.clip(value, -radius, radius)
        return math.exp(-clipped) * (1 - (value - clipped))

    rows = []
    for value in [-6., -4., 0., 4., 6.]:
        step = 1e-5
        weight = math.exp(-np.clip(value, -radius, radius))
        derivative = (continued_loss(value + step) - continued_loss(value - step)) / (2 * step)
        clipped_exp_derivative = (math.exp(-np.clip(value + step, -radius, radius)) -
                                  math.exp(-np.clip(value - step, -radius, radius))) / (2 * step)
        assert abs(derivative + weight) < 1e-6
        rows.append(dict(predictor=value, weight=weight, continued_loss=continued_loss(value),
                         continued_loss_derivative=derivative,
                         clipped_exponential_derivative=clipped_exp_derivative,
                         untruncated_exponential_derivative=-math.exp(-value)))
    return rows


def fixed_cutoff_check():
    """An exact Gaussian null example: weight randomness persists as n grows."""
    def weight(z):
        return max(0., .5 - max(abs(z) - 1., 0.) / 8.)

    squared_error = 2 * quad(lambda z: (weight(z) - .5) ** 2 * norm.pdf(z), 1, 12,
                            points=[5.], epsabs=1e-12)[0]
    assert squared_error > .001
    return dict(oracle_weight=.5, penalty_active_probability=float(2 * norm.sf(1)),
                mean_squared_distance_from_oracle=squared_error,
                sample_size_dependence="none: Z has N(0,1) law for every sample size")


def dimension_check():
    """Coordinatewise bounds do not bound the feature Gram operator norm."""
    rows = []
    # X_j = (Z + epsilon_j)/sqrt(2), independent Rademacher Z and epsilons.
    # Every X_j is bounded by sqrt(2); Gram = (I + 11')/2 has lambda_min=1/2.
    for dimension in [4, 20, 100, 200]:
        coefficient_error_squared = 1 / dimension  # delta_j = 1/dimension
        prediction_error_squared = .5 + .5 / dimension
        rows.append(dict(dimension=dimension, gram_min_eigenvalue=.5,
                         gram_max_eigenvalue=(dimension + 1) / 2,
                         coefficient_error_squared=coefficient_error_squared,
                         prediction_error_squared=prediction_error_squared,
                         required_lipschitz_constant=prediction_error_squared / coefficient_error_squared))
    return rows


def multisite_covariance_check():
    """Check the shared-target covariance and common-weight arm contrast."""
    target_treated = np.array([[.6, .8, .7], [.3, .4, .5], [.9, .6, .8], [.4, .5, .2]])
    target_control = np.array([[.2, .3, .1], [.4, .1, .2], [.5, .2, .4], [.1, .2, .3]])
    source_arms = [np.array([[.1, -.2], [-.3, .1], [.2, .2]]),
                   np.array([[.3, .1], [-.1, .3], [.4, -.2], [-.2, .2]])]
    source_weights = np.array([.2, .3])
    weights = np.r_[1 - source_weights.sum(), source_weights]
    target_contrast = target_treated - target_control
    sample_sizes = np.array([1000., 1000.])
    target_covariance = np.cov(target_contrast, rowvar=False, bias=True) / 1000
    source_variances = np.array([np.var(arms[:, 0] - arms[:, 1]) for arms in source_arms])
    estimator_covariance = target_covariance + np.diag(np.r_[0., source_variances / sample_sizes])
    matrix_variance = float(weights @ estimator_covariance @ weights)
    direct_variance = float(np.var(target_contrast @ weights) / 1000 +
                            np.sum(source_weights ** 2 * source_variances / sample_sizes))
    contrast_error = float(np.max(np.abs(target_contrast @ weights -
                                         (target_treated @ weights - target_control @ weights))))
    assert abs(matrix_variance - direct_variance) < 1e-15 and contrast_error < 1e-15
    return dict(matrix_variance=matrix_variance, direct_variance=direct_variance,
                common_weight_contrast_error=contrast_error,
                scope="fixed nuisance functions and weights; independent sites")


def summary_sufficiency_check():
    """Equal mean/variance summaries need not determine a variance influence."""
    first = np.array([-1., -1., 1., 1.])
    second = np.array([-math.sqrt(2), 0., 0., math.sqrt(2)])
    assert abs(np.mean(first) - np.mean(second)) < 1e-15
    assert abs(np.var(first) - np.var(second)) < 1e-15
    first_variance_if = (first - first.mean()) ** 2 - np.var(first)
    second_variance_if = (second - second.mean()) ** 2 - np.var(second)
    assert abs(np.mean(first_variance_if ** 2)) < 1e-15
    assert abs(np.mean(second_variance_if ** 2) - 1) < 1e-15
    return dict(first=first, second=second, common_mean=0., common_variance=1.,
                first_fourth_moment=float(np.mean(first ** 4)),
                second_fourth_moment=float(np.mean(second ** 4)),
                first_variance_influence_variance=float(np.mean(first_variance_if ** 2)),
                second_variance_influence_variance=float(np.mean(second_variance_if ** 2)))


def inner_outer_variance_check():
    """Initial and final source nuisances induce different TATE objectives.

    Both outcome models are correct, and the target anchor uses its true PS.
    Thus this isolates the variance-objective mismatch from target bias.
    """
    x = np.arange(-2., 3.)
    design = np.column_stack([np.ones(5), x])
    target_mass = np.array([.1, .2, .4, .2, .1])
    source_mass = np.array([.25, .15, .2, .1, .3])
    propensity = np.array([.2, .6, .3, .8, .4])
    source_variance_initial = source_variance_final = target_residual_variance = 0.
    arm_results = []
    for arm, alpha in [(0, np.array([-.9, .3])), (1, np.array([-.4, .7]))]:
        mean = expit(design @ alpha)
        derivative = mean * (1 - mean)
        source_arm_mass = source_mass * (propensity if arm else 1 - propensity)
        initial, _ = tilt_projection(design, target_mass, source_arm_mass, np.ones(5))
        final, _ = tilt_projection(design, target_mass, source_arm_mass, derivative)
        initial_variance = float(np.sum(source_arm_mass * np.exp(-2 * design @ initial) * derivative))
        final_variance = float(np.sum(source_arm_mass * np.exp(-2 * design @ final) * derivative))
        source_variance_initial += initial_variance
        source_variance_final += final_variance
        target_residual_variance += 2 * float(target_mass @ derivative)  # true target PS = 1/2
        arm_results.append(dict(arm=arm, alpha=alpha, gamma_initial=initial, gamma_final=final,
                                source_variance_initial=initial_variance,
                                source_variance_final=final_variance))
    initial_weight = target_residual_variance / (target_residual_variance + source_variance_initial)
    final_weight = target_residual_variance / (target_residual_variance + source_variance_final)
    assert abs(initial_weight - final_weight) > .01
    return dict(arms=arm_results, source_tate_variance_initial=source_variance_initial,
                source_tate_variance_final=source_variance_final,
                target_residual_variance=target_residual_variance,
                initial_objective_weight=initial_weight, final_objective_oracle_weight=final_weight,
                scope="population fixed-nuisance variance objectives, equal site sizes, no screening penalty")


def screening_estimand_check():
    """An inner-compatible source can have a biased final correction.

    The one-round initial OR is the correct target OR. Construct a source
    residual in the null space of the initial weighted OR score. The final OR
    therefore stays unchanged, but the final tilt changes its residual mean.
    """
    x = np.arange(-2., 3.)
    design = np.column_stack([np.ones(5), x])
    target_mass = np.array([.05, .15, .25, .25, .3])
    source_mass = np.array([.3, .25, .2, .15, .1])
    propensity = np.array([.8, .3, .7, .2, .5])
    alpha_initial = np.array([-1.5, 1.2])
    target_mean = expit(design @ alpha_initial)
    gamma_initial, _ = tilt_projection(design, target_mass, source_mass * propensity, np.ones(5))
    gamma_final, _ = tilt_projection(
        design, target_mass, source_mass * propensity, target_mean * (1 - target_mean))
    initial_mass = source_mass * propensity * np.exp(-design @ gamma_initial)
    final_mass = source_mass * propensity * np.exp(-design @ gamma_final)
    constraint = design.T * initial_mass
    residual = final_mass - constraint.T @ np.linalg.solve(constraint @ constraint.T,
                                                           constraint @ final_mass)
    residual /= np.max(np.abs(residual))
    bounds = np.where(residual > 0, (.99 - target_mean) / residual,
                      (.01 - target_mean) / residual)
    residual *= .8 * np.min(bounds)
    source_mean = target_mean + residual
    alpha_final, _ = logistic_projection(design, initial_mass, source_mean)
    final_prediction = expit(design @ alpha_final)
    truth = float(target_mass @ target_mean)
    inner_estimate = truth + float(initial_mass @ residual)
    outer_estimate = float(target_mass @ final_prediction +
                           final_mass @ (source_mean - final_prediction))
    assert np.all((source_mean > .01) & (source_mean < .99))
    assert np.max(np.abs(alpha_final - alpha_initial)) < 1e-9
    assert abs(inner_estimate - truth) < 1e-10 and abs(outer_estimate - truth) > .001
    return dict(x=x, target_mass=target_mass, source_mass=source_mass, propensity=propensity,
                target_mean=target_mean, source_mean=source_mean,
                alpha_initial=alpha_initial, alpha_final=alpha_final,
                gamma_initial=gamma_initial, gamma_final=gamma_final,
                inner_estimate=inner_estimate, outer_estimate=outer_estimate, truth=truth,
                inner_bias=inner_estimate - truth, outer_bias=outer_estimate - truth,
                scope="one-round population update; same TATE bias if control OR is constant and transportable")


def json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="New directory under /scratch.global")
    args = parser.parse_args()
    output = args.output.resolve()
    if Path("/scratch.global") not in output.parents:
        parser.error("Results must be stored under /scratch.global")
    if output.exists():
        parser.error("Output already exists; choose a new run directory")
    checks = {
        "finite_feature_balance": finite_balance_check,
        "initialization_alignment_and_product_remainder": alignment_and_remainder_check,
        "target_anchor_influence": target_anchor_check,
        "clipped_loss_derivative": clipping_check,
        "fixed_cutoff_oracle_weights": fixed_cutoff_check,
        "dimension_uniform_bound": dimension_check,
        "multisite_covariance_identity": multisite_covariance_check,
        "source_summary_sufficiency": summary_sufficiency_check,
        "inner_outer_variance_objective": inner_outer_variance_check,
        "inner_outer_screening_estimand": screening_estimand_check,
    }
    results = {name: check() for name, check in checks.items()}
    output.mkdir(parents=True)
    metadata = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    scope="Independent population fixtures, not production FACE Monte Carlo",
                    results=results)
    (output / "population_checks.json").write_text(
        json.dumps(metadata, indent=2, default=json_value, allow_nan=False) + "\n")
    with (output / "checks.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status"])
        writer.writeheader()
        writer.writerows(dict(check=name, status="expected_identity_or_counterexample_verified")
                         for name in checks)
    print("Verified {} independent population checks; results: {}".format(len(checks), output))


if __name__ == "__main__":
    main()
