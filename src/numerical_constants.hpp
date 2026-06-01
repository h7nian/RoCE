#ifndef NUMERICAL_CONSTANTS_HPP
#define NUMERICAL_CONSTANTS_HPP

// =============================================================================
// NUMERICAL STABILITY CONSTANTS (C++ side)
// =============================================================================
// Centralized numerical constants for FACE-HD algorithms.
// These constants define bounds for numerical stability across the C++ codebase.
//
// IMPORTANT: Keep in sync with R/constants.R which defines the R-side constants.
// Both files should use the same values for cross-language consistency.
// =============================================================================

namespace NumericalConstants {

    // ---- Linear predictor / eta bounds ----
    // Prevents exp() overflow/underflow in logistic and density ratio functions.
    // Used in: density ratio estimation, GLM response functions, CV validation.
    constexpr double ETA_CLIP_MIN = -50.0;
    constexpr double ETA_CLIP_MAX = 50.0;

    // ---- Density ratio weight bounds ----
    // Prevents extreme importance weights from destabilizing estimation.
    constexpr double WEIGHT_MIN = 1e-10;
    constexpr double WEIGHT_MAX = 1e10;

    // ---- Parameter bounds ----
    // Clips coefficient values to prevent overflow in coordinate descent.
    constexpr double PARAM_MAX = 100.0;

    // ---- Minimum variance floor ----
    // Prevents division by zero in variance calculations.
    constexpr double VAR_MIN = 1e-6;

    // ---- Active set threshold ----
    // Coefficients below this threshold are treated as zero in active set CD.
    constexpr double ACTIVE_SET_THRESHOLD = 1e-10;

    // ---- Hessian / general numerical floor ----
    // Minimum value for Hessian diagonals, variance denominators, and other
    // quantities that must stay positive.
    constexpr double HESSIAN_FLOOR = 1e-8;

    // ---- Ratio clipping (correction/variance) ----
    // Clips density ratio * residual terms in correction and variance calculations.
    constexpr double RATIO_CLIP_MIN = -1e6;
    constexpr double RATIO_CLIP_MAX = 1e6;

    // ---- Cross-validation settings ----
    // Convergence tolerance floor for CV (don't need full precision).
    constexpr double CV_TOL_FLOOR = 1e-4;
    // Maximum iterations for CV inner loops.
    constexpr int CV_MAX_ITER = 10000;

    // ---- Active set check frequency ----
    // Check all coordinates (not just active set) every N iterations.
    constexpr int ACTIVE_CHECK_FREQ = 5;

    // ---- Default truncation parameter ----
    // M_tau default for calibrated loss functions (eq:gamma_calibrated_loss in main.tex, truncation function).
    constexpr double M_TAU_DEFAULT = 10.0;

    // ---- Weight optimization bounds ----
    // Maximum absolute aggregation weight (resets to 0 if exceeded).
    constexpr double WEIGHT_MAX_ABS = 10.0;

    // ---- Log-sum-exp approximation threshold ----
    // For eta > this value, log(1+exp(eta)) ≈ eta (avoids exp overflow).
    constexpr double LOG_SUM_EXP_THRESHOLD = 20.0;

    // ---- Step size floor ----
    // Minimum step size in coordinate descent to prevent zero-step stalling.
    constexpr double STEP_SIZE_FLOOR = 1e-12;

    // ---- Probability floor ----
    // Clips predicted probabilities away from 0/1 in binomial log-likelihood.
    constexpr double PROBABILITY_FLOOR = 1e-15;
}

#endif // NUMERICAL_CONSTANTS_HPP
