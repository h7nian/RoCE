#ifndef CV_UTILS_HPP
#define CV_UTILS_HPP

#include "utils.h"
#include "numerical_constants.h"
#include <algorithm>
#include <atomic>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <numeric>  // for std::iota
#include <string>
#include <stdexcept>

#ifdef _OPENMP
#define ROCE_DO_PRAGMA_IMPL(directive) _Pragma(#directive)
#define ROCE_DO_PRAGMA(directive) ROCE_DO_PRAGMA_IMPL(directive)
#define ROCE_PARALLELIZE_CV_FOLDS(thread_count) \
    ROCE_DO_PRAGMA(omp parallel for schedule(static) num_threads(thread_count))
#else
#define ROCE_PARALLELIZE_CV_FOLDS(thread_count) (void)(thread_count);
#endif

// ============================================================================
// Shared utility functions used by density_ratio.cpp and outcome_model.cpp
// ============================================================================
//
// Notation mapping (variable names ↔ paper ↔ GLMUtils method):
//
//   h(η)   = response function = GLMUtils::response_function(η, link)
//   h'(η)  = response derivative = GLMUtils::response_derivative(η, link)
//
//   In the density ratio loss (eq:gamma_loss_final / eq:gamma_calibrated_loss):
//     ψ'(η) ≡ |h'(η)| = b''(θ) > 0  (variables: psi_prime, psi_p, psi_prime_all)
//     Note: std::abs() is applied at precomputation time because the paper's
//     ψ'(θ) = b''(θ) is the second cumulant derivative (always positive), but
//     h'(η) can be negative for non-canonical links (e.g., inverse: h'(η) = -1/η²).
//
//   In the outcome GLM fitting (fit_general_glm_cpp):
//     Uses GLMUtils::irls_working_weight(η, link, family) for Hessian diagonal,
//     which equals h'(η) for canonical links and h'(η)²/V(μ) for non-canonical.
//
//   Both refer to the same mathematical quantity — the derivative of the
//   inverse-link (response) function.  The naming difference traces back to
//   the paper notation (ψ for density ratio loss, h for GLM).
// ============================================================================

// Truncation used by the fold-specific calibration and evaluation scores.
// Following eq:gamma_calibrated_loss in main.tex (truncation function): τ(x) = sign(x)min{|x|, M_τ}
inline double truncation_function(double x, double M_tau = NumericalConstants::M_TAU_DEFAULT) {
    if (x > M_tau) return M_tau;
    if (x < -M_tau) return -M_tau;
    return x;
}

// ============================================================================
// Density ratio weight computation
// ============================================================================
// Shared helper for computing w(X; γ) = exp(−g(Z; γ)) with optional truncation.
// Used by fit_unified_outcome_cpp, select_lambda_cv_*_outcome_cpp, and variance.
//
//   Z_site_int: (n_treated x p_site) intercept-augmented site features
//   gamma_s:    (p_site) density ratio coefficients
//   calibrated: if true, apply \mathcal{T}(·) truncation (eq:alpha_calibrated_loss);
//               if false, apply raw ETA_CLIP (eq:alpha_loss_final)
//   M_tau:      truncation bound (used only when calibrated = true)
//
// Returns (n_treated x 1) weight vector, clipped to [WEIGHT_MIN, WEIGHT_MAX].
// ============================================================================
inline VectorXd compute_density_ratio_weights(
        const MatrixXd& Z_site_int, const VectorXd& gamma_s,
        bool calibrated = false, double M_tau = NumericalConstants::M_TAU_DEFAULT,
        bool use_weight_derivative = false) {
    if (Z_site_int.cols() != gamma_s.size() || !Z_site_int.allFinite() || !gamma_s.allFinite()) {
        throw std::invalid_argument("compute_density_ratio_weights: finite conformable design and coefficients are required.");
    }
    if (use_weight_derivative && (!calibrated || !std::isfinite(M_tau) || M_tau <= 0.0 ||
            M_tau > std::log(NumericalConstants::RATIO_CLIP_MAX))) {
        throw std::invalid_argument("Weight derivatives require a matching finite calibrated truncation radius.");
    }
    int n = Z_site_int.rows();
    VectorXd g_vals = Z_site_int * gamma_s;
    VectorXd original_predictors;
    if (use_weight_derivative) original_predictors = g_vals;

    if (calibrated) {
        // Truncation: τ(g) = sign(g)·min(|g|, M_τ)
        for (int i = 0; i < n; i++) {
            g_vals(i) = truncation_function(g_vals(i), M_tau);
        }
    } else {
        // Raw clipping for numerical stability (no truncation penalty)
        g_vals = g_vals.cwiseMax(NumericalConstants::ETA_CLIP_MIN)
                       .cwiseMin(NumericalConstants::ETA_CLIP_MAX);
    }

    VectorXd weights = (-g_vals.array()).exp();

    // Sanitize: replace NaN/Inf/non-positive with floor
    for (int i = 0; i < n; i++) {
        if (std::isnan(weights(i)) || std::isinf(weights(i)) || weights(i) <= 0) {
            weights(i) = NumericalConstants::VAR_MIN;
        }
    }

    // Clip to [WEIGHT_MIN, WEIGHT_MAX]
    weights = weights.cwiseMax(NumericalConstants::WEIGHT_MIN)
                     .cwiseMin(NumericalConstants::WEIGHT_MAX);
    if (use_weight_derivative) {
        for (int i = 0; i < n; ++i) {
            if (std::abs(original_predictors(i)) >= M_tau) weights(i) = 0.0;
        }
    }
    return weights;
}

// ============================================================================
// Internal CV Utilities (shared infrastructure for all 4 CV functions)
// ============================================================================
// Eliminates ~600 lines of duplicated scaffolding code by factoring out:
//   - Lambda path sorting for warm start
//   - Treated unit filtering and early-return logic
//   - Fold creation and train/val index splitting
//   - Matrix/vector row slicing for fold data
//   - Penalized outcome-GLM fits (coordinate descent and proximal Newton)
//   - Penalized density-ratio fits (coordinate descent and proximal Newton)
//   - Score aggregation and best-lambda selection
// ============================================================================

namespace CVUtils {

// Return the explicitly requested fold-level CV thread count. Source sites
// and treatment arms are already parallelized in R, so the default must stay
// at one to prevent accidental nested oversubscription. Production wrappers
// may opt in through ROCE_NUISANCE_CV_THREADS after requesting enough Slurm
// CPUs. Each thread owns one fold and writes one fold-score column, preserving
// the exact within-fold arithmetic and deterministic aggregate order.
inline int nuisance_cv_thread_count(int n_folds) {
    if (n_folds < 1) {
        throw std::runtime_error(
            "nuisance_cv_thread_count: n_folds must be positive."
        );
    }
    const char* raw_value = std::getenv("ROCE_NUISANCE_CV_THREADS");
    if (raw_value == nullptr || raw_value[0] == '\0') {
        return 1;
    }

    errno = 0;
    char* parse_end = nullptr;
    const long requested = std::strtol(raw_value, &parse_end, 10);
    if (errno != 0 || parse_end == raw_value || *parse_end != '\0' ||
        requested < 1 || requested > INT_MAX) {
        throw std::runtime_error(
            "ROCE_NUISANCE_CV_THREADS must be a positive integer."
        );
    }
#ifndef _OPENMP
    if (requested > 1) {
        throw std::runtime_error(
            "ROCE_NUISANCE_CV_THREADS exceeds one, but RoCE was compiled without OpenMP."
        );
    }
#endif
    return std::min(static_cast<int>(requested), n_folds);
}

// Sort lambda indices from LARGE to SMALL for effective warm starting
inline std::vector<int> sort_lambda_descending(const VectorXd& lambda_grid) {
    int n = lambda_grid.size();
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&lambda_grid](int a, int b) { return lambda_grid(a) > lambda_grid(b); });
    return order;
}

// Filter indices where A == a_val (delegates to filter_treated_indices in utils.h)
inline std::vector<int> filter_treated(const VectorXd& A, int a_val = 1) {
    return filter_treated_indices(A, a_val);
}

// Fold train/val index split
struct FoldSplit {
    std::vector<std::vector<int>> train;
    std::vector<std::vector<int>> val;
};

inline FoldSplit create_fold_splits(
        int n_samples, int n_folds,
        Rcpp::Nullable<Rcpp::NumericVector> cv_fold_id = R_NilValue) {
    // Keep the historical RNG path byte-for-byte unchanged when explicit IDs
    // are absent. Explicit IDs refer to the already arm-filtered row order.
    if (cv_fold_id.isNull()) {
        FoldSplit splits;
        splits.train.resize(n_folds);
        splits.val.resize(n_folds);
        VectorXd fold_ids = create_folds_cpp(n_samples, n_folds);
        for (int fold = 0; fold < n_folds; fold++) {
            for (int j = 0; j < n_samples; j++) {
                if (static_cast<int>(fold_ids(j)) == fold) {
                    splits.val[fold].push_back(j);
                } else {
                    splits.train[fold].push_back(j);
                }
            }
        }
        return splits;
    }

    if (n_samples < 1 || n_folds < 2 || n_samples < n_folds) {
        throw std::runtime_error(
            "create_fold_splits: explicit folds require n_samples >= n_folds >= 2."
        );
    }

    SEXP raw_ids = cv_fold_id.get();
    if ((TYPEOF(raw_ids) != REALSXP && TYPEOF(raw_ids) != INTSXP) ||
        Rf_isFactor(raw_ids) || Rf_getAttrib(raw_ids, R_DimSymbol) != R_NilValue) {
        throw std::runtime_error(
            "cv_fold_id must be a numeric, non-factor vector without dimensions."
        );
    }
    if (Rf_xlength(raw_ids) != static_cast<R_xlen_t>(n_samples)) {
        throw std::runtime_error(
            "cv_fold_id length must equal the number of arm-filtered rows."
        );
    }
    std::vector<int> validated_ids(n_samples);
    std::vector<int> validation_sizes(n_folds, 0);
    for (int j = 0; j < n_samples; j++) {
        const double value = TYPEOF(raw_ids) == REALSXP
            ? REAL(raw_ids)[j]
            : static_cast<double>(INTEGER(raw_ids)[j]);
        if (!std::isfinite(value) || std::floor(value) != value ||
            value < 1.0 || value > static_cast<double>(n_folds)) {
            throw std::runtime_error(
                "cv_fold_id labels must be finite integers in 1..n_folds."
            );
        }
        const int zero_based = static_cast<int>(value) - 1;
        validated_ids[j] = zero_based;
        validation_sizes[zero_based]++;
    }
    for (int fold = 0; fold < n_folds; fold++) {
        if (validation_sizes[fold] == 0 || validation_sizes[fold] == n_samples) {
            throw std::runtime_error(
                "cv_fold_id must give every fold nonempty validation and training sets."
            );
        }
    }
    FoldSplit splits;
    splits.train.resize(n_folds);
    splits.val.resize(n_folds);
    for (int fold = 0; fold < n_folds; fold++) {
        for (int j = 0; j < n_samples; j++) {
            if (validated_ids[j] == fold) {
                splits.val[fold].push_back(j);
            } else {
                splits.train[fold].push_back(j);
            }
        }
    }
    return splits;
}

// Add explicit-fold audit fields without changing the historical NULL schema.
// Fold scores remain equally averaged, so unequal validation-fold sizes do not
// induce observation-weighted CV selection.
inline List append_explicit_fold_audit(
        List result, const FoldSplit& folds,
        Rcpp::Nullable<Rcpp::NumericVector> cv_fold_id) {
    if (cv_fold_id.isNotNull()) {
        result["cv_fold_id"] =
            Rcpp::as<Rcpp::NumericVector>(cv_fold_id.get());
        Rcpp::IntegerVector validation_fold_sizes(folds.val.size());
        for (std::size_t fold = 0; fold < folds.val.size(); fold++) {
            validation_fold_sizes[fold] = folds.val[fold].size();
        }
        result["validation_fold_sizes"] = validation_fold_sizes;
    }
    return result;
}

// Slice rows from a matrix by index vector (delegates to subset_rows in utils.h)
inline MatrixXd slice_rows(const MatrixXd& M, const std::vector<int>& indices) {
    return subset_rows(M, indices);
}

// Slice elements from a vector by index vector (delegates to subset_elements in utils.h)
inline VectorXd slice_elements(const VectorXd& v, const std::vector<int>& indices) {
    return subset_elements(v, indices);
}

// Normalize weights to sum to n (preserve mean weighting)
inline VectorXd normalize_weights(const VectorXd& w) {
    double wsum = w.sum();
    if (wsum > 0) return w / wsum * w.size();
    return w;
}

// Aggregate unpenalized fold scores -> lambda.min/lambda.1se -> selected lambda.
inline List aggregate_cv_results(
        const MatrixXd& fold_scores, const VectorXd& lambda_grid,
        int n_lambda, int n_folds, int path_tail_skipped_fold_fits = 0) {
    if (path_tail_skipped_fold_fits < 0) {
        throw std::runtime_error(
            "aggregate_cv_results: skipped path-tail count must be nonnegative."
        );
    }
    VectorXd cv_scores(n_lambda);
    VectorXd cv_se(n_lambda);
    IntegerVector n_valid(n_lambda);
    for (int i = 0; i < n_lambda; i++) {
        StableAccumulator sum_acc;
        for (int fold = 0; fold < n_folds; fold++) {
            if (std::isfinite(fold_scores(i, fold))) {
                sum_acc.add(fold_scores(i, fold));
                n_valid[i]++;
            }
        }
        if (n_valid[i] != n_folds) {
            // A lambda is eligible only when its optimizer converged and its
            // validation loss was finite in every fold.  Partially valid
            // candidates remain explicitly invalid rather than being averaged
            // over a different fold set.
            cv_scores(i) = INFINITY;
            cv_se(i) = INFINITY;
            continue;
        }
        cv_scores(i) = sum_acc.value() / n_folds;

        if (n_folds > 1) {
            StableAccumulator var_acc;
            for (int fold = 0; fold < n_folds; fold++) {
                double diff = fold_scores(i, fold) - cv_scores(i);
                var_acc.add(diff * diff);
            }
            double fold_var = var_acc.value() / (n_folds - 1);
            cv_se(i) = std::sqrt(std::max(0.0, fold_var) / n_folds);
        } else {
            cv_se(i) = 0.0;
        }
    }

    int best_idx = -1;
    double min_score = INFINITY;
    for (int i = 0; i < n_lambda; i++) {
        if (std::isfinite(cv_scores(i)) && cv_scores(i) < min_score) {
            min_score = cv_scores(i);
            best_idx = i;
        }
    }
    if (best_idx < 0) {
        throw std::runtime_error(
            "aggregate_cv_results: no lambda converged with a finite validation loss across every CV fold."
        );
    }

    int idx_1se = best_idx;
    if (std::isfinite(min_score)) {
        double one_se_cutoff = min_score + cv_se(best_idx);
        double largest_lambda = lambda_grid(best_idx);
        for (int i = 0; i < n_lambda; i++) {
            if (std::isfinite(cv_scores(i)) &&
                cv_scores(i) <= one_se_cutoff &&
                lambda_grid(i) > largest_lambda) {
                largest_lambda = lambda_grid(i);
                idx_1se = i;
            }
        }
    }

    int valid_fold_fit_count = 0;
    int invalid_lambda_count = 0;
    for (int i = 0; i < n_lambda; i++) {
        valid_fold_fit_count += n_valid[i];
        if (n_valid[i] < n_folds) invalid_lambda_count++;
    }

    return List::create(Named("best_lambda") = lambda_grid(best_idx),
                       Named("best_idx") = best_idx + 1,
                       Named("lambda_min") = lambda_grid(best_idx),
                       Named("lambda_1se") = lambda_grid(idx_1se),
                       Named("idx_min") = best_idx + 1,
                       Named("idx_1se") = idx_1se + 1,
                       Named("cv_scores") = cv_scores,
                       Named("cv_se") = cv_se,
                       Named("n_valid_folds") = n_valid,
                       Named("invalid_fold_fits") =
                           n_lambda * n_folds - valid_fold_fit_count,
                       Named("invalid_lambdas") =
                           invalid_lambda_count,
                       Named("path_tail_skipped_fold_fits") =
                           path_tail_skipped_fold_fits,
                       Named("lambda_rule") = "min");
}

// Record one failed candidate along a descending lambda path. Returning a
// positive value means the patience threshold was reached; that many smaller
// candidate/fold fits can be left at their preallocated +Inf score and skipped.
inline int cv_path_tail_after_failure(
        int& consecutive_failures, int path_position, int n_lambda) {
    consecutive_failures++;
    if (consecutive_failures < NumericalConstants::CV_FAILURE_PATIENCE) {
        return 0;
    }
    return std::max(0, n_lambda - path_position - 1);
}

// ---- Penalized nuisance fits: coordinate descent and proximal Newton ----

// A small parameter step alone does not establish stationarity, particularly
// at a numerical coefficient bound. The intercept remains unpenalized.
inline double l1_score_residual(const VectorXd& coefficient,
                                const VectorXd& gradient, double lambda) {
    if (!coefficient.allFinite() || !gradient.allFinite()) return INFINITY;
    double residual = std::abs(gradient(0));
    for (int j = 1; j < coefficient.size(); ++j) {
        const double violation = coefficient(j) == 0.0 ?
            std::max(std::abs(gradient(j)) - lambda, 0.0) :
            std::abs(gradient(j) + std::copysign(lambda, coefficient(j)));
        residual = std::max(residual, violation);
    }
    return residual;
}

inline double composite_linear_change(const VectorXd& coefficient,
        const VectorXd& delta, const VectorXd& gradient, double lambda) {
    StableAccumulator change;
    change.add(gradient(0) * delta(0));
    for (int j = 1; j < coefficient.size(); ++j) {
        const double next = coefficient(j) + delta(j);
        if (coefficient(j) * next >= 0.0) {
            // Combine cancelling smooth and L1 slopes before multiplying by
            // a tiny Newton step. Subtracting two full L1 norms can erase the
            // model decrease near a high-precision stationary solution.
            const double sign_point = coefficient(j) == 0.0 ? next : coefficient(j);
            const double penalty_slope = sign_point == 0.0 ? 0.0 : std::copysign(lambda, sign_point);
            change.add((gradient(j) + penalty_slope) * delta(j));
        } else {
            change.add(gradient(j) * delta(j) +
                lambda * (std::abs(next) - std::abs(coefficient(j))));
        }
    }
    return change.value();
}

inline double glm_score_residual(const VectorXd& beta, const MatrixXd& X,
        const VectorXd& y, const VectorXd& weights, double lambda, LinkFunction link) {
    const VectorXd predictor = X * beta;
    VectorXd residual(y.size());
    for (int i = 0; i < y.size(); ++i) {
        residual(i) = weights(i) * GLMUtils::nll_gradient_residual(
            y(i), GLMUtils::response_function(predictor(i), link), link);
    }
    const VectorXd gradient = X.transpose() * residual / X.rows();
    return l1_score_residual(beta, gradient, lambda);
}

struct GLMFitResult {
    bool converged = false;
    int iterations = 0;
    double max_update = NA_REAL;
    double convergence_threshold = NA_REAL;
    int line_search_failures = 0;
    double kkt_residual = NA_REAL;
    double kkt_threshold = NA_REAL;
};

// Solver for every L1-penalized nuisance fit (outcome GLM and density ratio;
// CV folds and final refits).  The proximal-Newton path is the default;
// ROCE_NUISANCE_SOLVER=coordinate_descent selects the original coordinate
// descent, which minimizes the same objectives and is kept for paired
// equivalence audits against earlier production runs.
// R fitting calls install a scoped override before starting CV workers. All
// OpenMP workers read the same selection; process workers inherit it on fork.
// An atomic protects the shared setting without changing solver arithmetic.
inline std::atomic<int>& nuisance_solver_override() {
    static std::atomic<int> selection{-1};  // -1: environment, 0: CD, 1: Newton
    return selection;
}

inline bool nuisance_use_proximal_newton() {
    const int selection = nuisance_solver_override().load(std::memory_order_relaxed);
    if (selection >= 0) return selection == 1;
    static const bool use_newton = [] {
        const char* value = std::getenv("ROCE_NUISANCE_SOLVER");
        return value == nullptr || std::string(value) != "coordinate_descent";
    }();
    return use_newton;
}

inline const char* nuisance_solver_name() {
    return nuisance_use_proximal_newton() ? "proximal_newton" : "coordinate_descent";
}

// Shared proximal-Newton tuning.  The Levenberg-Marquardt damping is relative
// to the largest Hessian diagonal: a weighted Gram matrix is near-singular
// when p approaches the fitting sample size, so an undamped Newton direction
// can be arbitrarily long in flat directions and defeat the backtracking
// search.  Backtracking keeps iterates strictly inside the coefficient box, so
// a solution pinned at the bound sits just below PARAM_MAX; such a solution is
// reported as non-converged because coordinate descent never certifies
// convergence there either, which keeps the CV eligibility rule solver-invariant.
namespace ProximalNewton {
    constexpr int inner_max_sweeps = 1000;
    constexpr int max_backtracking_steps = 30;
    constexpr double backtracking_factor = 0.5;
    constexpr double armijo_constant = 1e-4;
    constexpr double damping_initial = 1e-3;
    constexpr double damping_minimum = 1e-6;
    constexpr double damping_maximum = 1e6;
    constexpr double damping_growth = 10.0;
    constexpr double damping_decay = 0.3;
    constexpr double bound_fraction = 0.99;
    // Inexact Newton: the quadratic model is solved only to a fraction of the
    // previous step length, tightening to the final tolerance as the outer
    // iteration converges.
    constexpr double inner_tol_initial = 1e-2;
    constexpr double inner_tol_step_fraction = 0.05;

    // Minimize grad'delta + delta'(H + ridge I) delta / 2 + lambda * ||theta + delta||_1
    // over delta by coordinate descent, where H = X' diag(curvature) X.
    // curvature_X = diag(curvature) X and hessian_diag = diag(H) are
    // precomputed by the caller; z = X delta is kept incrementally.
    // Sweeps alternate glmnet-style: one complete sweep fixes the active set
    // (nonzero coordinates plus the intercept), the active set is iterated to
    // tolerance, and a further complete sweep certifies the remaining
    // coordinates' optimality conditions.
    inline void inner_coordinate_descent(
            VectorXd& delta, VectorXd& z_vec, const VectorXd& theta,
            const VectorXd& gradient, const VectorXd& hessian_diag,
            const MatrixXd& X, const MatrixXd& curvature_X,
            double ridge, double lambda, double inner_tol) {
        const int p = theta.size();
        delta.setZero();
        z_vec.setZero();
        auto update_coordinate = [&](int j) {
            const double curvature_j = std::max(
                hessian_diag(j) + ridge, NumericalConstants::HESSIAN_FLOOR
            );
            const double partial = gradient(j) +
                curvature_X.col(j).dot(z_vec) + ridge * delta(j);
            const double current = theta(j) + delta(j);
            const double proposal = current - partial / curvature_j;
            const double updated = (j == 0) ? proposal :
                soft_threshold_cpp(proposal, lambda / curvature_j);
            const double change = updated - current;
            if (change != 0.0) {
                delta(j) += change;
                z_vec += change * X.col(j);
            }
            return std::abs(change);
        };

        std::vector<int> active_set;
        active_set.reserve(p);
        int sweeps = 0;
        while (sweeps < inner_max_sweeps) {
            double max_update = 0.0;
            sweeps++;
            active_set.clear();
            for (int j = 0; j < p; j++) {
                max_update = std::max(max_update, update_coordinate(j));
                if (j == 0 || theta(j) + delta(j) != 0.0) active_set.push_back(j);
            }
            if (max_update <= inner_tol) break;

            bool active_converged = false;
            while (sweeps < inner_max_sweeps) {
                double max_active_update = 0.0;
                sweeps++;
                for (int j : active_set) {
                    max_active_update = std::max(max_active_update, update_coordinate(j));
                }
                if (max_active_update <= inner_tol) {
                    active_converged = true;
                    break;
                }
            }
            if (!active_converged) break;
        }
    }

    inline double l1_penalty(const VectorXd& theta, double lambda) {
        double total = 0.0;
        for (int j = 1; j < theta.size(); j++) total += std::abs(theta(j));
        return lambda * total;
    }
}

// L1-penalized GLM coordinate descent (shared by refined & calibrated outcome
// CV). Returning convergence metadata prevents a finite CV path from silently
// accepting a lambda whose coordinate-descent refit exhausted its iteration
// budget.
inline GLMFitResult glm_coordinate_descent(
        VectorXd& beta, std::vector<bool>& active,
        const MatrixXd& X_train, const VectorXd& Y_train,
        const VectorXd& weights_train, int n_train,
        double lambda, LinkFunction link, GLMFamily family,
        double cv_tol, int cv_max_iter) {
    int p = beta.size();
    GLMFitResult result;

    for (int iter = 0; iter < cv_max_iter; iter++) {
        VectorXd beta_old = beta;
        const bool full_scan =
            iter % NumericalConstants::ACTIVE_CHECK_FREQ == 0 ||
            iter == cv_max_iter - 1;

        if (full_scan) {
            for (int j = 1; j < p; j++) {
                active[j] = (std::abs(beta(j)) > NumericalConstants::ACTIVE_SET_THRESHOLD) || (iter == 0);
            }
        }

        // Compute eta = X_train * beta once per outer iteration (O(np)),
        // then update incrementally after each coordinate change (O(n) per coordinate).
        // This reduces total complexity from O(np²) to O(np) per iteration.
        VectorXd eta = X_train * beta;

        for (int j = 0; j < p; j++) {
            if (j > 0 && !active[j] && !full_scan) continue;

            StableAccumulator grad_acc;
            StableAccumulator hess_acc;

            for (int i = 0; i < n_train; i++) {
                double x_ij = X_train(i, j);
                double w_i = weights_train(i);
                double mu_i = GLMUtils::response_function(eta(i), link);
                double ww_i = GLMUtils::irls_working_weight(eta(i), link, family);
                double r_i = GLMUtils::nll_gradient_residual(Y_train(i), mu_i, link);
                grad_acc.add(w_i * r_i * x_ij);
                hess_acc.add(w_i * ww_i * x_ij * x_ij);
            }

            double grad_j = grad_acc.value() / n_train;
            double hess_j = hess_acc.value() / n_train;
            hess_j = std::max(hess_j, NumericalConstants::HESSIAN_FLOOR);

            double step_size = 1.0 / hess_j;

            double old_beta_j = beta(j);
            beta(j) = sanitize_param(proximal_update(beta(j), grad_j, step_size, lambda, j));

            // Incremental eta update: O(n) instead of recomputing O(np)
            double delta_j = beta(j) - old_beta_j;
            if (delta_j != 0.0) {
                eta += delta_j * X_train.col(j);
                // A coordinate that was zero at the beginning of this full
                // scan may become non-zero after its KKT update.  Keep it in
                // the active set immediately; otherwise it would be skipped
                // until the next periodic full scan and the path solver could
                // declare convergence with a newly activated coefficient only
                // partially updated.
                if (j > 0) {
                    active[j] =
                        std::abs(beta(j)) > NumericalConstants::ACTIVE_SET_THRESHOLD;
                }
            }
        }

        // Convergence is certified only after a full KKT scan.  Between full
        // scans, an inactive coordinate is deliberately skipped and therefore
        // cannot be used to certify convergence of the complete problem.
        result.iterations = iter + 1;
        result.max_update = (beta - beta_old).lpNorm<Eigen::Infinity>();
        result.convergence_threshold = convergence_threshold_cpp(beta, cv_tol);
        if (beta.lpNorm<Eigen::Infinity>() >=
                NumericalConstants::PARAM_MAX * ProximalNewton::bound_fraction) break;
        if (full_scan && check_convergence_cpp(beta_old, beta, cv_tol) &&
                glm_score_residual(beta, X_train, Y_train, weights_train, lambda, link) <= cv_tol) {
            result.converged = true;
            break;
        }
    }
    return result;
}

// Proximal-Newton solver for the same L1-penalized weighted GLM objective
// (1/n) sum_i w_i loss(y_i, h(eta_i)) + lambda ||beta_{-0}||_1.  Each outer
// step builds the IRLS quadratic model at the current beta, minimizes model +
// L1 penalty by coordinate descent (no link evaluations inside the inner loop),
// and accepts the step by an Armijo backtracking search on the complete
// penalized objective.  The gradient (mu - y) x and working weight h'(eta) are
// exact derivatives of the loss for the canonical links used in production.
inline GLMFitResult glm_proximal_newton(
        VectorXd& beta, std::vector<bool>& active,
        const MatrixXd& X_train, const VectorXd& Y_train,
        const VectorXd& weights_train, int n_train,
        double lambda, LinkFunction link, GLMFamily family,
        double cv_tol, int cv_max_iter) {
    using namespace ProximalNewton;
    const int p = beta.size();
    GLMFitResult result;
    const double inner_tol_final = std::max(0.1 * cv_tol, 1e-12);
    double inner_tol = std::max(inner_tol_final, inner_tol_initial);

    auto loss_value = [&](const VectorXd& eta) {
        StableAccumulator acc;
        for (int i = 0; i < n_train; i++) {
            const double mu_i = GLMUtils::response_function(eta(i), link);
            acc.add(weights_train(i) * GLMUtils::calculate_loss(Y_train(i), mu_i, family) / n_train);
        }
        return acc.value();
    };

    VectorXd eta = X_train * beta;
    double objective = loss_value(eta) + l1_penalty(beta, lambda);

    VectorXd residual_vec(n_train), curvature_vec(n_train);
    VectorXd gradient(p), hessian_diag(p), delta(p), z_vec(n_train);
    MatrixXd curvature_X(n_train, p);
    double damping = damping_initial;

    for (int iter = 0; iter < cv_max_iter; iter++) {
        VectorXd beta_old = beta;

        for (int i = 0; i < n_train; i++) {
            const double mu_i = GLMUtils::response_function(eta(i), link);
            residual_vec(i) = weights_train(i) *
                GLMUtils::nll_gradient_residual(Y_train(i), mu_i, link) / n_train;
            curvature_vec(i) = weights_train(i) *
                GLMUtils::irls_working_weight(eta(i), link, family) / n_train;
        }
        for (int j = 0; j < p; j++) {
            StableAccumulator acc;
            for (int i = 0; i < n_train; i++) acc.add(X_train(i, j) * residual_vec(i));
            gradient(j) = acc.value();
        }
        curvature_X = curvature_vec.asDiagonal() * X_train;
        double hessian_scale = NumericalConstants::HESSIAN_FLOOR;
        for (int j = 0; j < p; j++) {
            hessian_diag(j) = curvature_X.col(j).dot(X_train.col(j));
            hessian_scale = std::max(hessian_scale, hessian_diag(j));
        }

        bool accepted = false;
        bool stationary = false;
        while (!accepted) {
            inner_coordinate_descent(
                delta, z_vec, beta, gradient, hessian_diag, X_train, curvature_X,
                damping * hessian_scale, lambda, inner_tol
            );
            const double model_decrease = composite_linear_change(beta, delta, gradient, lambda);
            if (delta.lpNorm<Eigen::Infinity>() == 0.0 || model_decrease >= 0.0) {
                stationary = true;
                break;
            }

            double step = 1.0;
            for (int backtrack = 0; backtrack < max_backtracking_steps; backtrack++) {
                VectorXd beta_candidate = beta + step * delta;
                if (beta_candidate.lpNorm<Eigen::Infinity>() > NumericalConstants::PARAM_MAX) {
                    step *= backtracking_factor;
                    continue;
                }
                VectorXd eta_candidate = eta + step * z_vec;
                const double objective_candidate =
                    loss_value(eta_candidate) + l1_penalty(beta_candidate, lambda);
                const double sufficient = objective + armijo_constant * step * model_decrease +
                    1e-12 * (1.0 + std::abs(objective));
                if (std::isfinite(objective_candidate) && objective_candidate <= sufficient) {
                    beta = beta_candidate;
                    eta = eta_candidate;
                    objective = objective_candidate;
                    accepted = true;
                    break;
                }
                step *= backtracking_factor;
            }

            if (!accepted) {
                if (damping >= damping_maximum) break;
                damping = std::min(damping * damping_growth, damping_maximum);
            } else if (step == 1.0) {
                damping = std::max(damping * damping_decay, damping_minimum);
            }
        }

        result.iterations = iter + 1;
        result.max_update = (beta - beta_old).lpNorm<Eigen::Infinity>();
        result.convergence_threshold = convergence_threshold_cpp(beta, cv_tol);
        if (stationary) {
            result.converged = glm_score_residual(
                beta, X_train, Y_train, weights_train, lambda, link) <= cv_tol;
            if (result.converged || inner_tol <= inner_tol_final) break;
            inner_tol = inner_tol_final;
            continue;
        }
        if (!accepted) {
            result.line_search_failures++;
            break;
        }
        // Reaching the coefficient bound is final (see below); stop iterating.
        if (beta.lpNorm<Eigen::Infinity>() >= NumericalConstants::PARAM_MAX * bound_fraction) {
            break;
        }
        if (check_convergence_cpp(beta_old, beta, cv_tol) &&
                glm_score_residual(beta, X_train, Y_train, weights_train, lambda, link) <= cv_tol) {
            result.converged = true;
            break;
        }
        inner_tol = std::min(
            inner_tol_initial,
            std::max(inner_tol_final, inner_tol_step_fraction * result.max_update)
        );
    }

    if (beta.lpNorm<Eigen::Infinity>() >= NumericalConstants::PARAM_MAX * bound_fraction) {
        result.converged = false;
    }
    for (int j = 1; j < p; j++) {
        active[j] = std::abs(beta(j)) > NumericalConstants::ACTIVE_SET_THRESHOLD;
    }
    return result;
}

inline GLMFitResult glm_penalized_fit(
        VectorXd& beta, std::vector<bool>& active,
        const MatrixXd& X_train, const VectorXd& Y_train,
        const VectorXd& weights_train, int n_train,
        double lambda, LinkFunction link, GLMFamily family,
        double cv_tol, int cv_max_iter) {
    GLMFitResult result;
    if (nuisance_use_proximal_newton()) {
        result = glm_proximal_newton(
            beta, active, X_train, Y_train, weights_train, n_train,
            lambda, link, family, cv_tol, cv_max_iter
        );
    } else {
        result = glm_coordinate_descent(
            beta, active, X_train, Y_train, weights_train, n_train,
            lambda, link, family, cv_tol, cv_max_iter);
    }
    result.kkt_residual = glm_score_residual(beta, X_train, Y_train, weights_train, lambda, link);
    result.kkt_threshold = cv_tol;
    result.converged = result.converged && result.kkt_residual <= result.kkt_threshold;
    return result;
}

struct DensityRatioFitResult {
    bool converged = false;
    int iterations = 0;
    double max_update = NA_REAL;
    double convergence_threshold = NA_REAL;
    int line_search_failures = 0;
    double kkt_residual = NA_REAL;
    double kkt_threshold = NA_REAL;
};

// A sufficient penalty bound below which no point can satisfy the existing
// absolute KKT tolerance. X contains an unpenalized intercept and the source
// score has form moment - X' w with nonnegative w. A previous iterate supplies
// a direction; this does not alter either the loss or the optimization start.
inline double density_cv_kkt_bound(const VectorXd& gamma, const MatrixXd& X,
                                  const VectorXd& moment, double tolerance) {
    if (X.rows() == 0 || X.cols() < 2 || gamma.size() != X.cols() ||
            moment.size() != X.cols() || !gamma.allFinite() || !X.allFinite() ||
            !moment.allFinite() || moment(0) <= 0 ||
            !std::isfinite(tolerance) || tolerance <= 0 ||
            (X.col(0).array() != 1.0).any()) return 0.0;
    const int slopes_count = X.cols() - 1;
    const double coefficient_scale = gamma.tail(slopes_count).cwiseAbs().maxCoeff();
    if (coefficient_scale == 0) return 0.0;
    VectorXd direction = gamma.tail(slopes_count) / coefficient_scale;
    direction /= direction.cwiseAbs().sum();
    const double slope_norm = direction.cwiseAbs().sum();
    const double scale = std::max(1.0, std::max(X.cwiseAbs().maxCoeff(),
                                             moment.cwiseAbs().maxCoeff()));
    const double padding = 64.0 * std::numeric_limits<double>::epsilon() *
        (static_cast<double>(X.rows()) + X.cols() + 1.0) * scale;
    const VectorXd projection = X.rightCols(slopes_count) * direction;
    if (!projection.allFinite() || !std::isfinite(padding)) return 0.0;
    const double intercept = -projection.minCoeff() + padding;
    const double roundoff = padding * (1.0 + moment.cwiseAbs().sum());
    const double bound = (-moment(0) * intercept -
        moment.tail(slopes_count).dot(direction) -
        tolerance * (std::abs(intercept) + slope_norm) - roundoff) / slope_norm;
    return std::isfinite(bound) ? std::max(0.0, bound) : 0.0;
}

// Tilting weight exp{-T_M(phi'gamma)} shared by every density-ratio loss, its
// score, and the influence function (main.tex, truncation paragraph in
// sec:nuisance). ETA_CLIP only guards exp() overflow when M_tau is infinite.
inline double tilt_weight(double linear_predictor, double M_tau) {
    double clipped = std::max(
        NumericalConstants::ETA_CLIP_MIN,
        std::min(NumericalConstants::ETA_CLIP_MAX,
                 truncation_function(linear_predictor, M_tau))
    );
    double weight = std::exp(-clipped);
    return std::max(
        NumericalConstants::WEIGHT_MIN,
        std::min(NumericalConstants::WEIGHT_MAX, weight)
    );
}

// Loss term whose derivative is -tilt_weight: exp(-g) inside the truncation
// radius and the tangent continuation outside it, so that coordinate descent
// minimizes exactly the objective whose score uses the truncated weights.
inline double tilt_loss(double linear_predictor, double M_tau) {
    double truncated = truncation_function(linear_predictor, M_tau);
    return tilt_weight(truncated, M_tau) * (1.0 - (linear_predictor - truncated));
}

inline double density_ratio_score_residual(const VectorXd& gamma,
        const MatrixXd& X, const VectorXd& psi_prime, double source_scale,
        const VectorXd& moment, double lambda, double radius) {
    const VectorXd predictor = X * gamma;
    VectorXd weight(X.rows());
    for (int i = 0; i < X.rows(); ++i) {
        weight(i) = source_scale * psi_prime(i) * tilt_weight(predictor(i), radius);
    }
    const VectorXd gradient = moment - X.transpose() * weight / X.rows();
    return l1_score_residual(gamma, gradient, lambda);
}

// L1-penalized density-ratio coordinate descent shared by CV and final refits.
// source_scale is the empirical arm fraction n_{s,a}/n_s.  Consequently,
// source_scale * mean_{i:A_i=a}(.) is the source-site empirical expectation
// n_s^{-1} sum_i I(A_i=a)(.) both for a complete sample and for an inner-CV
// training fold.  A monotone backtracking search replaces the former ad-hoc
// iteration-decaying cap, which could stall far above the requested tolerance.
inline DensityRatioFitResult density_ratio_coordinate_descent(
        VectorXd& gamma, std::vector<bool>& active,
        const MatrixXd& X_train, const VectorXd& psi_prime_train,
        double source_scale, const VectorXd& mean_grad_psi, double lambda,
        double M_tau, double convergence_tol, int max_iter,
        bool use_active_set = true) {
    int p = gamma.size();
    int n_train = X_train.rows();

    if (n_train <= 0 || psi_prime_train.size() != n_train) {
        throw std::runtime_error(
            "density_ratio_coordinate_descent: training rows and psi_prime must be non-empty and match."
        );
    }
    if (mean_grad_psi.size() != p) {
        throw std::runtime_error(
            "density_ratio_coordinate_descent: mean_grad_psi length must match gamma length."
        );
    }
    if (!std::isfinite(source_scale) || source_scale <= 0.0 || source_scale > 1.0) {
        throw std::runtime_error(
            "density_ratio_coordinate_descent: source_scale must be in (0, 1]."
        );
    }

    DensityRatioFitResult result;
    constexpr int max_backtracking_steps = 50;
    constexpr double backtracking_factor = 0.5;

    for (int iter = 0; iter < max_iter; iter++) {
        VectorXd gamma_old = gamma;
        bool sweep_line_search_failed = false;
        const bool full_scan = !use_active_set ||
            iter % NumericalConstants::ACTIVE_CHECK_FREQ == 0 ||
            iter == max_iter - 1;

        if (use_active_set && full_scan) {
            for (int j = 1; j < p; j++) {
                active[j] = (std::abs(gamma(j)) > NumericalConstants::ACTIVE_SET_THRESHOLD) || (iter == 0);
            }
        }

        // Compute g = X_train * gamma once per outer iteration (O(np)),
        // then update incrementally after each coordinate change (O(n) per coordinate).
        // This reduces total complexity from O(np²) to O(np) per iteration.
        VectorXd g_vec = X_train * gamma;
        VectorXd weight_vec(n_train);
        for (int i = 0; i < n_train; i++) {
            weight_vec(i) = tilt_weight(g_vec(i), M_tau);
        }

        for (int j = 0; j < p; j++) {
            if (use_active_set && j > 0 && !active[j] && !full_scan) continue;

            // Compute gradient and Hessian-diagonal approximation for coordinate j
            // grad_j = mean_grad_psi(j) - source_scale *
            //          mean_{A=a}[X_j * exp(-φ^T γ) * ψ']
            // hess_j is the corresponding diagonal curvature.
            double grad_j = mean_grad_psi(j);
            StableAccumulator grad_acc;
            StableAccumulator hess_acc;

            for (int i = 0; i < n_train; i++) {
                double exp_neg_g = weight_vec(i);
                double psi_p = psi_prime_train(i);

                // first derivative (gradient) contribution
                grad_acc.add(
                    -source_scale * X_train(i, j) * exp_neg_g * psi_p / n_train
                );
                // second derivative (Hessian diagonal) contribution; the
                // truncated loss is linear, hence flat, beyond the radius.
                double curvature =
                    (std::abs(g_vec(i)) < M_tau) ? exp_neg_g : 0.0;
                hess_acc.add(
                    source_scale * X_train(i, j) * X_train(i, j) *
                    curvature * psi_p / n_train
                );
            }
            grad_j += grad_acc.value();
            double step_size = 1.0 / std::max(
                hess_acc.value(), NumericalConstants::HESSIAN_FLOOR
            );

            // Propose a proximal Newton coordinate update, then backtrack until
            // the complete one-coordinate penalized objective does not increase.
            double old_gamma_j = gamma(j);
            bool accepted = false;
            VectorXd candidate_weights(n_train);
            for (int backtrack = 0; backtrack < max_backtracking_steps; backtrack++) {
                double candidate = proximal_update(
                    old_gamma_j, grad_j, step_size, lambda, j
                );
                if (!std::isfinite(candidate) ||
                        std::abs(candidate) > NumericalConstants::PARAM_MAX) {
                    step_size *= backtracking_factor;
                    continue;
                }

                double delta_j = candidate - old_gamma_j;
                if (delta_j == 0.0) {
                    accepted = true;
                    break;
                }

                StableAccumulator source_loss_change;
                for (int i = 0; i < n_train; i++) {
                    double new_predictor = g_vec(i) + delta_j * X_train(i, j);
                    candidate_weights(i) = tilt_weight(new_predictor, M_tau);
                    source_loss_change.add(
                        source_scale * psi_prime_train(i) *
                        (tilt_loss(new_predictor, M_tau) -
                         tilt_loss(g_vec(i), M_tau)) / n_train
                    );
                }
                double penalty_change = (j == 0) ? 0.0 :
                    lambda * (std::abs(candidate) - std::abs(old_gamma_j));
                double objective_change = mean_grad_psi(j) * delta_j +
                    source_loss_change.value() + penalty_change;
                double numerical_tolerance = 1e-12 *
                    (1.0 + std::abs(old_gamma_j) + std::abs(candidate));

                if (std::isfinite(objective_change) &&
                        objective_change <= numerical_tolerance) {
                    gamma(j) = candidate;
                    g_vec += delta_j * X_train.col(j);
                    weight_vec.swap(candidate_weights);
                    if (use_active_set && j > 0) {
                        active[j] = std::abs(gamma(j)) >
                            NumericalConstants::ACTIVE_SET_THRESHOLD;
                    }
                    accepted = true;
                    break;
                }
                step_size *= backtracking_factor;
                if (step_size < NumericalConstants::STEP_SIZE_FLOOR) break;
            }

            if (!accepted) {
                sweep_line_search_failed = true;
                result.line_search_failures++;
            }
        }

        result.iterations = iter + 1;
        result.max_update = (gamma - gamma_old).lpNorm<Eigen::Infinity>();
        result.convergence_threshold = convergence_threshold_cpp(
            gamma, convergence_tol
        );
        if (gamma.lpNorm<Eigen::Infinity>() >=
                NumericalConstants::PARAM_MAX * ProximalNewton::bound_fraction) break;
        if (!sweep_line_search_failed && full_scan &&
                check_convergence_cpp(gamma_old, gamma, convergence_tol) &&
                density_ratio_score_residual(gamma, X_train, psi_prime_train,
                    source_scale, mean_grad_psi, lambda, M_tau) <= convergence_tol) {
            result.converged = true;
            break;
        }
    }

    return result;
}

// Proximal-Newton solver for the same L1-penalized density-ratio objective.
// Each outer step builds the exact diagonal-weighted quadratic model of the
// smooth part at the current gamma, minimizes model + L1 penalty by coordinate
// descent (no exponentials inside the inner loop), and accepts the step by an
// Armijo backtracking search on the complete penalized objective.  Tilt
// weights and losses are the same functions as in the coordinate-descent path,
// so both solvers minimize the identical objective.
inline DensityRatioFitResult density_ratio_proximal_newton(
        VectorXd& gamma, std::vector<bool>& active,
        const MatrixXd& X_train, const VectorXd& psi_prime_train,
        double source_scale, const VectorXd& mean_grad_psi, double lambda,
        double M_tau, double convergence_tol, int max_iter,
        bool use_active_set = true) {
    const int p = gamma.size();
    const int n_train = X_train.rows();

    if (n_train <= 0 || psi_prime_train.size() != n_train) {
        throw std::runtime_error(
            "density_ratio_proximal_newton: training rows and psi_prime must be non-empty and match."
        );
    }
    if (mean_grad_psi.size() != p) {
        throw std::runtime_error(
            "density_ratio_proximal_newton: mean_grad_psi length must match gamma length."
        );
    }
    if (!std::isfinite(source_scale) || source_scale <= 0.0 || source_scale > 1.0) {
        throw std::runtime_error(
            "density_ratio_proximal_newton: source_scale must be in (0, 1]."
        );
    }

    using namespace ProximalNewton;
    DensityRatioFitResult result;
    const double inner_tol_final = std::max(0.1 * convergence_tol, 1e-12);
    double inner_tol = std::max(inner_tol_final, inner_tol_initial);

    auto source_loss_value = [&](const VectorXd& g) {
        StableAccumulator acc;
        for (int i = 0; i < n_train; i++) {
            acc.add(source_scale * psi_prime_train(i) * tilt_loss(g(i), M_tau) / n_train);
        }
        return acc.value();
    };

    VectorXd g_vec = X_train * gamma;
    double objective = mean_grad_psi.dot(gamma) + source_loss_value(g_vec) +
        l1_penalty(gamma, lambda);

    VectorXd weight_vec(n_train), curvature_vec(n_train);
    VectorXd gradient(p), hessian_diag(p), delta(p), z_vec(n_train);
    MatrixXd curvature_X(n_train, p);
    double damping = damping_initial;

    for (int iter = 0; iter < max_iter; iter++) {
        VectorXd gamma_old = gamma;

        for (int i = 0; i < n_train; i++) {
            const double tilt = source_scale * psi_prime_train(i) *
                tilt_weight(g_vec(i), M_tau) / n_train;
            weight_vec(i) = tilt;
            // The truncated loss is linear, hence flat, beyond the radius.
            curvature_vec(i) = (std::abs(g_vec(i)) < M_tau) ? tilt : 0.0;
        }
        // Compensated column sums keep an exactly balanced moment at an exact
        // zero gradient, as the coordinate-descent path does.
        for (int j = 0; j < p; j++) {
            StableAccumulator acc;
            for (int i = 0; i < n_train; i++) acc.add(X_train(i, j) * weight_vec(i));
            gradient(j) = mean_grad_psi(j) - acc.value();
        }
        curvature_X = curvature_vec.asDiagonal() * X_train;
        double hessian_scale = NumericalConstants::HESSIAN_FLOOR;
        for (int j = 0; j < p; j++) {
            hessian_diag(j) = curvature_X.col(j).dot(X_train.col(j));
            hessian_scale = std::max(hessian_scale, hessian_diag(j));
        }

        bool accepted = false;
        bool stationary = false;
        while (!accepted) {
            inner_coordinate_descent(
                delta, z_vec, gamma, gradient, hessian_diag, X_train, curvature_X,
                damping * hessian_scale, lambda, inner_tol
            );
            const double model_decrease = composite_linear_change(gamma, delta, gradient, lambda);
            if (delta.lpNorm<Eigen::Infinity>() == 0.0 || model_decrease >= 0.0) {
                stationary = true;
                break;
            }

            double step = 1.0;
            for (int backtrack = 0; backtrack < max_backtracking_steps; backtrack++) {
                VectorXd gamma_candidate = gamma + step * delta;
                if (gamma_candidate.lpNorm<Eigen::Infinity>() > NumericalConstants::PARAM_MAX) {
                    step *= backtracking_factor;
                    continue;
                }
                VectorXd g_candidate = g_vec + step * z_vec;
                const double objective_candidate = mean_grad_psi.dot(gamma_candidate) +
                    source_loss_value(g_candidate) + l1_penalty(gamma_candidate, lambda);
                const double sufficient = objective + armijo_constant * step * model_decrease +
                    1e-12 * (1.0 + std::abs(objective));
                if (std::isfinite(objective_candidate) && objective_candidate <= sufficient) {
                    gamma = gamma_candidate;
                    g_vec = g_candidate;
                    objective = objective_candidate;
                    accepted = true;
                    break;
                }
                step *= backtracking_factor;
            }

            if (!accepted) {
                if (damping >= damping_maximum) break;
                damping = std::min(damping * damping_growth, damping_maximum);
            } else if (step == 1.0) {
                damping = std::max(damping * damping_decay, damping_minimum);
            }
        }

        result.iterations = iter + 1;
        result.max_update = (gamma - gamma_old).lpNorm<Eigen::Infinity>();
        result.convergence_threshold = convergence_threshold_cpp(gamma, convergence_tol);
        if (stationary) {
            result.converged = density_ratio_score_residual(gamma, X_train, psi_prime_train,
                source_scale, mean_grad_psi, lambda, M_tau) <= convergence_tol;
            if (result.converged || inner_tol <= inner_tol_final) break;
            inner_tol = inner_tol_final;
            continue;
        }
        if (!accepted) {
            result.line_search_failures++;
            break;
        }
        // Reaching the coefficient bound is final (see below); stop iterating.
        if (gamma.lpNorm<Eigen::Infinity>() >= NumericalConstants::PARAM_MAX * bound_fraction) {
            break;
        }
        if (check_convergence_cpp(gamma_old, gamma, convergence_tol) &&
                density_ratio_score_residual(gamma, X_train, psi_prime_train,
                    source_scale, mean_grad_psi, lambda, M_tau) <= convergence_tol) {
            result.converged = true;
            break;
        }
        inner_tol = std::min(
            inner_tol_initial,
            std::max(inner_tol_final, inner_tol_step_fraction * result.max_update)
        );
    }

    // A solution pinned at the coefficient bound is the exponential-tilting
    // analogue of logistic separation: the truncated loss is linear beyond the
    // radius, so the penalized objective keeps decreasing toward the bound.
    if (gamma.lpNorm<Eigen::Infinity>() >= NumericalConstants::PARAM_MAX * bound_fraction) {
        result.converged = false;
    }

    if (use_active_set) {
        for (int j = 1; j < p; j++) {
            active[j] = std::abs(gamma(j)) > NumericalConstants::ACTIVE_SET_THRESHOLD;
        }
    }
    return result;
}

inline DensityRatioFitResult density_ratio_penalized_fit(
        VectorXd& gamma, std::vector<bool>& active,
        const MatrixXd& X_train, const VectorXd& psi_prime_train,
        double source_scale, const VectorXd& mean_grad_psi, double lambda,
        double M_tau, double convergence_tol, int max_iter,
        bool use_active_set = true) {
    DensityRatioFitResult result;
    if (nuisance_use_proximal_newton()) {
        result = density_ratio_proximal_newton(
            gamma, active, X_train, psi_prime_train, source_scale, mean_grad_psi,
            lambda, M_tau, convergence_tol, max_iter, use_active_set
        );
    } else {
        result = density_ratio_coordinate_descent(
            gamma, active, X_train, psi_prime_train, source_scale, mean_grad_psi,
            lambda, M_tau, convergence_tol, max_iter, use_active_set);
    }
    result.kkt_residual = density_ratio_score_residual(gamma, X_train, psi_prime_train,
        source_scale, mean_grad_psi, lambda, M_tau);
    result.kkt_threshold = convergence_tol;
    result.converged = result.converged && result.kkt_residual <= result.kkt_threshold;
    return result;
}

// CV-only acceleration. Certified candidates follow the caller's ordinary
// failure path, including warm-start restoration and failure-patience counting.
// A failed fit can supply a certificate without becoming an optimization start.
inline DensityRatioFitResult density_ratio_cv_fit(
        VectorXd& gamma, std::vector<bool>& active,
        const MatrixXd& X_train, const VectorXd& psi_prime_train,
        double source_scale, const VectorXd& moment, double lambda,
        double radius, double tolerance, int max_iter,
        bool use_certificate, double& certified_bound, int& certified_failures) {
    if (use_certificate && lambda < certified_bound) {
        ++certified_failures;
        return DensityRatioFitResult();
    }
    DensityRatioFitResult result = density_ratio_penalized_fit(
        gamma, active, X_train, psi_prime_train, source_scale, moment,
        lambda, radius, tolerance, max_iter);
    if (use_certificate) {
        certified_bound = std::max(certified_bound,
            density_cv_kkt_bound(gamma, X_train, moment, tolerance));
    }
    return result;
}

// Density ratio validation loss: grad^T γ + source_scale * Ẽ_val[exp(-φ^T γ) ψ']
//
// source_scale = n_treated / n_source rescales the treated-arm validation average
// (Σ_{treated in val} ... / n_val) into the full-source empirical expectation
// Ẽ_{s_j}[I(A=1) exp(-φ^T γ) ψ'] = (n_treated / n_source) * mean_{treated}[...] used by
// the training objective (density_ratio_penalized_fit uses the same arm fraction times
// a treated-fold mean) and by main.tex (eq:gamma_init,
// eq:gamma_calibrated_loss). Without this factor the
// validation loss is a treated-arm average that over-weights the exp term by
// 1/arm_fraction relative to the linear grad^T γ term, biasing λ_γ selection and
// leaving the estimated tilt γ̂ mis-regularized (degrades C2 coverage at large n).
inline double density_ratio_val_loss(const VectorXd& gamma, const VectorXd& mean_grad_psi,
                                     const MatrixXd& X_val, const VectorXd& psi_prime_val,
                                     double source_scale, double M_tau) {
    int n_val = X_val.rows();
    long double linear = static_cast<long double>(mean_grad_psi.dot(gamma));
    long double source_sum = 0.0L;
    for (int j = 0; j < n_val; j++) {
        source_sum += static_cast<long double>(
            tilt_loss(X_val.row(j).dot(gamma), M_tau) * psi_prime_val(j));
    }
    return static_cast<double>(
        linear + static_cast<long double>(source_scale) * source_sum / n_val);
}

// GLM validation loss: Σ w_j * loss(y_j, μ_j) / n_val
inline double glm_val_loss(const VectorXd& beta, const MatrixXd& X_val, const VectorXd& Y_val,
                           const VectorXd& weights_val, GLMFamily family, LinkFunction link) {
    int n_val = X_val.rows();
    StableAccumulator loss_acc;
    for (int j = 0; j < n_val; j++) {
        double eta_val = X_val.row(j).dot(beta);
        double mu_val = GLMUtils::response_function(eta_val, link);
        double loss_j = GLMUtils::calculate_loss(Y_val(j), mu_val, family);
        loss_acc.add(weights_val(j) * loss_j);
    }
    return loss_acc.value() / n_val;
}

} // namespace CVUtils

#endif // CV_UTILS_HPP
