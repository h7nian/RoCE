#ifndef CV_UTILS_HPP
#define CV_UTILS_HPP

#include "utils.hpp"
#include "numerical_constants.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>  // for std::iota
#include <stdexcept>

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

// Truncation function for numerical stability in two-layer cross-fitting
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
        bool calibrated = false, double M_tau = NumericalConstants::M_TAU_DEFAULT) {
    int n = Z_site_int.rows();
    VectorXd g_vals = Z_site_int * gamma_s;

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
    return weights.cwiseMax(NumericalConstants::WEIGHT_MIN)
                  .cwiseMin(NumericalConstants::WEIGHT_MAX);
}

// Adaptive step size computation for coordinate descent
// Uses gradient magnitude to scale step size, with Lipschitz constant estimation
// step_size = base_step / (1 + scale * |grad|)
// This ensures stability for large gradients while maintaining reasonable speed for small gradients
inline double compute_adaptive_step_size(double grad_magnitude, int iter, 
                                          double base_step = 0.1, double min_step = 0.001) {
    // Decay based on iteration for convergence guarantee
    double iter_decay = 1.0 / std::sqrt(iter + 1);
    
    // Gradient-based scaling for stability
    // Larger gradients -> smaller steps (prevents overshooting)
    double grad_scale = 1.0 / (1.0 + 0.5 * grad_magnitude);
    
    // Combined step size with minimum floor
    double step = base_step * iter_decay * grad_scale;
    return std::max(step, min_step);
}

// ============================================================================
// Internal CV Utilities (shared infrastructure for all 4 CV functions)
// ============================================================================
// Eliminates ~600 lines of duplicated scaffolding code by factoring out:
//   - Lambda path sorting for warm start
//   - Treated unit filtering and early-return logic
//   - Fold creation and train/val index splitting
//   - Matrix/vector row slicing for fold data
//   - GLM coordinate descent inner loop (outcome models)
//   - Density ratio coordinate descent inner loop (DR models)
//   - Score aggregation and best-lambda selection
// ============================================================================

namespace CVUtils {

// Sort lambda indices from LARGE to SMALL for effective warm starting
inline std::vector<int> sort_lambda_descending(const VectorXd& lambda_grid) {
    int n = lambda_grid.size();
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&lambda_grid](int a, int b) { return lambda_grid(a) > lambda_grid(b); });
    return order;
}

// Filter indices where A == a_val (delegates to filter_treated_indices in utils.hpp)
inline std::vector<int> filter_treated(const VectorXd& A, int a_val = 1) {
    return filter_treated_indices(A, a_val);
}

// Fold train/val index split
struct FoldSplit {
    std::vector<std::vector<int>> train;
    std::vector<std::vector<int>> val;
};

inline FoldSplit create_fold_splits(int n_samples, int n_folds) {
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

// Slice rows from a matrix by index vector (delegates to subset_rows in utils.hpp)
inline MatrixXd slice_rows(const MatrixXd& M, const std::vector<int>& indices) {
    return subset_rows(M, indices);
}

// Slice elements from a vector by index vector (delegates to subset_elements in utils.hpp)
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
inline List aggregate_cv_results(const MatrixXd& fold_scores, const VectorXd& lambda_grid,
                                 int n_lambda, int n_folds) {
    VectorXd cv_scores(n_lambda);
    VectorXd cv_se(n_lambda);
    IntegerVector n_valid(n_lambda);
    for (int i = 0; i < n_lambda; i++) {
        StableAccumulator sum_acc;
        for (int fold = 0; fold < n_folds; fold++) {
            if (!std::isfinite(fold_scores(i, fold))) {
                throw std::runtime_error(
                    "aggregate_cv_results: non-finite validation loss at lambda index " +
                    std::to_string(i + 1) + ", fold " + std::to_string(fold + 1) +
                    "; CV lambda selection requires every validation fold to be finite."
                );
            }
            sum_acc.add(fold_scores(i, fold));
        }
        n_valid[i] = n_folds;
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
        throw std::runtime_error("aggregate_cv_results: all lambda values produced invalid CV scores.");
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

    return List::create(Named("best_lambda") = lambda_grid(best_idx),
                       Named("best_idx") = best_idx + 1,
                       Named("lambda_min") = lambda_grid(best_idx),
                       Named("lambda_1se") = lambda_grid(idx_1se),
                       Named("idx_min") = best_idx + 1,
                       Named("idx_1se") = idx_1se + 1,
                       Named("cv_scores") = cv_scores,
                       Named("cv_se") = cv_se,
                       Named("n_valid_folds") = n_valid,
                       Named("lambda_rule") = "min");
}

// ---- Coordinate Descent Inner Loops ----

// L1-penalized GLM coordinate descent (shared by refined & calibrated outcome CV)
inline void glm_cd_update(VectorXd& beta, std::vector<bool>& active,
                          const MatrixXd& X_train, const VectorXd& Y_train,
                          const VectorXd& weights_train, int n_train,
                          double lambda, LinkFunction link, GLMFamily family,
                          double cv_tol, int cv_max_iter) {
    int p = beta.size();

    for (int iter = 0; iter < cv_max_iter; iter++) {
        VectorXd beta_old = beta;

        if (iter % NumericalConstants::ACTIVE_CHECK_FREQ == 0) {
            for (int j = 1; j < p; j++) {
                active[j] = (std::abs(beta(j)) > NumericalConstants::ACTIVE_SET_THRESHOLD) || (iter == 0);
            }
        }

        // Compute eta = X_train * beta once per outer iteration (O(np)),
        // then update incrementally after each coordinate change (O(n) per coordinate).
        // This reduces total complexity from O(np²) to O(np) per iteration.
        VectorXd eta = X_train * beta;

        for (int j = 0; j < p; j++) {
            if (j > 0 && !active[j] && iter % NumericalConstants::ACTIVE_CHECK_FREQ != 0) continue;

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
            }
        }

        if (check_convergence_cpp(beta_old, beta, cv_tol)) break;
    }
}

// L1-penalized density ratio coordinate descent (shared by refined & calibrated DR CV)
// Optimized with incremental g = X*gamma update (same trick as glm_cd_update),
// reducing per-iteration complexity from O(np²) to O(np).
inline void density_ratio_cd_update(VectorXd& gamma, std::vector<bool>& active,
                                    const MatrixXd& X_train, const VectorXd& psi_prime_train,
                                    int n_total, const VectorXd& mean_grad_psi, double lambda,
                                    double cv_tol, int cv_max_iter) {
    int p = gamma.size();
    int n_train = X_train.rows();

    if (mean_grad_psi.size() != p) {
        throw std::runtime_error(
            "density_ratio_cd_update: mean_grad_psi length must match gamma length."
        );
    }

    for (int iter = 0; iter < cv_max_iter; iter++) {
        VectorXd gamma_old = gamma;

        if (iter % NumericalConstants::ACTIVE_CHECK_FREQ == 0) {
            for (int j = 1; j < p; j++) {
                active[j] = (std::abs(gamma(j)) > NumericalConstants::ACTIVE_SET_THRESHOLD) || (iter == 0);
            }
        }

        // Compute g = X_train * gamma once per outer iteration (O(np)),
        // then update incrementally after each coordinate change (O(n) per coordinate).
        // This reduces total complexity from O(np²) to O(np) per iteration.
        VectorXd g_vec = X_train * gamma;

        for (int j = 0; j < p; j++) {
            if (j > 0 && !active[j] && iter % NumericalConstants::ACTIVE_CHECK_FREQ != 0) continue;

            // Compute gradient and Hessian-diagonal approximation for coordinate j
            // grad_j = mean_grad_psi(j) - E_s[ X_j * exp(-φ^T γ) * ψ' ] / n_total
            // hess_j = E_s[ X_j^2 * exp(-φ^T γ) * ψ' ] / n_total  (used for Newton-like step)
            double grad_j = mean_grad_psi(j);
            StableAccumulator grad_acc;
            StableAccumulator hess_acc;

            for (int i = 0; i < n_train; i++) {
                double current_g = std::max(NumericalConstants::ETA_CLIP_MIN, std::min(NumericalConstants::ETA_CLIP_MAX, g_vec(i)));
                double exp_neg_g = std::exp(-current_g);
                exp_neg_g = std::max(NumericalConstants::WEIGHT_MIN, std::min(NumericalConstants::WEIGHT_MAX, exp_neg_g));
                double psi_p = psi_prime_train(i);

                // first derivative (gradient) contribution
                grad_acc.add(-X_train(i, j) * exp_neg_g * psi_p / n_total);
                // second derivative (Hessian diagonal) contribution
                hess_acc.add(X_train(i, j) * X_train(i, j) * exp_neg_g * psi_p / n_total);
            }
            grad_j += grad_acc.value();
            double hess_j = hess_acc.value();

            // Hessian-based step size (with floor for stability)
            hess_j = std::max(hess_j, NumericalConstants::HESSIAN_FLOOR);
            double step_size_hessian = 1.0 / hess_j;

            // Combine with adaptive step as an upper bound for safety (prevents huge steps)
            double step_size_adaptive = compute_adaptive_step_size(std::abs(grad_j), iter);
            double step_size = std::min(step_size_hessian, step_size_adaptive * 10.0);
            step_size = std::max(step_size, NumericalConstants::STEP_SIZE_FLOOR);

            // Proximal gradient update with soft-thresholding (no penalty on intercept)
            double old_gamma_j = gamma(j);
            gamma(j) = sanitize_param(proximal_update(gamma(j), grad_j, step_size, lambda, j));

            // Incremental g update: O(n) instead of recomputing O(np)
            double delta_j = gamma(j) - old_gamma_j;
            if (delta_j != 0.0) {
                g_vec += delta_j * X_train.col(j);
            }
        }

        if (check_convergence_cpp(gamma_old, gamma, cv_tol)) break;
    }
}

// Density ratio validation loss: grad^T γ + source_scale * Ẽ_val[exp(-φ^T γ) ψ']
//
// source_scale = n_treated / n_source rescales the treated-arm validation average
// (Σ_{treated in val} ... / n_val) into the full-source empirical expectation
// Ẽ_{s_j}[I(A=1) exp(-φ^T γ) ψ'] = (n_treated / n_source) * mean_{treated}[...] used by
// the training objective (density_ratio_cd_update normalizes by the full source size n)
// and by main.tex (eq:gamma_init, eq:gamma_calibrated_loss). Without this factor the
// validation loss is a treated-arm average that over-weights the exp term by
// 1/arm_fraction relative to the linear grad^T γ term, biasing λ_γ selection and
// leaving the estimated tilt γ̂ mis-regularized (degrades C2 coverage at large n).
inline double density_ratio_val_loss(const VectorXd& gamma, const VectorXd& mean_grad_psi,
                                     const MatrixXd& X_val, const VectorXd& psi_prime_val,
                                     double source_scale) {
    int n_val = X_val.rows();
    long double linear = static_cast<long double>(mean_grad_psi.dot(gamma));
    long double source_sum = 0.0L;
    for (int j = 0; j < n_val; j++) {
        double eta_gamma = X_val.row(j).dot(gamma);
        eta_gamma = std::max(NumericalConstants::ETA_CLIP_MIN, std::min(NumericalConstants::ETA_CLIP_MAX, eta_gamma));
        double exp_neg_g = std::exp(-eta_gamma);
        exp_neg_g = std::max(NumericalConstants::WEIGHT_MIN, std::min(NumericalConstants::WEIGHT_MAX, exp_neg_g));
        source_sum += static_cast<long double>(exp_neg_g * psi_prime_val(j));
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
