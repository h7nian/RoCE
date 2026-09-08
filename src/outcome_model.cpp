// ============================================================================
// outcome_model.cpp — Outcome model (α) fitting, GLM utilities, and CV
// ============================================================================
// Functions:
//   - fit_unified_outcome_cpp:                         Refined/Calibrated α fitting
//   - fit_general_glm_cpp:                             General GLM fitting (internal)
//   - select_lambda_cv_general_refined_outcome_cpp:    CV for refined α
//   - select_lambda_cv_calibrated_outcome_cpp:         CV for calibrated α
//   - calculate_glm_gradient_cpp:                      GLM gradient computation
//   - predict_glm_cpp:                                 GLM prediction (vectorized)
// ============================================================================

#include "optimization.h"
#include "cv_utils.h"

namespace {

// Return the same stopping diagnostics for outcome fits that are already
// exposed by the density-ratio solvers.  Keeping these in the fit object lets
// the R layer distinguish a well-converged solution from a finite parameter
// vector that stopped at the iteration or coefficient boundary.
List outcome_fit_result(const VectorXd& alpha, bool converged,
                        int iterations, double max_update,
                        double convergence_threshold,
                        int family_int, int link_int) {
    return List::create(
        Named("alpha") = alpha,
        Named("converged") = converged,
        Named("iterations") = iterations,
        Named("max_update") = max_update,
        Named("convergence_threshold") = convergence_threshold,
        // The current GLM coordinate update has no backtracking stage, so a
        // line-search failure is structurally impossible rather than missing.
        Named("line_search_failures") = 0,
        Named("max_abs_coefficient") = alpha.lpNorm<Eigen::Infinity>(),
        Named("family") = family_int,
        Named("link") = link_int
    );
}

}  // namespace

// Forward declaration — fit_general_glm_cpp is called by fit_unified_outcome_cpp
List fit_general_glm_cpp(const MatrixXd& X, const VectorXd& Y, const VectorXd& weights,
                        int family_int, int link_int, double lambda, 
                        int max_iter, double tol,
                        const VectorXd& warm_start);

// UNIFIED Refined/Calibrated loss function for outcome model (α) parameters
// Following main.tex:
//   - Untruncated counterpart of eq:alpha_calibrated_loss: Refined loss (calibrated = false)
//   - eq:alpha_calibrated_loss in main.tex: Calibrated loss with truncation \mathcal{T}(·) (calibrated = true)
// Supports multiple GLM families via family_int and link_int parameters
// IMPORTANT: Z_site is used for density ratio weights, W_outcome for outcome model
// [[Rcpp::export]]
List fit_unified_outcome_cpp(const MatrixXd& W_outcome, const VectorXd& Y_source, 
                             const VectorXd& A_source, const VectorXd& gamma_s,
                             int family_int, int link_int, double lambda, 
                             int max_iter, double tol, int A_val,
                             bool calibrated, double M_tau,
                             const MatrixXd& Z_site,
                             const VectorXd& warm_start) {
    
    if (A_val != 0 && A_val != 1) {
        throw std::runtime_error("fit_unified_outcome_cpp: A_val must be 0 or 1.");
    }

    int n = W_outcome.rows();
    int p_outcome = W_outcome.cols() + 1; // +1 for intercept
    
    if (Z_site.rows() != n || Z_site.cols() <= 0) {
        throw std::runtime_error("fit_unified_outcome_cpp: Z_site must be non-empty and match W_outcome rows.");
    }
    
    // Filter to the requested treatment arm A == A_val. This may be treatment
    // or control, so use arm_* names instead of treated_*.
    std::vector<int> arm_idx;
    for (int i = 0; i < n; i++) {
        if (A_source(i) == A_val) {
            arm_idx.push_back(i);
        }
    }
    
    if (arm_idx.empty()) {
        throw std::runtime_error("fit_unified_outcome_cpp: no observations with A_val in this fold.");
    }
    
    int n_arm = arm_idx.size();
    MatrixXd X_arm(n_arm, W_outcome.cols());
    VectorXd Y_arm(n_arm);
    
    for (int i = 0; i < n_arm; i++) {
        X_arm.row(i) = W_outcome.row(arm_idx[i]);
        Y_arm(i) = Y_source(arm_idx[i]);
    }
    
    // Calculate density ratio weights for units in the requested arm
    // Following main.tex: w_a(X_i; γ̂) = I(A=a)/exp(g(Z_i; γ̂))
    // CRITICAL: Use Z_site (not W_outcome) for density ratio calculation!
    MatrixXd Z_arm_for_dr = prepend_intercept(subset_rows(Z_site, arm_idx));
    VectorXd weights = compute_density_ratio_weights(Z_arm_for_dr, gamma_s, calibrated, M_tau);
    
    // fit_general_glm_cpp averages over the arm-only matrix passed below.
    // Scale weights by n_arm / n_source so the data-fit term equals the
    // paper's full-source expectation:
    //   (1/n_source) Σ I(A=a) w_i loss_i.
    // Without this factor, the same lambda would be too weak whenever the
    // treatment arm is a fraction of the source sample.
    // Only guard against degenerate case where all weights are zero.
    double weight_sum = weights.sum();
    if (weight_sum <= 0) {
        throw std::runtime_error("fit_unified_outcome_cpp: density-ratio weights sum to zero; cannot fit method-aligned weighted outcome loss.");
    }
    weights *= static_cast<double>(n_arm) / static_cast<double>(n);
    
    // Use general GLM fitting with unnormalized weights
    List result = fit_general_glm_cpp(X_arm, Y_arm, weights, family_int, link_int, lambda, max_iter, tol, warm_start);
    
    return result;
}


// General GLM fitting function - called internally by fit_unified_outcome_cpp
// Follows the untruncated refined counterpart of eq:alpha_calibrated_loss in main.tex
// [[Rcpp::export]]
List fit_general_glm_cpp(const MatrixXd& X, const VectorXd& Y, const VectorXd& weights,
                        int family_int, int link_int, double lambda, 
                        int max_iter, double tol,
                        const VectorXd& warm_start) {
    
    GLMFamily family = static_cast<GLMFamily>(family_int);
    LinkFunction link = static_cast<LinkFunction>(link_int);
    
    int n = X.rows();
    MatrixXd X_int = prepend_intercept(X);
    int p = X_int.cols();
    
    // Initialize beta: warm-start from previous solution if provided
    VectorXd beta = (warm_start.size() == p) ? warm_start : VectorXd::Zero(p);
    
    bool converged = false;
    int iterations = 0;
    double max_update = NA_REAL;
    double convergence_threshold = NA_REAL;

    // Coordinate descent for the refined (untruncated) weighted GLM objective
    for (int iter = 0; iter < max_iter; iter++) {
        VectorXd beta_old = beta;
        
        // Compute eta = X_int * beta once per outer iteration (O(np)),
        // then update incrementally after each coordinate change (O(n) per coordinate).
        // This reduces total complexity from O(np²) to O(np) per iteration.
        VectorXd eta = X_int * beta;
        
        for (int j = 0; j < p; j++) {
            // Calculate weighted gradient and Hessian for coordinate j
            StableAccumulator grad_acc;
            StableAccumulator hess_acc;
            
            for (int i = 0; i < n; i++) {
                double x_ij = X_int(i, j);
                double w_i = weights(i);
                double eta_i = eta(i);
                double y_i = Y(i);
                
                // Refined weighted GLM objective: ℓ(β) = Ẽ_s[w(X;α̂) NLL(Y, h(φ^T β))] + λ||β||_1
                // Gradient w.r.t. β_j: Ẽ_s[w(X;α̂) · r_i · φ_j] where r_i = ∂NLL/∂η = -(y-μ)h'(η)/V(μ)
                double mu_i = GLMUtils::response_function(eta_i, link);
                
                // NLL gradient residual: (μ-y) for canonical links, (y-μ) for inverse link
                double r_i = GLMUtils::nll_gradient_residual(y_i, mu_i, link);
                grad_acc.add(w_i * r_i * x_ij);
                
                // Weighted Hessian diagonal: IRLS working weight h'(η)² / V(μ)
                // For canonical links this equals h'(η); for non-canonical (inverse)
                // this ensures positive-definiteness.
                double ww_i = GLMUtils::irls_working_weight(eta_i, link, family);
                hess_acc.add(w_i * ww_i * x_ij * x_ij);
            }
            
            // Average over sample size
            double grad_j = grad_acc.value() / n;
            double hess_j = hess_acc.value() / n;
            
            // Ensure positive definite Hessian
            hess_j = std::max(hess_j, NumericalConstants::HESSIAN_FLOOR);
            
            // Newton-Raphson step size from Hessian diagonal
            double step_size = 1.0 / hess_j;
            
            // Save old value for incremental eta update
            double old_beta_j = beta(j);
            
            // Proximal gradient update with soft-thresholding (no penalty on intercept)
            beta(j) = sanitize_param(proximal_update(beta(j), grad_j, step_size, lambda, j));
            
            // Incremental eta update: O(n) instead of recomputing O(np)
            double delta_j = beta(j) - old_beta_j;
            if (delta_j != 0.0) {
                eta += delta_j * X_int.col(j);
            }
        }
        
        iterations = iter + 1;
        max_update = (beta - beta_old).lpNorm<Eigen::Infinity>();
        convergence_threshold = convergence_threshold_cpp(beta, tol);

        // Check convergence
        if (check_convergence_cpp(beta_old, beta, tol)) {
            converged = true;
            break;
        }
    }

    return outcome_fit_result(
        beta, converged, iterations, max_update, convergence_threshold,
        family_int, link_int
    );
}

// ============================================================================
// GLMNET-STYLE CV: General Refined Outcome Model (GLM)
// Optimizations: Warm Start + Pathwise Descent + Pre-computed Folds + Active Set
// ============================================================================
// [[Rcpp::export]]
List select_lambda_cv_general_refined_outcome_cpp(const MatrixXd& W_outcome, const VectorXd& Y_source,
                                                 const VectorXd& A_source, const VectorXd& gamma_s,
                                                 int family_int, int link_int, const VectorXd& lambda_grid, 
                                                 int n_folds, int max_iter, double tol, int A_val,
                                                 const MatrixXd& Z_site,
                                                 Rcpp::Nullable<Rcpp::NumericVector> cv_fold_id = R_NilValue) {
    
    int n = W_outcome.rows();
    int p = W_outcome.cols() + 1;
    int n_lambda = lambda_grid.size();

    GLMFamily family = static_cast<GLMFamily>(family_int);
    LinkFunction link = static_cast<LinkFunction>(link_int);

    if (Z_site.rows() != n || Z_site.cols() <= 0) {
        throw std::runtime_error("select_lambda_cv_general_refined_outcome_cpp: Z_site must be non-empty and match W_outcome rows.");
    }

    auto lambda_order = CVUtils::sort_lambda_descending(lambda_grid);
    auto arm_idx = CVUtils::filter_treated(A_source, A_val);
    if (arm_idx.empty()) {
        throw std::runtime_error("select_lambda_cv_general_refined_outcome_cpp: no observations with A_val.");
    }
    int n_arm = arm_idx.size();
    if (n_arm < n_folds) {
        throw std::runtime_error("select_lambda_cv_general_refined_outcome_cpp: not enough A_val observations for requested CV folds.");
    }
    auto folds = CVUtils::create_fold_splits(n_arm, n_folds, cv_fold_id);

    // Pre-compute outcome features with intercept
    MatrixXd X_arm_int = prepend_intercept(subset_rows(W_outcome, arm_idx));
    VectorXd Y_arm = subset_elements(Y_source, arm_idx);

    // Pre-compute weights from density ratio (refined = not calibrated)
    MatrixXd Z_arm_for_dr = prepend_intercept(subset_rows(Z_site, arm_idx));
    VectorXd weights_all = compute_density_ratio_weights(Z_arm_for_dr, gamma_s, false);

    // Pre-allocate fold data
    std::vector<MatrixXd> X_train_folds(n_folds), X_val_folds(n_folds);
    std::vector<VectorXd> Y_train_folds(n_folds), Y_val_folds(n_folds);
    std::vector<VectorXd> w_train_folds(n_folds), w_val_folds(n_folds);
    for (int fold = 0; fold < n_folds; fold++) {
        X_train_folds[fold] = CVUtils::slice_rows(X_arm_int, folds.train[fold]);
        Y_train_folds[fold] = CVUtils::slice_elements(Y_arm, folds.train[fold]);
        // The arm-training mean estimates the conditional arm expectation.
        // Use the full source arm fraction, as in the final refit; including
        // the CV training fraction here would change the effective lambda.
        double train_scale = static_cast<double>(n_arm) / static_cast<double>(n);
        double val_scale = static_cast<double>(folds.val[fold].size()) / static_cast<double>(n);
        w_train_folds[fold] = CVUtils::slice_elements(weights_all, folds.train[fold]) * train_scale;
        X_val_folds[fold] = CVUtils::slice_rows(X_arm_int, folds.val[fold]);
        Y_val_folds[fold] = CVUtils::slice_elements(Y_arm, folds.val[fold]);
        w_val_folds[fold] = CVUtils::slice_elements(weights_all, folds.val[fold]) * val_scale;
    }

    // CV loop: outer=fold, inner=lambda
    MatrixXd fold_scores(n_lambda, n_folds);
    fold_scores.fill(INFINITY);
    std::vector<int> path_tail_skipped(n_folds, 0);

    const int cv_threads = CVUtils::nuisance_cv_thread_count(n_folds);
    ROCE_PARALLELIZE_CV_FOLDS(cv_threads)
    for (int fold = 0; fold < n_folds; fold++) {
        if (folds.train[fold].empty() || folds.val[fold].empty()) continue;

        int n_train = folds.train[fold].size();
        VectorXd beta = VectorXd::Zero(p);
        std::vector<bool> active(p, true);
        int consecutive_failures = 0;

        for (int li = 0; li < n_lambda; li++) {
            int lambda_idx = lambda_order[li];
            double lambda = lambda_grid(lambda_idx);
            double cv_tol = std::max(tol, NumericalConstants::CV_TOL_FLOOR);
            int cv_max_iter = std::min(max_iter, NumericalConstants::CV_MAX_ITER);

            VectorXd beta_before = beta;
            CVUtils::GLMCDResult fit = CVUtils::glm_cd_update(
                beta, active, X_train_folds[fold], Y_train_folds[fold],
                w_train_folds[fold], n_train, lambda, link, family,
                cv_tol, cv_max_iter
            );
            if (fit.converged) {
                consecutive_failures = 0;
                fold_scores(lambda_idx, fold) = CVUtils::glm_val_loss(
                    beta, X_val_folds[fold], Y_val_folds[fold],
                    w_val_folds[fold], family, link
                );
            } else {
                beta = beta_before;
                std::fill(active.begin(), active.end(), true);
                path_tail_skipped[fold] = CVUtils::cv_path_tail_after_failure(
                    consecutive_failures, li, n_lambda
                );
                if (path_tail_skipped[fold] > 0) break;
            }
        }
    }

    int skipped_fold_fits = std::accumulate(
        path_tail_skipped.begin(), path_tail_skipped.end(), 0
    );
    return CVUtils::append_explicit_fold_audit(
        CVUtils::aggregate_cv_results(
            fold_scores, lambda_grid, n_lambda, n_folds, skipped_fold_fits
        ), folds, cv_fold_id
    );
}

// ============================================================================
// GLMNET-STYLE CV: Calibrated Outcome Model
// Optimizations: Warm Start + Pathwise Descent + Pre-computed Folds + Active Set
// ============================================================================
// [[Rcpp::export]]
List select_lambda_cv_calibrated_outcome_cpp(const MatrixXd& W_outcome, const VectorXd& Y_source,
                                             const VectorXd& A_source, const VectorXd& gamma_init,
                                             const VectorXd& lambda_grid, int n_folds,
                                             int max_iter, double tol, int A_val, double M_tau,
                                             const MatrixXd& Z_site,
                                             int family_int, int link_int,
                                             Rcpp::Nullable<Rcpp::NumericVector> cv_fold_id = R_NilValue) {
    
    int n = W_outcome.rows();
    int n_lambda = lambda_grid.size();
    int p = W_outcome.cols() + 1;

    if (Z_site.rows() != n || Z_site.cols() <= 0) {
        throw std::runtime_error("select_lambda_cv_calibrated_outcome_cpp: Z_site must be non-empty and match W_outcome rows.");
    }

    auto lambda_order = CVUtils::sort_lambda_descending(lambda_grid);
    auto arm_idx = CVUtils::filter_treated(A_source, A_val);
    if (arm_idx.empty()) {
        throw std::runtime_error("select_lambda_cv_calibrated_outcome_cpp: no observations with A_val.");
    }
    int n_arm = arm_idx.size();
    if (n_arm < n_folds) {
        throw std::runtime_error("select_lambda_cv_calibrated_outcome_cpp: not enough A_val observations for requested CV folds.");
    }
    auto folds = CVUtils::create_fold_splits(n_arm, n_folds, cv_fold_id);

    // Pre-compute outcome features with intercept
    MatrixXd X_arm_int = prepend_intercept(subset_rows(W_outcome, arm_idx));
    VectorXd Y_arm = subset_elements(Y_source, arm_idx);

    // Pre-compute truncated weights (calibrated = true)
    MatrixXd Z_arm_for_dr = prepend_intercept(subset_rows(Z_site, arm_idx));
    VectorXd weights_all = compute_density_ratio_weights(Z_arm_for_dr, gamma_init, true, M_tau);

    // Pre-allocate fold data
    std::vector<MatrixXd> X_train_folds(n_folds), X_val_folds(n_folds);
    std::vector<VectorXd> Y_train_folds(n_folds), Y_val_folds(n_folds);
    std::vector<VectorXd> w_train_folds(n_folds), w_val_folds(n_folds);
    for (int fold = 0; fold < n_folds; fold++) {
        X_train_folds[fold] = CVUtils::slice_rows(X_arm_int, folds.train[fold]);
        Y_train_folds[fold] = CVUtils::slice_elements(Y_arm, folds.train[fold]);
        // Match the full-refit loss scale without an extra CV training fraction.
        // Validation weights below retain observation-count weighting across folds.
        double train_scale = static_cast<double>(n_arm) / static_cast<double>(n);
        double val_scale = static_cast<double>(folds.val[fold].size()) / static_cast<double>(n);
        w_train_folds[fold] = CVUtils::slice_elements(weights_all, folds.train[fold]) * train_scale;
        X_val_folds[fold] = CVUtils::slice_rows(X_arm_int, folds.val[fold]);
        Y_val_folds[fold] = CVUtils::slice_elements(Y_arm, folds.val[fold]);
        w_val_folds[fold] = CVUtils::slice_elements(weights_all, folds.val[fold]) * val_scale;
    }

    // CV loop — use specified GLM family/link for calibrated outcome
    GLMFamily family = static_cast<GLMFamily>(family_int);
    LinkFunction link = static_cast<LinkFunction>(link_int);
    MatrixXd fold_scores(n_lambda, n_folds);
    fold_scores.fill(INFINITY);
    std::vector<int> path_tail_skipped(n_folds, 0);

    const int cv_threads = CVUtils::nuisance_cv_thread_count(n_folds);
    ROCE_PARALLELIZE_CV_FOLDS(cv_threads)
    for (int fold = 0; fold < n_folds; fold++) {
        if (folds.train[fold].empty() || folds.val[fold].empty()) continue;

        int n_train = folds.train[fold].size();
        VectorXd alpha = VectorXd::Zero(p);
        std::vector<bool> active(p, true);
        int consecutive_failures = 0;

        for (int li = 0; li < n_lambda; li++) {
            int lambda_idx = lambda_order[li];
            double lambda = lambda_grid(lambda_idx);
            double cv_tol = std::max(tol, NumericalConstants::CV_TOL_FLOOR);
            int cv_max_iter = std::min(max_iter, NumericalConstants::CV_MAX_ITER);

            VectorXd alpha_before = alpha;
            CVUtils::GLMCDResult fit = CVUtils::glm_cd_update(
                alpha, active, X_train_folds[fold], Y_train_folds[fold],
                w_train_folds[fold], n_train, lambda, link, family,
                cv_tol, cv_max_iter
            );
            if (fit.converged) {
                consecutive_failures = 0;
                fold_scores(lambda_idx, fold) = CVUtils::glm_val_loss(
                    alpha, X_val_folds[fold], Y_val_folds[fold],
                    w_val_folds[fold], family, link
                );
            } else {
                alpha = alpha_before;
                std::fill(active.begin(), active.end(), true);
                path_tail_skipped[fold] = CVUtils::cv_path_tail_after_failure(
                    consecutive_failures, li, n_lambda
                );
                if (path_tail_skipped[fold] > 0) break;
            }
        }
    }

    int skipped_fold_fits = std::accumulate(
        path_tail_skipped.begin(), path_tail_skipped.end(), 0
    );
    return CVUtils::append_explicit_fold_audit(
        CVUtils::aggregate_cv_results(
            fold_scores, lambda_grid, n_lambda, n_folds, skipped_fold_fits
        ), folds, cv_fold_id
    );
}


// GLM gradient calculation (default: BINOMIAL with LOGIT link)
// [[Rcpp::export]]
MatrixXd calculate_glm_gradient_cpp(const MatrixXd& W_outcome, const VectorXd& beta,
                                   int family_int, int link_int) {
    // family_int accepted for API consistency but only link determines gradients
    LinkFunction link = static_cast<LinkFunction>(link_int);
    
    int n = W_outcome.rows();
    
    // Build intercept-augmented matrix and compute linear predictors via matrix multiply
    MatrixXd X_int = prepend_intercept(W_outcome);
    VectorXd eta_vec = X_int * beta;
    
    MatrixXd gradients(n, X_int.cols());
    for (int i = 0; i < n; i++) {
        // Calculate response function derivative h'(η)
        double h_prime = GLMUtils::response_derivative(eta_vec(i), link);
        // ∇_β h(φ(X); β) = h'(φ(X)^T β) * φ(X)
        gradients.row(i) = h_prime * X_int.row(i);
    }
    
    return gradients;
}

// Column-means of GLM gradient — avoids n×p intermediate matrix and R round-trip.
//
// Computes  (1/n) Σᵢ h'(φ(Xᵢ)ᵀβ) · φ̃(Xᵢ)  directly as a (p+1)-vector,
// where φ̃ = [1, φ(X)] (intercept-augmented features).
//
// Called in run_crossfit (two-round and one-round modes) to compute
// mean_grad_psi_init = colMeans(∇_α ψ(W; α̂_init)) for each target fold k2.
// Replacing colMeans(calculate_glm_gradient_cpp(...)) eliminates:
//   (a) n×(p+1) matrix allocation in C++ and its transfer back to R, and
//   (b) R's colMeans dispatching over the n-row result.
// For p=200 and n_fold=200 this removes ~320 KB per call; called K×(K_f-1) times.
// [[Rcpp::export]]
VectorXd mean_glm_gradient_cpp(const MatrixXd& W_outcome, const VectorXd& beta,
                               int family_int, int link_int) {
    LinkFunction link = static_cast<LinkFunction>(link_int);

    int n = W_outcome.rows();

    // Build intercept-augmented matrix and linear predictors
    MatrixXd X_int = prepend_intercept(W_outcome);
    VectorXd eta_vec = X_int * beta;

    // Accumulate weighted column sums directly — no n×(p+1) matrix needed
    // NOTE: Uses |h'(η)| because the paper's ψ'(θ) = b''(θ) > 0 always.
    // For canonical links (identity, logit, log), h'(η) > 0 so abs is a no-op.
    // For inverse link, h'(η) = -1/η² < 0 but b''(θ) = 1/η² > 0.
    VectorXd col_sums = VectorXd::Zero(X_int.cols());
    for (int i = 0; i < n; i++) {
        double h_prime = std::abs(GLMUtils::response_derivative(eta_vec(i), link));
        col_sums += h_prime * X_int.row(i).transpose();
    }
    return col_sums / static_cast<double>(n);
}

// GLM prediction (default: BINOMIAL with LOGIT link)
// Vectorized: uses Eigen array operations for common link functions
// [[Rcpp::export]]
VectorXd predict_glm_cpp(const MatrixXd& W_outcome, const VectorXd& beta,
                        int family_int, int link_int) {
    // family_int accepted for API consistency but only link determines predictions
    LinkFunction link = static_cast<LinkFunction>(link_int);
    
    // Build intercept-augmented matrix and compute linear predictors via matrix multiply
    MatrixXd X_int = prepend_intercept(W_outcome);
    VectorXd eta_vec = X_int * beta;
    
    // Vectorized response function for common link types
    switch (link) {
        case LinkFunction::IDENTITY:
            return eta_vec;  // h(η) = η
        case LinkFunction::LOGIT: {
            // h(η) = 1/(1 + exp(-η)) — vectorized with Eigen
            Eigen::ArrayXd eta_clipped = eta_vec.array()
                .max(NumericalConstants::ETA_CLIP_MIN)
                .min(NumericalConstants::ETA_CLIP_MAX);
            return (1.0 + (-eta_clipped).exp()).inverse().matrix();
        }
        default: {
            Rcpp::stop("Unsupported GLM link code: only identity and logit are supported.");
        }
    }
}
