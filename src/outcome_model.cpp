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

#include "optimization.hpp"
#include "cv_utils.hpp"

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
    
    int n = W_outcome.rows();
    int p_outcome = W_outcome.cols() + 1; // +1 for intercept
    
    if (Z_site.rows() != n || Z_site.cols() <= 0) {
        throw std::runtime_error("fit_unified_outcome_cpp: Z_site must be non-empty and match W_outcome rows.");
    }
    
    // Filter to treatment group only
    std::vector<int> treated_idx;
    for (int i = 0; i < n; i++) {
        if (A_source(i) == A_val) {
            treated_idx.push_back(i);
        }
    }
    
    if (treated_idx.empty()) {
        return List::create(Named("alpha") = VectorXd::Zero(p_outcome),
                           Named("converged") = false,
                           Named("iterations") = 0,
                           Named("family") = family_int,
                           Named("link") = link_int);
    }
    
    int n_treated = treated_idx.size();
    MatrixXd X_treated(n_treated, W_outcome.cols());
    VectorXd Y_treated(n_treated);
    
    for (int i = 0; i < n_treated; i++) {
        X_treated.row(i) = W_outcome.row(treated_idx[i]);
        Y_treated(i) = Y_source(treated_idx[i]);
    }
    
    // Calculate density ratio weights for treated units
    // Following main.tex: w(X_i; γ̂) = I(A=1)/exp(g(Z_i; γ̂))
    // CRITICAL: Use Z_site (not W_outcome) for density ratio calculation!
    MatrixXd Z_treated_for_dr = prepend_intercept(subset_rows(Z_site, treated_idx));
    VectorXd weights = compute_density_ratio_weights(Z_treated_for_dr, gamma_s, calibrated, M_tau);
    
    // Use raw (unnormalized) weights to match the paper's empirical expectation
    // Ẽ_{s_j}[w(X;γ̂) loss(Y, ψ(φ^Tα))], where the 1/n normalization is handled
    // inside fit_general_glm_cpp. Normalizing weights would change the effective
    // regularization strength relative to the data-fit term.
    // Only guard against degenerate case where all weights are zero.
    double weight_sum = weights.sum();
    if (weight_sum <= 0) {
        weights.setOnes();
    }
    
    // Use general GLM fitting with unnormalized weights
    List result = fit_general_glm_cpp(X_treated, Y_treated, weights, family_int, link_int, lambda, max_iter, tol, warm_start);
    
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
        
        // Check convergence
        if (check_convergence_cpp(beta_old, beta, tol)) {
            return List::create(Named("alpha") = beta,
                               Named("converged") = true,
                               Named("iterations") = iter + 1,
                               Named("family") = family_int,
                               Named("link") = link_int);
        }
    }
    
    return List::create(Named("alpha") = beta,
                       Named("converged") = false,
                       Named("iterations") = max_iter,
                       Named("family") = family_int,
                       Named("link") = link_int);
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
                                                 const MatrixXd& Z_site) {
    
    int n = W_outcome.rows();
    int p = W_outcome.cols() + 1;
    int n_lambda = lambda_grid.size();

    GLMFamily family = static_cast<GLMFamily>(family_int);
    LinkFunction link = static_cast<LinkFunction>(link_int);

    if (Z_site.rows() != n || Z_site.cols() <= 0) {
        throw std::runtime_error("select_lambda_cv_general_refined_outcome_cpp: Z_site must be non-empty and match W_outcome rows.");
    }

    auto lambda_order = CVUtils::sort_lambda_descending(lambda_grid);
    auto treated_idx = CVUtils::filter_treated(A_source, A_val);
    if (treated_idx.empty()) return CVUtils::make_early_return(lambda_grid, n_lambda, true);
    int n_treated = treated_idx.size();
    if (n_treated < n_folds) return CVUtils::make_early_return(lambda_grid, n_lambda, false);

    // Pre-compute outcome features with intercept
    MatrixXd X_treated_int = prepend_intercept(subset_rows(W_outcome, treated_idx));
    VectorXd Y_treated = subset_elements(Y_source, treated_idx);

    // Pre-compute weights from density ratio (refined = not calibrated)
    MatrixXd Z_treated_for_dr = prepend_intercept(subset_rows(Z_site, treated_idx));
    VectorXd weights_all = compute_density_ratio_weights(Z_treated_for_dr, gamma_s, false);

    auto folds = CVUtils::create_fold_splits(n_treated, n_folds);

    // Pre-allocate fold data
    std::vector<MatrixXd> X_train_folds(n_folds), X_val_folds(n_folds);
    std::vector<VectorXd> Y_train_folds(n_folds), Y_val_folds(n_folds);
    std::vector<VectorXd> w_train_folds(n_folds), w_val_folds(n_folds);
    for (int fold = 0; fold < n_folds; fold++) {
        X_train_folds[fold] = CVUtils::slice_rows(X_treated_int, folds.train[fold]);
        Y_train_folds[fold] = CVUtils::slice_elements(Y_treated, folds.train[fold]);
        // Use raw (unnormalized) weights to match fit_unified_outcome_cpp / fit_general_glm_cpp.
        // Normalizing here would shift the effective λ relative to the final fit.
        w_train_folds[fold] = CVUtils::slice_elements(weights_all, folds.train[fold]);
        X_val_folds[fold] = CVUtils::slice_rows(X_treated_int, folds.val[fold]);
        Y_val_folds[fold] = CVUtils::slice_elements(Y_treated, folds.val[fold]);
        w_val_folds[fold] = CVUtils::slice_elements(weights_all, folds.val[fold]);
    }

    // CV loop: outer=fold, inner=lambda
    MatrixXd fold_scores(n_lambda, n_folds);
    fold_scores.fill(INFINITY);

    for (int fold = 0; fold < n_folds; fold++) {
        if (folds.train[fold].empty() || folds.val[fold].empty()) continue;

        int n_train = folds.train[fold].size();
        VectorXd beta = VectorXd::Zero(p);
        std::vector<bool> active(p, true);

        for (int li = 0; li < n_lambda; li++) {
            int lambda_idx = lambda_order[li];
            double lambda = lambda_grid(lambda_idx);
            double cv_tol = std::max(tol, NumericalConstants::CV_TOL_FLOOR);
            int cv_max_iter = std::min(max_iter, NumericalConstants::CV_MAX_ITER);

            CVUtils::glm_cd_update(beta, active, X_train_folds[fold], Y_train_folds[fold],
                                   w_train_folds[fold], n_train, lambda, link, family, cv_tol, cv_max_iter);

            fold_scores(lambda_idx, fold) = CVUtils::glm_val_loss(
                beta, X_val_folds[fold], Y_val_folds[fold], w_val_folds[fold], family, link);
        }
    }

    return CVUtils::aggregate_cv_results(fold_scores, lambda_grid, n_lambda, n_folds);
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
                                             int family_int, int link_int) {
    
    int n = W_outcome.rows();
    int n_lambda = lambda_grid.size();
    int p = W_outcome.cols() + 1;

    if (Z_site.rows() != n || Z_site.cols() <= 0) {
        throw std::runtime_error("select_lambda_cv_calibrated_outcome_cpp: Z_site must be non-empty and match W_outcome rows.");
    }

    auto lambda_order = CVUtils::sort_lambda_descending(lambda_grid);
    auto treated_idx = CVUtils::filter_treated(A_source, A_val);
    if (treated_idx.empty()) return CVUtils::make_early_return(lambda_grid, n_lambda, true);
    int n_treated = treated_idx.size();
    if (n_treated < n_folds) return CVUtils::make_early_return(lambda_grid, n_lambda, false);

    // Pre-compute outcome features with intercept
    MatrixXd X_treated_int = prepend_intercept(subset_rows(W_outcome, treated_idx));
    VectorXd Y_treated = subset_elements(Y_source, treated_idx);

    // Pre-compute truncated weights (calibrated = true)
    MatrixXd Z_treated_for_dr = prepend_intercept(subset_rows(Z_site, treated_idx));
    VectorXd weights_all = compute_density_ratio_weights(Z_treated_for_dr, gamma_init, true, M_tau);

    auto folds = CVUtils::create_fold_splits(n_treated, n_folds);

    // Pre-allocate fold data
    std::vector<MatrixXd> X_train_folds(n_folds), X_val_folds(n_folds);
    std::vector<VectorXd> Y_train_folds(n_folds), Y_val_folds(n_folds);
    std::vector<VectorXd> w_train_folds(n_folds), w_val_folds(n_folds);
    for (int fold = 0; fold < n_folds; fold++) {
        X_train_folds[fold] = CVUtils::slice_rows(X_treated_int, folds.train[fold]);
        Y_train_folds[fold] = CVUtils::slice_elements(Y_treated, folds.train[fold]);
        // Use raw (unnormalized) weights to match fit_unified_outcome_cpp / fit_general_glm_cpp.
        w_train_folds[fold] = CVUtils::slice_elements(weights_all, folds.train[fold]);
        X_val_folds[fold] = CVUtils::slice_rows(X_treated_int, folds.val[fold]);
        Y_val_folds[fold] = CVUtils::slice_elements(Y_treated, folds.val[fold]);
        w_val_folds[fold] = CVUtils::slice_elements(weights_all, folds.val[fold]);
    }

    // CV loop — use specified GLM family/link for calibrated outcome
    GLMFamily family = static_cast<GLMFamily>(family_int);
    LinkFunction link = static_cast<LinkFunction>(link_int);
    MatrixXd fold_scores(n_lambda, n_folds);
    fold_scores.fill(INFINITY);

    for (int fold = 0; fold < n_folds; fold++) {
        if (folds.train[fold].empty() || folds.val[fold].empty()) continue;

        int n_train = folds.train[fold].size();
        VectorXd alpha = VectorXd::Zero(p);
        std::vector<bool> active(p, true);

        for (int li = 0; li < n_lambda; li++) {
            int lambda_idx = lambda_order[li];
            double lambda = lambda_grid(lambda_idx);
            double cv_tol = std::max(tol, NumericalConstants::CV_TOL_FLOOR);
            int cv_max_iter = std::min(max_iter, NumericalConstants::CV_MAX_ITER);

            CVUtils::glm_cd_update(alpha, active, X_train_folds[fold], Y_train_folds[fold],
                                   w_train_folds[fold], n_train, lambda, link, family, cv_tol, cv_max_iter);

            fold_scores(lambda_idx, fold) = CVUtils::glm_val_loss(
                alpha, X_val_folds[fold], Y_val_folds[fold], w_val_folds[fold], family, link);
        }
    }

    return CVUtils::aggregate_cv_results(fold_scores, lambda_grid, n_lambda, n_folds);
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
