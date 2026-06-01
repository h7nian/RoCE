// ============================================================================
// density_ratio.cpp — Density ratio (γ) fitting and CV functions
// ============================================================================
// Functions:
//   - fit_unified_density_ratio_cpp:                   Refined/Calibrated γ fitting
//   - fit_initial_density_ratio_cpp:                   Initial γ fitting (no ψ')
//   - select_lambda_cv_density_ratio_cpp:              CV for refined γ
//   - select_lambda_cv_initial_density_ratio_cpp:      CV for initial γ
//   - select_lambda_cv_calibrated_density_ratio_cpp:   CV for calibrated γ
// ============================================================================

#include "optimization.hpp"
#include "cv_utils.hpp"

// ============================================================================
// NAMING CONVENTION NOTE (C++ ←→ R mapping)
// ============================================================================
// C++ Z_site     = R Z_site     (site model features for density ratio)
// C++ W_outcome  = R W_outcome  (outcome model features)
// C++ A_source   = R A          (treatment indicator for source site)
// C++ mean_phi   = R mean_phi   (target site mean features Ẽ_t[φ(X)])
// ============================================================================

// UNIFIED Refined/Calibrated loss function for density ratio (γ) parameters
// Following main.tex:
//   - Untruncated counterpart of eq:gamma_calibrated_loss: Refined loss (calibrated = false)
//   - eq:gamma_calibrated_loss in main.tex: Calibrated loss with truncation \mathcal{T}(·) (calibrated = true)
// 
// Parameters:
//   - Z_site: Features for site/treatment model (Z_site) - used for γ parameterization
//   - W_outcome: Features for outcome model - used for ψ' computation
//
// SAMPLE SIZE NORMALIZATION:
// The gradient contribution from source data is normalized by n (total sample size), NOT n_treated.
// This matches the empirical expectation definition in main.tex:
//   Ẽ_{s_j}[I(A=a) f(X)] = (1/n_s) * Σ_{i=1}^{n_s} I(A_i=a) * f(X_i)
// The indicator function I(A_i=a) zeros out contributions from the other arm,
// but the normalizing constant is the full sample size n_s.
// [[Rcpp::export]]
List fit_unified_density_ratio_cpp(const MatrixXd& Z_site, const VectorXd& A_source, 
                                   const VectorXd& mean_grad_psi, const VectorXd& alpha_init,
                                   double lambda, int max_iter, double tol, 
                                   bool calibrated, double M_tau,
                                   const MatrixXd& W_outcome,
                                   int A_val,
                                   int family_int, int link_int,
                                   const VectorXd& warm_start) {
    
    if (A_val != 0 && A_val != 1) {
        throw std::runtime_error("fit_unified_density_ratio_cpp: A_val must be 0 or 1.");
    }

    int n = Z_site.rows();
    int p_site = Z_site.cols() + 1; // +1 for intercept
    
    if (W_outcome.rows() != n || W_outcome.cols() <= 0) {
        throw std::runtime_error("fit_unified_density_ratio_cpp: W_outcome must be non-empty and match Z_site rows.");
    }
    int p_outcome = W_outcome.cols() + 1;
    
    if (mean_grad_psi.size() != p_site) {
        throw std::runtime_error(
            "fit_unified_density_ratio_cpp: mean_grad_psi length must match "
            "ncol(Z_site) + 1. Compute E_t[h'(W alpha) * Z_tilde] when W_outcome "
            "and Z_site differ."
        );
    }
    if (alpha_init.size() != p_outcome) {
        throw std::runtime_error("fit_unified_density_ratio_cpp: alpha_init length must match ncol(W_outcome) + 1.");
    }
    
    // Add intercept term to Z_site (for γ)
    MatrixXd Z_site_int = prepend_intercept(Z_site);
    
    // Use only units with A_i == A_val for density ratio estimation
    std::vector<int> treated_idx = filter_treated_indices(A_source, A_val);
    
    if (treated_idx.empty()) {
        throw std::runtime_error("fit_unified_density_ratio_cpp: no observations with A_val in this fold.");
    }
    
    int n_treated = treated_idx.size();
    MatrixXd Z_site_treated = subset_rows(Z_site_int, treated_idx);
    
    // Initialize gamma: warm-start from previous solution if provided
    VectorXd gamma = (warm_start.size() == p_site) ? warm_start : VectorXd::Zero(p_site);
    
    // Store outcome covariates for treated units (for ψ' computation)
    MatrixXd W_outcome_treated(n_treated, p_outcome - 1);
    for (int i = 0; i < n_treated; i++) {
        W_outcome_treated.row(i) = W_outcome.row(treated_idx[i]);
    }
    
    // Pre-compute ψ'(η_i) for all treated units — this is CONSTANT across all
    // coordinate descent iterations since it depends only on alpha_init (not gamma).
    // Eliminates O(n_treated * p_outcome) VectorXd allocation + dot product per (j,i) pair.
    MatrixXd W_outcome_int = prepend_intercept(W_outcome_treated);
    VectorXd eta_outcome = W_outcome_int * alpha_init;
    
    // Convert GLM family/link integers to enums for response_derivative
    LinkFunction link = static_cast<LinkFunction>(link_int);
    
    VectorXd psi_prime_precomputed(n_treated);
    // Pre-compute ψ'(η_i) = |h'(η_i)| = b''(θ_i) > 0 for all link functions.
    // The absolute value is needed because the paper's ψ'(θ) = b''(θ) is always
    // positive, but response_derivative returns h'(η) which can be negative for
    // non-canonical parameterizations (e.g., inverse link: h'(η) = -1/η²).
    for (int i = 0; i < n_treated; i++) {
        double eta_i = eta_outcome(i);
        if (calibrated) {
            eta_i = truncation_function(eta_i, M_tau);
        }
        psi_prime_precomputed(i) = std::abs(GLMUtils::response_derivative(eta_i, link));
    }
    
    // Main iteration loop for density ratio estimation:
    // untruncated refined form (calibrated=false) or eq:gamma_calibrated_loss (calibrated=true).
    for (int iter = 0; iter < max_iter; iter++) {
        VectorXd gamma_old = gamma;
        
        // Compute g = Z_site_treated * gamma once per outer iteration (O(n_treated * p)),
        // then update incrementally after each coordinate change (O(n_treated) per coordinate).
        // This reduces total complexity from O(n_treated * p²) to O(n_treated * p) per iteration.
        VectorXd g_vec = Z_site_treated * gamma;
        
        // Coordinate descent: update γ_j one by one
        for (int j = 0; j < p_site; j++) {
            // First term: mean_grad_psi[j] represents
            // E_t[h'(W(X)^T alpha_init) * Z_tilde_j(X)], matching γ's basis.
            double grad_j = mean_grad_psi(j);
            
            // Second term: -Ẽ_{s_j}[I(A=1) φ_site(X)_j exp(-φ_site(X)^T γ) ψ'(η)]
            // where η = φ_outcome(X)^T α (using W_outcome features)
            // NOTE: ψ'(η_i) is pre-computed since it depends only on alpha_init (constant)
            StableAccumulator grad_acc;
            StableAccumulator hess_acc;
            for (int i = 0; i < n_treated; i++) {
                // exp(-φ_site(X)^T γ) uses Z_site features — read from incremental g_vec
                double current_g = std::max(NumericalConstants::ETA_CLIP_MIN, std::min(NumericalConstants::ETA_CLIP_MAX, g_vec(i)));
                double exp_neg_g = std::exp(-current_g);
                
                if (std::isnan(exp_neg_g) || std::isinf(exp_neg_g)) {
                    throw std::runtime_error("fit_unified_density_ratio_cpp: non-finite density-ratio weight encountered after logit clipping.");
                }
                exp_neg_g = std::max(NumericalConstants::WEIGHT_MIN, std::min(NumericalConstants::WEIGHT_MAX, exp_neg_g));
                
                // Use pre-computed ψ'(η_i) — avoids VectorXd allocation and dot product per (j,i)
                double psi_prime_i = psi_prime_precomputed(i);
                
                // φ_site(X)_j from Z_site features
                double phi_ij = Z_site_treated(i, j);
                grad_acc.add(-phi_ij * exp_neg_g * psi_prime_i / n);
                // Hessian diagonal: H_jj = Ẽ[φ_j² exp(-g) ψ'] / n
                hess_acc.add(phi_ij * phi_ij * exp_neg_g * psi_prime_i / n);
            }
            grad_j += grad_acc.value();
            double hess_j = hess_acc.value();
            
            // Hessian-based step size (matching CV function density_ratio_cd_update)
            hess_j = std::max(hess_j, NumericalConstants::HESSIAN_FLOOR);
            double step_size_hessian = 1.0 / hess_j;
            // Cap with adaptive step for safety (prevents huge steps early on)
            double step_size_adaptive = compute_adaptive_step_size(std::abs(grad_j), iter);
            double step_size = std::min(step_size_hessian, step_size_adaptive * 10.0);
            step_size = std::max(step_size, NumericalConstants::STEP_SIZE_FLOOR);
            
            // Proximal gradient update with soft-thresholding (no penalty on intercept)
            double old_gamma_j = gamma(j);
            gamma(j) = sanitize_param(proximal_update(gamma(j), grad_j, step_size, lambda, j));
            
            // Incremental g update: O(n_treated) instead of recomputing O(n_treated * p)
            double delta_j = gamma(j) - old_gamma_j;
            if (delta_j != 0.0) {
                g_vec += delta_j * Z_site_treated.col(j);
            }
        }
        
        // Check convergence (L∞ norm — consistent with outcome model and check_convergence_cpp)
        if (check_convergence_cpp(gamma_old, gamma, tol)) {
            return List::create(Named("gamma") = gamma,
                               Named("converged") = true,
                               Named("iterations") = iter + 1);
        }
    }
    
    return List::create(Named("gamma") = gamma,
                       Named("converged") = false,
                       Named("iterations") = max_iter);
}

// C++ version of INITIAL density ratio estimation (γ_init)
// Following eq:gamma_init in main.tex: INITIAL Loss for γ_{s_j,1}
// ℓ(γ_{s_j,1}) = (Ẽ_t[φ(X)])^T γ_{s_j,1} + Ẽ_{s_j}[I(A=1) exp(-φ(X)^T γ_{s_j,1})] + λ_γ ||γ_{s_j,1}||_1
// NOTE: This does NOT include ψ'(α_init) term - only uses mean_phi from target
//
// SAMPLE SIZE NORMALIZATION:
// The gradient contribution from source data is normalized by n (total sample size), NOT n_treated.
// This is correct because the empirical expectation Ẽ_{s_j}[I(A=1) f(X)] is defined as:
//   Ẽ_{s_j}[I(A=1) f(X)] = (1/n_s) * Σ_{i=1}^{n_s} I(A_i=1) * f(X_i)
// where the indicator I(A_i=1) is implicitly handled by summing only over treated units,
// but the normalization uses the full sample size n_s to match the definition.
// [[Rcpp::export]]
List fit_initial_density_ratio_cpp(const MatrixXd& Z_site, const VectorXd& A_source, 
                                   const VectorXd& mean_phi,
                                   double lambda, int max_iter, double tol,
                                   int A_val,
                                   const VectorXd& warm_start) {
    
    int n = Z_site.rows();
    MatrixXd Z_site_int = prepend_intercept(Z_site);
    int p = Z_site_int.cols();
    
    // Use only units with A_i == A_val
    auto treated_idx = filter_treated_indices(A_source, A_val);
    
    if (treated_idx.empty()) {
        return List::create(Named("gamma") = VectorXd::Zero(p),
                           Named("converged") = false,
                           Named("iterations") = 0);
    }
    
    int n_treated = treated_idx.size();
    MatrixXd X_treated = subset_rows(Z_site_int, treated_idx);
    
    // Initialize gamma: warm-start from previous solution if provided
    VectorXd gamma = (warm_start.size() == p) ? warm_start : VectorXd::Zero(p);
    
    // Main iteration loop for INITIAL density ratio estimation
    // Following eq:gamma_init in main.tex: NO ψ' term, only mean_phi and exp(-g)
    // FIXED: Use n (total sample size) for normalization, not n_treated
    // This matches the theory: Ẽ_{s_j}[I(A=1) f(X)] = (1/n_s) Σ I(A_i=1) f(X_i)
    for (int iter = 0; iter < max_iter; iter++) {
        VectorXd gamma_old = gamma;
        
        // Compute g = X_treated * gamma once per outer iteration (O(n_treated * p)),
        // then update incrementally after each coordinate change (O(n_treated) per coordinate).
        // This reduces total complexity from O(n_treated * p²) to O(n_treated * p) per iteration.
        VectorXd g_vec = X_treated * gamma;
        
        // Coordinate descent: update γ_j one by one
        for (int j = 0; j < p; j++) {
            // First term: Ẽ_t[φ(X)]_j = mean_phi[j]
            double grad_j = mean_phi(j);
            StableAccumulator grad_acc;
            StableAccumulator hess_acc;
            
            // Second term: -Ẽ_{s_j}[I(A=1) φ(X)_j exp(-φ(X)^T γ)]
            // NO ψ' term here - this is the key difference from refined loss
            // NOTE: Divide by n (total samples) not n_treated to match theory
            for (int i = 0; i < n_treated; i++) {
                double current_g = std::max(NumericalConstants::ETA_CLIP_MIN, std::min(NumericalConstants::ETA_CLIP_MAX, g_vec(i)));
                double exp_neg_g = std::exp(-current_g);
                
                if (std::isnan(exp_neg_g) || std::isinf(exp_neg_g)) {
                    exp_neg_g = 1.0;
                }
                exp_neg_g = std::max(NumericalConstants::WEIGHT_MIN, std::min(NumericalConstants::WEIGHT_MAX, exp_neg_g));
                
                double phi_ij = X_treated(i, j);
                grad_acc.add(-phi_ij * exp_neg_g / n);  // /n (total sample size) matches Ẽ_{s_j} convention
                // Hessian diagonal: H_jj = Ẽ[φ_j² exp(-g)] / n
                hess_acc.add(phi_ij * phi_ij * exp_neg_g / n);
            }
            grad_j += grad_acc.value();
            double hess_j = hess_acc.value();
            
            // Hessian-based step size (matching CV function density_ratio_cd_update)
            hess_j = std::max(hess_j, NumericalConstants::HESSIAN_FLOOR);
            double step_size_hessian = 1.0 / hess_j;
            // Cap with adaptive step for safety (prevents huge steps early on)
            double step_size_adaptive = compute_adaptive_step_size(std::abs(grad_j), iter);
            double step_size = std::min(step_size_hessian, step_size_adaptive * 10.0);
            step_size = std::max(step_size, NumericalConstants::STEP_SIZE_FLOOR);
            
            // Proximal gradient update with soft-thresholding (no penalty on intercept)
            double old_gamma_j = gamma(j);
            gamma(j) = sanitize_param(proximal_update(gamma(j), grad_j, step_size, lambda, j));
            
            // Incremental g update: O(n_treated) instead of recomputing O(n_treated * p)
            double delta_j = gamma(j) - old_gamma_j;
            if (delta_j != 0.0) {
                g_vec += delta_j * X_treated.col(j);
            }
        }
        
        // Check convergence (L∞ norm — consistent with outcome model and check_convergence_cpp)
        if (check_convergence_cpp(gamma_old, gamma, tol)) {
            return List::create(Named("gamma") = gamma,
                               Named("converged") = true,
                               Named("iterations") = iter + 1);
        }
    }
    
    return List::create(Named("gamma") = gamma,
                       Named("converged") = false,
                       Named("iterations") = max_iter);
}

// ============================================================================
// GLMNET-STYLE CV: Density Ratio Model
// Optimizations: Warm Start + Pathwise Descent + Pre-computed Folds + Active Set
// ============================================================================
// [[Rcpp::export]]
List select_lambda_cv_density_ratio_cpp(const MatrixXd& Z_site, const VectorXd& A_source, 
                                       const VectorXd& mean_grad_psi, const VectorXd& alpha_init,
                                       const VectorXd& lambda_grid, int n_folds,
                                       int max_iter, double tol,
                                       int A_val,
                                       int family_int, int link_int) {
    
    int n = Z_site.rows();
    int n_lambda = lambda_grid.size();
    int p = Z_site.cols() + 1;

    auto lambda_order = CVUtils::sort_lambda_descending(lambda_grid);
    auto treated_idx = CVUtils::filter_treated(A_source, A_val);
    if (treated_idx.empty()) {
        throw std::runtime_error("select_lambda_cv_density_ratio_cpp: no observations with A_val.");
    }
    int n_treated = treated_idx.size();
    if (n_treated < n_folds) {
        throw std::runtime_error("select_lambda_cv_density_ratio_cpp: not enough A_val observations for requested CV folds.");
    }

    // Pre-compute treated data with intercept
    MatrixXd X_treated_int = prepend_intercept(subset_rows(Z_site, treated_idx));

    // Pre-compute |ψ'(φ^T α_init)| = b''(θ) > 0 for treated units
    LinkFunction link = static_cast<LinkFunction>(link_int);
    VectorXd psi_prime_all(n_treated);
    for (int i = 0; i < n_treated; i++) {
        double eta_alpha = X_treated_int.row(i).dot(alpha_init);
        psi_prime_all(i) = std::abs(GLMUtils::response_derivative(eta_alpha, link));
    }

    auto folds = CVUtils::create_fold_splits(n_treated, n_folds);

    // Pre-allocate fold data
    std::vector<MatrixXd> X_train_folds(n_folds), X_val_folds(n_folds);
    std::vector<VectorXd> pp_train(n_folds), pp_val(n_folds);
    for (int fold = 0; fold < n_folds; fold++) {
        X_train_folds[fold] = CVUtils::slice_rows(X_treated_int, folds.train[fold]);
        pp_train[fold] = CVUtils::slice_elements(psi_prime_all, folds.train[fold]);
        X_val_folds[fold] = CVUtils::slice_rows(X_treated_int, folds.val[fold]);
        pp_val[fold] = CVUtils::slice_elements(psi_prime_all, folds.val[fold]);
    }

    // CV loop: outer=fold, inner=lambda (warm start across lambda path)
    MatrixXd fold_scores(n_lambda, n_folds);
    fold_scores.fill(INFINITY);

    for (int fold = 0; fold < n_folds; fold++) {
        if (folds.train[fold].empty() || folds.val[fold].empty()) continue;

        VectorXd gamma = VectorXd::Zero(p);
        std::vector<bool> active(p, true);

        for (int li = 0; li < n_lambda; li++) {
            int lambda_idx = lambda_order[li];
            double lambda = lambda_grid(lambda_idx);
            double cv_tol = std::max(tol, NumericalConstants::CV_TOL_FLOOR);
            int cv_max_iter = std::min(max_iter, NumericalConstants::CV_MAX_ITER);

            CVUtils::density_ratio_cd_update(gamma, active, X_train_folds[fold],
                                             pp_train[fold], n, mean_grad_psi,
                                             lambda, cv_tol, cv_max_iter);

            fold_scores(lambda_idx, fold) = CVUtils::density_ratio_val_loss(
                gamma, mean_grad_psi, X_val_folds[fold], pp_val[fold],
                static_cast<double>(n_treated) / n);
        }
    }

    return CVUtils::aggregate_cv_results(fold_scores, lambda_grid, n_lambda, n_folds);
}

// ============================================================================
// GLMNET-STYLE CV: INITIAL Density Ratio Model (no outcome model term)
// The initial density ratio loss (eq:gamma_init in main.tex) is a special case of the
// refined density ratio loss with ψ'(·) ≡ 1 (no outcome model weighting).
// ℓ(γ) = (Ẽ_t[φ(X)])^T γ + Ẽ_{s_j}[I(A=1) exp(-φ(X)^T γ)] + λ||γ||_1
// ============================================================================
// [[Rcpp::export]]
List select_lambda_cv_initial_density_ratio_cpp(const MatrixXd& Z_site, const VectorXd& A_source,
                                                const VectorXd& mean_phi,
                                                const VectorXd& lambda_grid, int n_folds,
                                                int max_iter, double tol,
                                                int A_val = 1) {
    
    int n = Z_site.rows();
    int n_lambda = lambda_grid.size();
    int p = Z_site.cols() + 1;

    auto lambda_order = CVUtils::sort_lambda_descending(lambda_grid);
    auto treated_idx = CVUtils::filter_treated(A_source, A_val);
    if (treated_idx.empty()) {
        throw std::runtime_error("select_lambda_cv_initial_density_ratio_cpp: no observations with A_val.");
    }
    int n_treated = treated_idx.size();
    if (n_treated < n_folds) {
        throw std::runtime_error("select_lambda_cv_initial_density_ratio_cpp: not enough A_val observations for requested CV folds.");
    }

    // Pre-compute treated data with intercept
    MatrixXd X_treated_int = prepend_intercept(subset_rows(Z_site, treated_idx));

    // For initial density ratio, ψ'(·) ≡ 1 (no outcome model dependency)
    VectorXd psi_prime_all = VectorXd::Ones(n_treated);

    auto folds = CVUtils::create_fold_splits(n_treated, n_folds);

    // Pre-allocate fold data
    std::vector<MatrixXd> X_train_folds(n_folds), X_val_folds(n_folds);
    std::vector<VectorXd> pp_train(n_folds), pp_val(n_folds);
    for (int fold = 0; fold < n_folds; fold++) {
        X_train_folds[fold] = CVUtils::slice_rows(X_treated_int, folds.train[fold]);
        pp_train[fold] = CVUtils::slice_elements(psi_prime_all, folds.train[fold]);
        X_val_folds[fold] = CVUtils::slice_rows(X_treated_int, folds.val[fold]);
        pp_val[fold] = CVUtils::slice_elements(psi_prime_all, folds.val[fold]);
    }

    // CV loop: outer=fold, inner=lambda (warm start across lambda path)
    MatrixXd fold_scores(n_lambda, n_folds);
    fold_scores.fill(INFINITY);

    for (int fold = 0; fold < n_folds; fold++) {
        if (folds.train[fold].empty() || folds.val[fold].empty()) continue;

        VectorXd gamma = VectorXd::Zero(p);
        std::vector<bool> active(p, true);

        for (int li = 0; li < n_lambda; li++) {
            int lambda_idx = lambda_order[li];
            double lambda = lambda_grid(lambda_idx);
            double cv_tol = std::max(tol, NumericalConstants::CV_TOL_FLOOR);
            int cv_max_iter = std::min(max_iter, NumericalConstants::CV_MAX_ITER);

            CVUtils::density_ratio_cd_update(gamma, active, X_train_folds[fold],
                                             pp_train[fold], n, mean_phi,
                                             lambda, cv_tol, cv_max_iter);

            fold_scores(lambda_idx, fold) = CVUtils::density_ratio_val_loss(
                gamma, mean_phi, X_val_folds[fold], pp_val[fold],
                static_cast<double>(n_treated) / n);
        }
    }

    return CVUtils::aggregate_cv_results(fold_scores, lambda_grid, n_lambda, n_folds);
}

// ============================================================================
// GLMNET-STYLE CV: Calibrated Density Ratio Model  
// Optimizations: Warm Start + Pathwise Descent + Pre-computed Folds + Active Set
// ============================================================================
// [[Rcpp::export]]
List select_lambda_cv_calibrated_density_ratio_cpp(const MatrixXd& Z_site, const VectorXd& A_source,
                                                   const VectorXd& mean_grad_psi, const VectorXd& alpha_init,
                                                   const VectorXd& lambda_grid, int n_folds, 
                                                   int max_iter, double tol, double M_tau,
                                                   const MatrixXd& W_outcome,
                                                   int A_val,
                                                   int family_int, int link_int) {
    
    int n = Z_site.rows();
    int n_lambda = lambda_grid.size();
    int p_site = Z_site.cols() + 1;

    if (W_outcome.rows() != n || W_outcome.cols() <= 0) {
        throw std::runtime_error("select_lambda_cv_calibrated_density_ratio_cpp: W_outcome must be non-empty and match Z_site rows.");
    }

    auto lambda_order = CVUtils::sort_lambda_descending(lambda_grid);
    auto treated_idx = CVUtils::filter_treated(A_source, A_val);
    if (treated_idx.empty()) {
        throw std::runtime_error("select_lambda_cv_calibrated_density_ratio_cpp: no observations with A_val.");
    }
    int n_treated = treated_idx.size();
    if (n_treated < n_folds) {
        throw std::runtime_error("select_lambda_cv_calibrated_density_ratio_cpp: not enough A_val observations for requested CV folds.");
    }

    // Prepare site features with intercept
    MatrixXd Z_site_int = prepend_intercept(subset_rows(Z_site, treated_idx));

    // Prepare outcome features and compute |ψ'| = |h'(η)| = b''(θ) > 0 with truncation
    LinkFunction link = static_cast<LinkFunction>(link_int);
    MatrixXd W_out_int = prepend_intercept(subset_rows(W_outcome, treated_idx));
    VectorXd eta_outcome = W_out_int * alpha_init;
    VectorXd psi_prime_all(n_treated);
    for (int i = 0; i < n_treated; i++) {
        double eta_truncated = truncation_function(eta_outcome(i), M_tau);
        psi_prime_all(i) = std::abs(GLMUtils::response_derivative(eta_truncated, link));
    }

    auto folds = CVUtils::create_fold_splits(n_treated, n_folds);

    // Pre-allocate fold data
    std::vector<MatrixXd> X_train_folds(n_folds), X_val_folds(n_folds);
    std::vector<VectorXd> pp_train(n_folds), pp_val(n_folds);
    for (int fold = 0; fold < n_folds; fold++) {
        X_train_folds[fold] = CVUtils::slice_rows(Z_site_int, folds.train[fold]);
        pp_train[fold] = CVUtils::slice_elements(psi_prime_all, folds.train[fold]);
        X_val_folds[fold] = CVUtils::slice_rows(Z_site_int, folds.val[fold]);
        pp_val[fold] = CVUtils::slice_elements(psi_prime_all, folds.val[fold]);
    }

    // CV loop
    MatrixXd fold_scores(n_lambda, n_folds);
    fold_scores.fill(INFINITY);

    for (int fold = 0; fold < n_folds; fold++) {
        if (folds.train[fold].empty() || folds.val[fold].empty()) continue;

        VectorXd gamma = VectorXd::Zero(p_site);
        std::vector<bool> active(p_site, true);

        for (int li = 0; li < n_lambda; li++) {
            int lambda_idx = lambda_order[li];
            double lambda = lambda_grid(lambda_idx);
            double cv_tol = std::max(tol, NumericalConstants::CV_TOL_FLOOR);
            int cv_max_iter = std::min(max_iter, NumericalConstants::CV_MAX_ITER);

            CVUtils::density_ratio_cd_update(gamma, active, X_train_folds[fold],
                                             pp_train[fold], n, mean_grad_psi,
                                             lambda, cv_tol, cv_max_iter);

            fold_scores(lambda_idx, fold) = CVUtils::density_ratio_val_loss(
                gamma, mean_grad_psi, X_val_folds[fold], pp_val[fold],
                static_cast<double>(n_treated) / n);
        }
    }

    return CVUtils::aggregate_cv_results(fold_scores, lambda_grid, n_lambda, n_folds);
}
