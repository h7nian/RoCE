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

#include "optimization.h"
#include "cv_utils.h"

namespace {

// Keep the final-iteration stopping diagnostics alongside every density-ratio
// fit.  A bare converged/iterations flag cannot distinguish a fit that missed
// the tolerance by a few percent from one that is genuinely stalled.  The R
// layer propagates these scalar diagnostics to the per-setting simulation QC.
List density_ratio_fit_result(const VectorXd& gamma, bool converged,
                              int iterations, double max_update,
                              double convergence_threshold,
                              int line_search_failures = 0) {
    return List::create(
        Named("gamma") = gamma,
        Named("converged") = converged,
        Named("iterations") = iterations,
        Named("max_update") = max_update,
        Named("convergence_threshold") = convergence_threshold,
        Named("line_search_failures") = line_search_failures,
        Named("max_abs_coefficient") = gamma.lpNorm<Eigen::Infinity>()
    );
}

}  // namespace

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
    
    // The refined (two-round) loss keeps its untruncated tilt: its CV path is
    // also run with an infinite radius from R.
    double tilt_radius = calibrated ? M_tau : std::numeric_limits<double>::infinity();
    std::vector<bool> active(p_site, true);
    CVUtils::DensityRatioCDResult fit = CVUtils::density_ratio_cd_update(
        gamma, active, Z_site_treated, psi_prime_precomputed,
        static_cast<double>(n_treated) / n, mean_grad_psi, lambda,
        tilt_radius, tol, max_iter, false
    );
    return density_ratio_fit_result(
        gamma, fit.converged, fit.iterations, fit.max_update,
        fit.convergence_threshold, fit.line_search_failures
    );
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
                                   int A_val, double M_tau,
                                   const VectorXd& warm_start) {

    int n = Z_site.rows();
    MatrixXd Z_site_int = prepend_intercept(Z_site);
    int p = Z_site_int.cols();
    
    // Use only units with A_i == A_val
    auto treated_idx = filter_treated_indices(A_source, A_val);
    
    if (treated_idx.empty()) {
        return density_ratio_fit_result(
            VectorXd::Zero(p), false, 0, NA_REAL, NA_REAL
        );
    }
    
    int n_treated = treated_idx.size();
    MatrixXd X_treated = subset_rows(Z_site_int, treated_idx);
    
    // Initialize gamma: warm-start from previous solution if provided
    VectorXd gamma = (warm_start.size() == p) ? warm_start : VectorXd::Zero(p);
    
    VectorXd psi_prime = VectorXd::Ones(n_treated);
    std::vector<bool> active(p, true);
    CVUtils::DensityRatioCDResult fit = CVUtils::density_ratio_cd_update(
        gamma, active, X_treated, psi_prime,
        static_cast<double>(n_treated) / n, mean_phi, lambda,
        M_tau, tol, max_iter, false
    );
    return density_ratio_fit_result(
        gamma, fit.converged, fit.iterations, fit.max_update,
        fit.convergence_threshold, fit.line_search_failures
    );
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
                                       int family_int, int link_int,
                                       Rcpp::Nullable<Rcpp::NumericVector> cv_fold_id = R_NilValue) {
    
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
    auto folds = CVUtils::create_fold_splits(n_treated, n_folds, cv_fold_id);

    // Pre-compute treated data with intercept
    MatrixXd X_treated_int = prepend_intercept(subset_rows(Z_site, treated_idx));

    // Pre-compute |ψ'(φ^T α_init)| = b''(θ) > 0 for treated units
    LinkFunction link = static_cast<LinkFunction>(link_int);
    VectorXd psi_prime_all(n_treated);
    for (int i = 0; i < n_treated; i++) {
        double eta_alpha = X_treated_int.row(i).dot(alpha_init);
        psi_prime_all(i) = std::abs(GLMUtils::response_derivative(eta_alpha, link));
    }

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
    std::vector<int> path_tail_skipped(n_folds, 0);

    const int cv_threads = CVUtils::nuisance_cv_thread_count(n_folds);
    ROCE_PARALLELIZE_CV_FOLDS(cv_threads)
    for (int fold = 0; fold < n_folds; fold++) {
        if (folds.train[fold].empty() || folds.val[fold].empty()) continue;

        VectorXd gamma = VectorXd::Zero(p);
        std::vector<bool> active(p, true);
        int consecutive_failures = 0;

        for (int li = 0; li < n_lambda; li++) {
            int lambda_idx = lambda_order[li];
            double lambda = lambda_grid(lambda_idx);
            double cv_tol = std::max(tol, NumericalConstants::CV_TOL_FLOOR);
            int cv_max_iter = std::min(max_iter, NumericalConstants::CV_MAX_ITER);
            VectorXd gamma_before = gamma;

            CVUtils::DensityRatioCDResult fit =
                CVUtils::density_ratio_cd_update(
                    gamma, active, X_train_folds[fold], pp_train[fold],
                    static_cast<double>(n_treated) / n, mean_grad_psi,
                    lambda, std::numeric_limits<double>::infinity(),
                    cv_tol, cv_max_iter
                );
            if (fit.converged) {
                consecutive_failures = 0;
                fold_scores(lambda_idx, fold) = CVUtils::density_ratio_val_loss(
                    gamma, mean_grad_psi, X_val_folds[fold], pp_val[fold],
                    static_cast<double>(n_treated) / n,
                    std::numeric_limits<double>::infinity());
            } else {
                // Never warm-start the next candidate from an unconverged
                // solution.  The +Inf score records this candidate as invalid.
                gamma = gamma_before;
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
                                                int A_val, double M_tau,
                                                Rcpp::Nullable<Rcpp::NumericVector> cv_fold_id = R_NilValue) {

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
    auto folds = CVUtils::create_fold_splits(n_treated, n_folds, cv_fold_id);

    // Pre-compute treated data with intercept
    MatrixXd X_treated_int = prepend_intercept(subset_rows(Z_site, treated_idx));

    // For initial density ratio, ψ'(·) ≡ 1 (no outcome model dependency)
    VectorXd psi_prime_all = VectorXd::Ones(n_treated);

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
    std::vector<int> path_tail_skipped(n_folds, 0);

    const int cv_threads = CVUtils::nuisance_cv_thread_count(n_folds);
    ROCE_PARALLELIZE_CV_FOLDS(cv_threads)
    for (int fold = 0; fold < n_folds; fold++) {
        if (folds.train[fold].empty() || folds.val[fold].empty()) continue;

        VectorXd gamma = VectorXd::Zero(p);
        std::vector<bool> active(p, true);
        int consecutive_failures = 0;

        for (int li = 0; li < n_lambda; li++) {
            int lambda_idx = lambda_order[li];
            double lambda = lambda_grid(lambda_idx);
            double cv_tol = std::max(tol, NumericalConstants::CV_TOL_FLOOR);
            int cv_max_iter = std::min(max_iter, NumericalConstants::CV_MAX_ITER);
            VectorXd gamma_before = gamma;

            CVUtils::DensityRatioCDResult fit =
                CVUtils::density_ratio_cd_update(
                    gamma, active, X_train_folds[fold], pp_train[fold],
                    static_cast<double>(n_treated) / n, mean_phi,
                    lambda, M_tau, cv_tol, cv_max_iter
                );
            if (fit.converged) {
                consecutive_failures = 0;
                fold_scores(lambda_idx, fold) = CVUtils::density_ratio_val_loss(
                    gamma, mean_phi, X_val_folds[fold], pp_val[fold],
                    static_cast<double>(n_treated) / n, M_tau);
            } else {
                gamma = gamma_before;
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
                                                   int family_int, int link_int,
                                                   Rcpp::Nullable<Rcpp::NumericVector> cv_fold_id = R_NilValue) {
    
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
    auto folds = CVUtils::create_fold_splits(n_treated, n_folds, cv_fold_id);

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
    std::vector<int> path_tail_skipped(n_folds, 0);

    const int cv_threads = CVUtils::nuisance_cv_thread_count(n_folds);
    ROCE_PARALLELIZE_CV_FOLDS(cv_threads)
    for (int fold = 0; fold < n_folds; fold++) {
        if (folds.train[fold].empty() || folds.val[fold].empty()) continue;

        VectorXd gamma = VectorXd::Zero(p_site);
        std::vector<bool> active(p_site, true);
        int consecutive_failures = 0;

        for (int li = 0; li < n_lambda; li++) {
            int lambda_idx = lambda_order[li];
            double lambda = lambda_grid(lambda_idx);
            double cv_tol = std::max(tol, NumericalConstants::CV_TOL_FLOOR);
            int cv_max_iter = std::min(max_iter, NumericalConstants::CV_MAX_ITER);
            VectorXd gamma_before = gamma;

            CVUtils::DensityRatioCDResult fit =
                CVUtils::density_ratio_cd_update(
                    gamma, active, X_train_folds[fold], pp_train[fold],
                    static_cast<double>(n_treated) / n, mean_grad_psi,
                    lambda, M_tau, cv_tol, cv_max_iter
                );
            if (fit.converged) {
                consecutive_failures = 0;
                fold_scores(lambda_idx, fold) = CVUtils::density_ratio_val_loss(
                    gamma, mean_grad_psi, X_val_folds[fold], pp_val[fold],
                    static_cast<double>(n_treated) / n, M_tau);
            } else {
                gamma = gamma_before;
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
