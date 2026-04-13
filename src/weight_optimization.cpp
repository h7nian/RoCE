// ============================================================================
// weight_optimization.cpp — Weight optimization and aggregation functions
// ============================================================================
// Functions:
//   - optimize_weights_cpp:              Optimal weight vector η via proximal gradient
//   - calculate_aggregated_estimate_cpp: Aggregated point estimate μ̂
//   - calculate_aggregated_variance_cpp: Aggregated variance V̂
// ============================================================================

#include "optimization.hpp"
#include "cv_utils.hpp"
#include <Eigen/Eigenvalues>
#include <cstdlib>
#include <string>

namespace {

inline MatrixXd symmetrize_matrix(const MatrixXd& M) {
    return 0.5 * (M + M.transpose());
}

inline bool is_weight_debug_enabled() {
    const char* debug_env = std::getenv("FACEC_DEBUG_WEIGHTS");
    return (debug_env != nullptr) && (std::string(debug_env) == "1");
}

struct WeightSmoothContext {
    double V_ot;
    double n_t;
    const VectorXd& V_t;
    const VectorXd& V_s;
    const VectorXd& n_s;
    const VectorXd& C_ot;
    const MatrixXd& C_cross_sym;
    bool has_cross;
    double psd_ridge;
};

inline double compute_smooth_gradient_coordinate(
        int j, double eta_j, double sum_eta, double dotC, double Cv_j,
        const WeightSmoothContext& ctx) {
    double n_s_j = std::max(ctx.n_s(j), 1.0);

    double grad = -2.0 * (1.0 - sum_eta) * ctx.V_ot / ctx.n_t +
                  2.0 * eta_j * (ctx.V_t(j) / ctx.n_t + ctx.V_s(j) / n_s_j) +
                  2.0 * ctx.C_ot(j) * (1.0 - sum_eta) / ctx.n_t -
                  2.0 * dotC / ctx.n_t;

    if (ctx.has_cross) {
        grad += 2.0 * (Cv_j - ctx.C_cross_sym(j, j) * eta_j) / ctx.n_t;
    }

    if (ctx.psd_ridge > 0.0) {
        grad += ctx.psd_ridge * eta_j;
    }

    return grad;
}

inline double compute_diagonal_curvature_coordinate(
        int j, const WeightSmoothContext& ctx) {
    double n_s_j = std::max(ctx.n_s(j), 1.0);
    double hess_diag = 2.0 * (ctx.V_ot / ctx.n_t + ctx.V_t(j) / ctx.n_t + ctx.V_s(j) / n_s_j)
                     - 4.0 * ctx.C_ot(j) / ctx.n_t;

    if (ctx.psd_ridge > 0.0) {
        hess_diag += ctx.psd_ridge;
    }

    return std::max(hess_diag, NumericalConstants::VAR_MIN);
}

// Finite-sample safeguard: the plug-in quadratic form in eq:final_opt can be
// slightly indefinite due to estimation noise (and because we pass only the
// off-diagonal cross-site covariances via C_cross with diag set to 0).
//
// We add the minimal ridge r >= 0 such that the smooth Hessian + r I is
// positive semidefinite. This keeps the proximal coordinate updates well-posed
// while remaining asymptotically negligible when the plug-in Hessian converges.
inline double compute_psd_ridge_for_weight_hessian(
        const VectorXd& V_t, const VectorXd& V_s, const VectorXd& n_s,
        double V_ot, double n_t, const VectorXd& C_ot,
        const MatrixXd& C_cross_sym, bool has_cross) {

    int K = V_t.size();
    if (K <= 0) return 0.0;

    VectorXd ones = VectorXd::Ones(K);

    // Diagonal variance contribution: 2 * diag( V_t/n_t + V_s/n_s )
    VectorXd diagD(K);
    for (int j = 0; j < K; j++) {
        double n_s_j = std::max(n_s(j), 1.0);
        diagD(j) = V_t(j) / n_t + V_s(j) / n_s_j;
    }
    MatrixXd H = 2.0 * diagD.asDiagonal();

    // Rank-1 term from (1 - sum eta)^2 * V_ot/n_t
    H.noalias() += (2.0 * V_ot / n_t) * (ones * ones.transpose());

    // Target-only / source-assisted covariance term
    H.noalias() -= (2.0 / n_t) * (C_ot * ones.transpose() + ones * C_ot.transpose());

    // Cross-site covariance contribution (off-diagonal; diag is typically 0)
    if (has_cross) {
        H.noalias() += (2.0 / n_t) * C_cross_sym;
    }

    // Ensure symmetry before eigen-decomposition
    H = symmetrize_matrix(H);

    Eigen::SelfAdjointEigenSolver<MatrixXd> solver(H);
    if (solver.info() != Eigen::Success) {
        return 0.0;
    }

    double min_eig = solver.eigenvalues().minCoeff();
    constexpr double psd_eps = 1e-10;
    if (min_eig < psd_eps) {
        return psd_eps - min_eig;
    }
    return 0.0;
}

}  // namespace

// C++ version of optimize_weights function
// [[Rcpp::export]]
List optimize_weights_cpp(const VectorXd& estimates, const VectorXd& V_t, 
                         const VectorXd& V_s, const VectorXd& n_s, 
                         double V_ot, double n_t, const VectorXd& C_ot,
                         double lambda, double mu_ot, int max_iter, double tol,
                         const MatrixXd& C_cross,
                         const VectorXd& warm_start) {
    
    int K = estimates.size();
    // Warm-start from a previous solution if provided (consistent with
    // fit_unified_density_ratio_cpp and fit_unified_outcome_cpp).
    VectorXd eta = (warm_start.size() == K) ? warm_start : VectorXd::Zero(K);
    
    // Pre-calculate penalty terms
    VectorXd penalties = lambda * (mu_ot - estimates.array()).square();
    
    // Pre-check cross-site matrix validity once outside the iteration loop
    bool has_cross = (C_cross.rows() == K && C_cross.cols() == K);
    MatrixXd C_cross_sym = MatrixXd::Zero(0, 0);
    if (has_cross) {
        C_cross_sym = symmetrize_matrix(C_cross);
    }

    // PSD ridge safeguard for the smooth quadratic objective
    double psd_ridge = compute_psd_ridge_for_weight_hessian(
        V_t, V_s, n_s, V_ot, n_t, C_ot, C_cross_sym, has_cross);

    const WeightSmoothContext smooth_ctx{
        V_ot, n_t, V_t, V_s, n_s, C_ot, C_cross_sym, has_cross, psd_ridge
    };

    if (is_weight_debug_enabled() && psd_ridge > 0.0) {
        Rcpp::Rcout << "[optimize_weights_cpp] PSD ridge activated: " << psd_ridge << "\n";
    }

    for (int iter = 0; iter < max_iter; iter++) {
        VectorXd eta_old = eta;

        // ----------------------------------------------------------------
        // Incremental aggregates (updated O(1) per coordinate instead of
        // re-computing O(K) reductions inside every coordinate step):
        //   S    = Σ_j η_j         (used as "1−S" in gradient)
        //   DotC = η · C_ot        (used in gradient cross-term)
        //   Cv   = C_cross · η     (cross-site gradient contribution;
        //                           updated with O(K) column-AXPY per step)
        // ----------------------------------------------------------------
        double S    = eta.sum();
        double DotC = eta.dot(C_ot);
        // Initialise Cv as zero; fill from C_cross * eta only when available.
        // Two-statement form avoids Eigen expression-template type mismatch
        // in the ternary operator (lazy Product vs CwiseNullaryOp).
        VectorXd Cv = VectorXd::Zero(K);
        if (has_cross) Cv = C_cross_sym * eta;
        
        for (int j = 0; j < K; j++) {
            double grad_smooth = compute_smooth_gradient_coordinate(
                j, eta(j), S, DotC, Cv(j), smooth_ctx
            );

            // Proximal gradient step (soft-thresholding)
            // Diagonal Hessian: d²V/dη_j² includes covariance correction
            // Full derivation: 2*(V_ot/n_t + V_t(j)/n_t + V_s(j)/n_s_j) - 4*C_ot(j)/n_t
            // The -4*C_ot(j)/n_t comes from d²/dη_j² of the cross-term
            //   2*(1-S)*Σ_k η_k C_ot(k)/n_t  in the variance objective.
            // Cross-site C_cross contributes 0 to diagonal Hessian (only off-diag).
            double hess_diag = compute_diagonal_curvature_coordinate(j, smooth_ctx);
            double step_size = 1.0 / hess_diag;
            
            double old_eta_j  = eta(j);
            double update_val = old_eta_j - step_size * grad_smooth;
            eta(j) = soft_threshold_cpp(update_val, step_size * penalties(j));
            
            // Enhanced numerical stability for weight optimization
            if (std::isnan(eta(j)) || std::isinf(eta(j)) || std::abs(eta(j)) > NumericalConstants::WEIGHT_MAX_ABS) {
                eta(j) = 0.0;
            }

            // Incremental aggregate updates — O(1) for S/DotC, O(K) for Cv
            double delta_j = eta(j) - old_eta_j;
            if (delta_j != 0.0) {
                S    += delta_j;
                DotC += delta_j * C_ot(j);
                if (has_cross) {
                    Cv += delta_j * C_cross_sym.col(j);
                }
            }
        }
        
        if (check_convergence_cpp(eta_old, eta, tol)) {
            return List::create(Named("weights") = eta, 
                               Named("converged") = true, 
                               Named("iterations") = iter + 1);
        }
    }
    
    return List::create(Named("weights") = eta, 
                       Named("converged") = false, 
                       Named("iterations") = max_iter);
}

// C++ version of aggregated estimate calculation
// [[Rcpp::export]]
double calculate_aggregated_estimate_cpp(double mu_ot, const VectorXd& mu_ts, 
                                        const VectorXd& eta) {
    
    double target_weight = 1.0 - eta.sum();
    double aggregated = target_weight * mu_ot + eta.dot(mu_ts);
    
    return aggregated;
}

// Calculate aggregated variance with finite sample corrections and cross-site covariances
// Unified function that handles both with and without penalty term
// [[Rcpp::export]]
double calculate_aggregated_variance_cpp(const VectorXd& eta, const VectorXd& V_t, 
                                   const VectorXd& V_s, const VectorXd& n_s, 
                                   double V_ot, double n_t, const VectorXd& C_ot,
                                   const MatrixXd& C_cross,
                                   double lambda, double mu_ot) {
    
    int K = eta.size();
    StableAccumulator eta_sum_acc;
    for (int j = 0; j < K; j++) eta_sum_acc.add(eta(j));
    double sum_eta = eta_sum_acc.value();
    
    // FIXED: Don't early-return for sum_eta <= 0
    // When using unconstrained weights (clip_weights=FALSE), sum_eta can be 0 or negative
    // The variance formula is still mathematically valid in these cases:
    // - sum_eta = 0: estimate is target-only, variance = V_ot/n_t (recovered by formula)
    // - sum_eta < 0: rare but possible, formula still computes correct variance
    
    StableAccumulator var_acc;
    
    // Add target-only variance component with proper weighting
    // When sum_eta = 0, this becomes V_ot/n_t (target-only)
    var_acc.add(std::pow(1.0 - sum_eta, 2) * V_ot / n_t);
    
    // Add source variance and target-source covariance components  
    // Following eq:final_opt in main.tex: eta^2 * V_hat / N for each variance component
    for (int j = 0; j < K; j++) {
        // Source variance component (guard against division by zero)
        double n_s_j = std::max(n_s(j), 1.0);  // Prevent division by zero
        var_acc.add(std::pow(eta(j), 2) * (V_t(j) / n_t + V_s(j) / n_s_j));
        
        // Target-source covariance component
        var_acc.add(2.0 * eta(j) * (1.0 - sum_eta) * C_ot(j) / n_t);
    }
    
    // Add cross-site covariance terms: (2/N_t) * sum_{j<k} eta_j * eta_k * C_{t,s_j,s_k}
    if (C_cross.rows() == K && C_cross.cols() == K) {
        for (int j = 0; j < K; j++) {
            for (int k = j + 1; k < K; k++) {
                var_acc.add(2.0 * eta(j) * eta(k) * C_cross(j, k) / n_t);
            }
        }
    }
    
    // No penalty term in this version
    
    // Ensure non-negative variance with conservative minimum
    double variance = var_acc.value();
    return std::max(variance, NumericalConstants::VAR_MIN);
}

// ============================================================================
// Batch lambda grid evaluation for weight optimization
// ============================================================================
// Evaluates the aggregated variance for every lambda in a grid in a single
// C++ call, eliminating the per-lambda R→C++ round-trips that dominate
// select_lambda_cv_crossfit when the grid contains ~100 candidates.
//
// Pathwise warm-start: lambdas are traversed from LARGEST to SMALLEST
// (most-regularized → least-regularized).  The solution at a large lambda
// (η ≈ 0) is a good initialiser for the next smaller lambda, so the
// optimizer converges in far fewer iterations — same trick used by all four
// CV functions in cv_utils.hpp.
//
// Results are stored in original-grid order so the R caller's lambda → index
// mapping is preserved.
//
// Inputs are assumed to be already validated/clipped by the R caller
// (mirror the pre-processing inside optimize_weights / calculate_aggregated_variance).
// Returns a NumericVector of the same length as lambda_grid; entries are
// the aggregated variance for the corresponding lambda, or R_PosInf on failure.
// [[Rcpp::export]]
NumericVector evaluate_lambda_grid_weights_cpp(
        const VectorXd& estimates, const VectorXd& V_t, const VectorXd& V_s,
        const VectorXd& n_s_vec, double V_ot, double n_t, const VectorXd& C_ot,
        double mu_ot, const MatrixXd& C_cross, const VectorXd& lambda_grid,
        int max_iter, double tol) {

    int n_lambda = lambda_grid.size();
    NumericVector out(n_lambda, R_PosInf);

    if (n_lambda == 0 || estimates.size() == 0) return out;

    // Sort lambda indices DESCENDING (largest first) for pathwise warm-start.
    // Mirrors CVUtils::sort_lambda_descending used in all four C++ CV functions.
    std::vector<int> lambda_order(n_lambda);
    std::iota(lambda_order.begin(), lambda_order.end(), 0);
    std::sort(lambda_order.begin(), lambda_order.end(),
              [&lambda_grid](int a, int b) { return lambda_grid(a) > lambda_grid(b); });

    // eta carries the solution from the previous (larger) lambda as warm-start
    // for the current (smaller) lambda.  Starts at zero (correct for λ=∞).
    VectorXd eta_warm = VectorXd::Zero(estimates.size());

    for (int li = 0; li < n_lambda; li++) {
        int orig_idx = lambda_order[li];
        double lam   = lambda_grid(orig_idx);
        if (!std::isfinite(lam) || lam < 0.0) continue;

        // Weight optimization with pathwise warm-start (intra-C++ — no R boundary)
        List wt_res = optimize_weights_cpp(estimates, V_t, V_s, n_s_vec,
                                           V_ot, n_t, C_ot,
                                           lam, mu_ot, max_iter, tol,
                                           C_cross, eta_warm);

        bool converged = Rcpp::as<bool>(wt_res["converged"]);
        if (!converged) continue;

        // Extract weights; propagate as warm-start for next lambda
        VectorXd eta = wt_res["weights"];
        if (eta.size() == 0) continue;
        eta_warm = eta;   // warm-start for next (smaller) lambda

        // Variance computation (intra-C++ call)
        double var = calculate_aggregated_variance_cpp(eta, V_t, V_s, n_s_vec,
                                                       V_ot, n_t, C_ot, C_cross,
                                                       lam, mu_ot);
        out[orig_idx] = (std::isfinite(var) && var > 0.0) ? var : R_PosInf;
    }

    return out;
}
