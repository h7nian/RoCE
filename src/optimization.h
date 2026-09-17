#ifndef OPTIMIZATION_HPP
#define OPTIMIZATION_HPP

#include "utils.h"
#include "numerical_constants.h"

// ============================================================================
// Weight optimization
// ============================================================================

List optimize_weights_cpp(const VectorXd& estimates, const VectorXd& V_t,
                          const VectorXd& V_s, const VectorXd& n_s, 
                          double V_ot, double n_t, const VectorXd& C_ot,
                          double lambda, double mu_ot, int max_iter, double tol,
                          const MatrixXd& C_cross = MatrixXd::Zero(0,0),
                          const VectorXd& warm_start = VectorXd());

// ============================================================================
// Density ratio fitting (defined in density_ratio.cpp)
// ============================================================================

// Forward declarations only — defaults are specified in the .cpp definitions
List fit_unified_density_ratio_cpp(const MatrixXd& Z_site, const VectorXd& A_source,
                                   const VectorXd& mean_grad_psi, const VectorXd& alpha_init,
                                   double lambda, int max_iter, double tol,
                                   bool calibrated, double M_tau,
                                   const MatrixXd& W_outcome, int A_val,
                                   int family_int, int link_int,
                                   const VectorXd& warm_start);

List fit_initial_density_ratio_cpp(const MatrixXd& Z_site, const VectorXd& A_source,
                                   const VectorXd& mean_phi,
                                   double lambda, int max_iter, double tol,
                                   int A_val, double M_tau, const VectorXd& warm_start);

// ============================================================================
// Outcome model fitting (defined in outcome_model.cpp)
// ============================================================================

List fit_unified_outcome_cpp(const MatrixXd& W_outcome, const VectorXd& Y_source,
                             const VectorXd& A_source, const VectorXd& gamma_s,
                             int family_int, int link_int, double lambda,
                             int max_iter, double tol, int A_val,
                             bool calibrated, double M_tau,
                             const MatrixXd& Z_site, const VectorXd& warm_start);

// ============================================================================
// CV lambda selection
// ============================================================================

List select_lambda_cv_density_ratio_cpp(const MatrixXd& Z_site, const VectorXd& A_source, 
                                       const VectorXd& mean_grad_psi, const VectorXd& alpha_init,
                                       const VectorXd& lambda_grid, int n_folds,
                                       int max_iter, double tol,
                                       int A_val,
                                       int family_int, int link_int,
                                       Rcpp::Nullable<Rcpp::NumericVector> cv_fold_id);

List select_lambda_cv_initial_density_ratio_cpp(const MatrixXd& Z_site, const VectorXd& A_source,
                                                const VectorXd& mean_phi,
                                                const VectorXd& lambda_grid, int n_folds,
                                                int max_iter, double tol,
                                                int A_val, double M_tau,
                                                Rcpp::Nullable<Rcpp::NumericVector> cv_fold_id);

List select_lambda_cv_general_refined_outcome_cpp(const MatrixXd& W_outcome, const VectorXd& Y_source,
                                                 const VectorXd& A_source, const VectorXd& gamma_s,
                                                 int family_int, int link_int, const VectorXd& lambda_grid, 
                                                 int n_folds, int max_iter, double tol, int A_val,
                                                 const MatrixXd& Z_site,
                                                 Rcpp::Nullable<Rcpp::NumericVector> cv_fold_id);

List select_lambda_cv_calibrated_density_ratio_cpp(const MatrixXd& Z_site, const VectorXd& A_source,
                                                   const VectorXd& mean_grad_psi, const VectorXd& alpha_init,
                                                   const VectorXd& lambda_grid, int n_folds,
                                                   int max_iter, double tol, double M_tau,
                                                   const MatrixXd& W_outcome,
                                                   int A_val,
                                                   int family_int, int link_int,
                                                   Rcpp::Nullable<Rcpp::NumericVector> cv_fold_id);

List select_lambda_cv_calibrated_outcome_cpp(const MatrixXd& W_outcome, const VectorXd& Y_source,
                                             const VectorXd& A_source, const VectorXd& gamma_init,
                                             const VectorXd& lambda_grid, int n_folds,
                                             int max_iter, double tol, int A_val, double M_tau,
                                             const MatrixXd& Z_site,
                                             int family_int, int link_int,
                                             Rcpp::Nullable<Rcpp::NumericVector> cv_fold_id);

// ============================================================================
// GLM fitting
// ============================================================================

List fit_general_glm_cpp(const MatrixXd& X, const VectorXd& Y, const VectorXd& weights,
                        int family_int, int link_int, double lambda, 
                        int max_iter, double tol,
                        const VectorXd& warm_start = VectorXd());

// GLM utility functions (unified — default to BINOMIAL/LOGIT)
MatrixXd calculate_glm_gradient_cpp(const MatrixXd& W_outcome, const VectorXd& beta,
                                   int family_int = 1, int link_int = 1);
VectorXd predict_glm_cpp(const MatrixXd& W_outcome, const VectorXd& beta,
                        int family_int = 1, int link_int = 1);

// Column-means of GLM gradient — avoids n×p matrix round-trip through R.
// Returns colMeans( h'(Xβ) · X ) directly as a (p+1)-vector.
// Replaces: colMeans(calculate_glm_gradient_cpp(W, beta, fi, li))
VectorXd mean_glm_gradient_cpp(const MatrixXd& W_outcome, const VectorXd& beta,
                               int family_int = 1, int link_int = 1);

// ============================================================================
// Variance calculation functions
// ============================================================================

List calculate_correction_term_cpp(const MatrixXd& Z_site, const VectorXd& A_source,
                                  const VectorXd& Y_source, const VectorXd& gamma_s,
                                  const VectorXd& alpha_ts,
                                  const MatrixXd& W_outcome,
                                  double M_tau = NumericalConstants::M_TAU_DEFAULT,
                                  int family_int = 1, int link_int = 1,
                                  int A_val = 1);

List calculate_source_variance_cpp(const MatrixXd& Z_site, const VectorXd& A_source,
                                  const VectorXd& Y_source, const VectorXd& gamma_s,
                                  const VectorXd& alpha_ts, double delta_ts,
                                  const MatrixXd& W_outcome,
                                  double M_tau = NumericalConstants::M_TAU_DEFAULT,
                                  int family_int = 1, int link_int = 1,
                                  int A_val = 1);

List calculate_target_variance_cpp(const MatrixXd& target_W_outcome, const VectorXd& alpha_ts,
                                  double M_ts,
                                  int family_int = 1, int link_int = 1);

double calculate_covariance_term_cpp(const VectorXd& varphi_ot, const VectorXd& zeta_components);

double calculate_aggregated_estimate_cpp(double mu_ot, const VectorXd& mu_ts, 
                                        const VectorXd& eta);

// Aggregated variance functions declared in utils.h, defined in weight_optimization.cpp

#endif // OPTIMIZATION_HPP
