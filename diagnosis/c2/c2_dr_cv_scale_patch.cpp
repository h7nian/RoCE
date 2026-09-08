// Diagnosis-only density-ratio CV kernels with validation source scaling.
//
// This file intentionally lives under diagnosis/ and is loaded with
// Rcpp::sourceCpp() only by C2 probe jobs. It does not modify RoCE package
// source or exported symbols.

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/cv_utils.h"

using namespace Rcpp;
using namespace Eigen;

namespace {

double diag_abs_response_derivative(double eta, int link_int) {
    if (link_int == 0) return 1.0;
    if (link_int == 1) {
        eta = std::max(NumericalConstants::ETA_CLIP_MIN,
                       std::min(NumericalConstants::ETA_CLIP_MAX, eta));
        double mu = 1.0 / (1.0 + std::exp(-eta));
        return mu * (1.0 - mu);
    }
    throw std::runtime_error("c2_dr_cv_scale_patch: unsupported link_int.");
}

double scaled_density_ratio_val_loss(const VectorXd& gamma,
                                     const VectorXd& mean_grad_psi,
                                     const MatrixXd& X_val,
                                     const VectorXd& psi_prime_val,
                                     double source_scale) {
    int n_val = X_val.rows();
    if (n_val <= 0) {
        throw std::runtime_error("scaled_density_ratio_val_loss: empty validation fold.");
    }
    if (!std::isfinite(source_scale) || source_scale <= 0.0) {
        throw std::runtime_error("scaled_density_ratio_val_loss: invalid source_scale.");
    }

    long double source_mean = 0.0L;
    long double c = 0.0L;
    for (int j = 0; j < n_val; j++) {
        double eta_gamma = X_val.row(j).dot(gamma);
        eta_gamma = std::max(NumericalConstants::ETA_CLIP_MIN,
                             std::min(NumericalConstants::ETA_CLIP_MAX, eta_gamma));
        double exp_neg_g = std::exp(-eta_gamma);
        exp_neg_g = std::max(NumericalConstants::WEIGHT_MIN,
                             std::min(NumericalConstants::WEIGHT_MAX, exp_neg_g));
        long double add = static_cast<long double>(exp_neg_g * psi_prime_val(j)) / n_val;
        long double t = source_mean + add;
        if (std::fabs(source_mean) >= std::fabs(add)) {
            c += (source_mean - t) + add;
        } else {
            c += (add - t) + source_mean;
        }
        source_mean = t;
    }

    return static_cast<double>(
        static_cast<long double>(mean_grad_psi.dot(gamma)) +
        static_cast<long double>(source_scale) * (source_mean + c)
    );
}

int cv_train_denominator(int n_source, int n_train_treated,
                         int n_treated, const std::string& train_scale) {
    if (train_scale == "production") return n_source;
    if (train_scale == "arm_fraction") {
        double arm_fraction = static_cast<double>(n_treated) / n_source;
        return std::max(1, static_cast<int>(std::round(n_train_treated / arm_fraction)));
    }
    throw std::runtime_error("c2_dr_cv_scale_patch: train_scale must be production or arm_fraction.");
}

List add_patch_metadata(List out, double source_scale, const std::string& train_scale) {
    out["cv_source_scale"] = source_scale;
    out["cv_train_scale"] = train_scale;
    out["cv_scale_patch"] = "validation_source_arm_fraction";
    return out;
}

}  // namespace

// [[Rcpp::export]]
List c2_select_lambda_cv_initial_density_ratio_scaled_cpp(
        const MatrixXd& Z_site, const VectorXd& A_source,
        const VectorXd& mean_phi,
        const VectorXd& lambda_grid, int n_folds,
        int max_iter, double tol,
        int A_val,
        std::string train_scale = "production") {

    int n = Z_site.rows();
    int n_lambda = lambda_grid.size();
    int p = Z_site.cols() + 1;

    auto lambda_order = CVUtils::sort_lambda_descending(lambda_grid);
    auto treated_idx = CVUtils::filter_treated(A_source, A_val);
    if (treated_idx.empty()) {
        throw std::runtime_error("c2_select_lambda_cv_initial_density_ratio_scaled_cpp: no observations with A_val.");
    }
    int n_treated = treated_idx.size();
    if (n_treated < n_folds) {
        throw std::runtime_error("c2_select_lambda_cv_initial_density_ratio_scaled_cpp: not enough A_val observations for requested CV folds.");
    }

    double source_scale = static_cast<double>(n_treated) / n;
    MatrixXd X_treated_int = prepend_intercept(subset_rows(Z_site, treated_idx));
    VectorXd psi_prime_all = VectorXd::Ones(n_treated);
    auto folds = CVUtils::create_fold_splits(n_treated, n_folds);

    std::vector<MatrixXd> X_train_folds(n_folds), X_val_folds(n_folds);
    std::vector<VectorXd> pp_train(n_folds), pp_val(n_folds);
    for (int fold = 0; fold < n_folds; fold++) {
        X_train_folds[fold] = CVUtils::slice_rows(X_treated_int, folds.train[fold]);
        pp_train[fold] = CVUtils::slice_elements(psi_prime_all, folds.train[fold]);
        X_val_folds[fold] = CVUtils::slice_rows(X_treated_int, folds.val[fold]);
        pp_val[fold] = CVUtils::slice_elements(psi_prime_all, folds.val[fold]);
    }

    MatrixXd fold_scores(n_lambda, n_folds);
    fold_scores.fill(INFINITY);

    for (int fold = 0; fold < n_folds; fold++) {
        if (folds.train[fold].empty() || folds.val[fold].empty()) continue;

        VectorXd gamma = VectorXd::Zero(p);
        std::vector<bool> active(p, true);
        int n_train_scale = cv_train_denominator(
            n, X_train_folds[fold].rows(), n_treated, train_scale
        );

        for (int li = 0; li < n_lambda; li++) {
            int lambda_idx = lambda_order[li];
            double lambda = lambda_grid(lambda_idx);
            double cv_tol = std::max(tol, NumericalConstants::CV_TOL_FLOOR);
            int cv_max_iter = std::min(max_iter, NumericalConstants::CV_MAX_ITER);

            CVUtils::density_ratio_cd_update(
                gamma, active, X_train_folds[fold], pp_train[fold],
                n_train_scale, mean_phi, lambda, cv_tol, cv_max_iter
            );

            fold_scores(lambda_idx, fold) = scaled_density_ratio_val_loss(
                gamma, mean_phi, X_val_folds[fold], pp_val[fold], source_scale
            );
        }
    }

    List out = CVUtils::aggregate_cv_results(fold_scores, lambda_grid, n_lambda, n_folds);
    return add_patch_metadata(out, source_scale, train_scale);
}

// [[Rcpp::export]]
List c2_select_lambda_cv_calibrated_density_ratio_scaled_cpp(
        const MatrixXd& Z_site, const VectorXd& A_source,
        const VectorXd& mean_grad_psi, const VectorXd& alpha_init,
        const VectorXd& lambda_grid, int n_folds,
        int max_iter, double tol, double M_tau,
        const MatrixXd& W_outcome,
        int A_val,
        int family_int, int link_int,
        std::string train_scale = "production") {

    int n = Z_site.rows();
    int n_lambda = lambda_grid.size();
    int p_site = Z_site.cols() + 1;

    if (W_outcome.rows() != n || W_outcome.cols() <= 0) {
        throw std::runtime_error("c2_select_lambda_cv_calibrated_density_ratio_scaled_cpp: W_outcome must be non-empty and match Z_site rows.");
    }
    if (family_int != 0 && family_int != 1) {
        throw std::runtime_error("c2_select_lambda_cv_calibrated_density_ratio_scaled_cpp: unsupported family_int.");
    }

    auto lambda_order = CVUtils::sort_lambda_descending(lambda_grid);
    auto treated_idx = CVUtils::filter_treated(A_source, A_val);
    if (treated_idx.empty()) {
        throw std::runtime_error("c2_select_lambda_cv_calibrated_density_ratio_scaled_cpp: no observations with A_val.");
    }
    int n_treated = treated_idx.size();
    if (n_treated < n_folds) {
        throw std::runtime_error("c2_select_lambda_cv_calibrated_density_ratio_scaled_cpp: not enough A_val observations for requested CV folds.");
    }

    double source_scale = static_cast<double>(n_treated) / n;
    MatrixXd Z_site_int = prepend_intercept(subset_rows(Z_site, treated_idx));
    MatrixXd W_out_int = prepend_intercept(subset_rows(W_outcome, treated_idx));
    VectorXd eta_outcome = W_out_int * alpha_init;
    VectorXd psi_prime_all(n_treated);
    for (int i = 0; i < n_treated; i++) {
        double eta_i = eta_outcome(i);
        if (std::isfinite(M_tau)) eta_i = truncation_function(eta_i, M_tau);
        psi_prime_all(i) = diag_abs_response_derivative(eta_i, link_int);
    }

    auto folds = CVUtils::create_fold_splits(n_treated, n_folds);
    std::vector<MatrixXd> X_train_folds(n_folds), X_val_folds(n_folds);
    std::vector<VectorXd> pp_train(n_folds), pp_val(n_folds);
    for (int fold = 0; fold < n_folds; fold++) {
        X_train_folds[fold] = CVUtils::slice_rows(Z_site_int, folds.train[fold]);
        pp_train[fold] = CVUtils::slice_elements(psi_prime_all, folds.train[fold]);
        X_val_folds[fold] = CVUtils::slice_rows(Z_site_int, folds.val[fold]);
        pp_val[fold] = CVUtils::slice_elements(psi_prime_all, folds.val[fold]);
    }

    MatrixXd fold_scores(n_lambda, n_folds);
    fold_scores.fill(INFINITY);

    for (int fold = 0; fold < n_folds; fold++) {
        if (folds.train[fold].empty() || folds.val[fold].empty()) continue;

        VectorXd gamma = VectorXd::Zero(p_site);
        std::vector<bool> active(p_site, true);
        int n_train_scale = cv_train_denominator(
            n, X_train_folds[fold].rows(), n_treated, train_scale
        );

        for (int li = 0; li < n_lambda; li++) {
            int lambda_idx = lambda_order[li];
            double lambda = lambda_grid(lambda_idx);
            double cv_tol = std::max(tol, NumericalConstants::CV_TOL_FLOOR);
            int cv_max_iter = std::min(max_iter, NumericalConstants::CV_MAX_ITER);

            CVUtils::density_ratio_cd_update(
                gamma, active, X_train_folds[fold], pp_train[fold],
                n_train_scale, mean_grad_psi, lambda, cv_tol, cv_max_iter
            );

            fold_scores(lambda_idx, fold) = scaled_density_ratio_val_loss(
                gamma, mean_grad_psi, X_val_folds[fold], pp_val[fold], source_scale
            );
        }
    }

    List out = CVUtils::aggregate_cv_results(fold_scores, lambda_grid, n_lambda, n_folds);
    return add_patch_metadata(out, source_scale, train_scale);
}
