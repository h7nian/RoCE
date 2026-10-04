// Diagnostic only: reproduce the frozen calibrated-OR CV path and flush each
// candidate's timing/KKT result so long-running fits can be located precisely.
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp14)]]
#include "cv_utils.h"
#include <chrono>
#include <fstream>
#include <iomanip>

// [[Rcpp::export]]
Rcpp::List profile_outcome_grid_cpp(
        const MatrixXd& W, const VectorXd& Y, const VectorXd& A,
        const VectorXd& gamma, const MatrixXd& Z, const VectorXd& lambdas,
        int arm, int n_folds, double radius, double tolerance, int max_iter,
        int family_int, int link_int, bool use_weight_derivative,
        const std::string& progress_path,
        Rcpp::Nullable<Rcpp::NumericVector> cv_fold_id = R_NilValue) {
    auto selected = CVUtils::filter_treated(A, arm);
    auto folds = CVUtils::create_fold_splits(selected.size(), n_folds, cv_fold_id);
    MatrixXd design = prepend_intercept(subset_rows(W, selected));
    VectorXd outcome = subset_elements(Y, selected);
    MatrixXd weight_design = prepend_intercept(subset_rows(Z, selected));
    VectorXd weights = compute_density_ratio_weights(
        weight_design, gamma, true, radius, use_weight_derivative);
    const int count = lambdas.size(), p = design.cols();
    const int cv_max_iter = std::min(max_iter, NumericalConstants::CV_MAX_ITER);
    const double cv_tol = std::max(tolerance, NumericalConstants::CV_TOL_FLOOR);
    const double train_scale = static_cast<double>(selected.size()) / W.rows();
    const GLMFamily family = static_cast<GLMFamily>(family_int);
    const LinkFunction link = static_cast<LinkFunction>(link_int);
    auto order = CVUtils::sort_lambda_descending(lambdas);
    MatrixXd scores = MatrixXd::Constant(count, n_folds, INFINITY);
    Rcpp::NumericVector fold_ids(selected.size());
    std::ofstream progress(progress_path);
    if (!progress) Rcpp::stop("Cannot create the outcome-CV progress file.");
    progress << "fold,lambda_index,lambda,attempted,converged,seconds,iterations,kkt,max_coefficient\n";
    progress << std::setprecision(17);
    int skipped = 0;
    for (int fold = 0; fold < n_folds; ++fold) {
        MatrixXd train = CVUtils::slice_rows(design, folds.train[fold]);
        MatrixXd validation = CVUtils::slice_rows(design, folds.val[fold]);
        VectorXd train_y = CVUtils::slice_elements(outcome, folds.train[fold]);
        VectorXd validation_y = CVUtils::slice_elements(outcome, folds.val[fold]);
        VectorXd train_w = CVUtils::slice_elements(weights, folds.train[fold]) * train_scale;
        const double val_scale = static_cast<double>(folds.val[fold].size()) / W.rows();
        VectorXd validation_w = CVUtils::slice_elements(weights, folds.val[fold]) * val_scale;
        for (int i : folds.val[fold]) fold_ids[i] = fold + 1;
        VectorXd alpha = VectorXd::Zero(p);
        std::vector<bool> active(p, true);
        int failures = 0;
        bool skip = false;
        for (int position = 0; position < count; ++position) {
            const int index = order[position];
            progress << fold + 1 << ',' << index + 1 << ',' << lambdas[index] << ',';
            if (skip) {
                progress << "0,0,0,0,NA,NA\n";
                ++skipped;
                continue;
            }
            VectorXd before = alpha;
            auto start = std::chrono::steady_clock::now();
            auto fit = CVUtils::glm_penalized_fit(alpha, active, train, train_y,
                train_w, train.rows(), lambdas[index], link, family, cv_tol, cv_max_iter);
            const double seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
            progress << "1," << fit.converged << ',' << seconds << ',' << fit.iterations
                     << ',' << fit.kkt_residual << ',' << alpha.cwiseAbs().maxCoeff() << '\n';
            progress.flush();
            if (fit.converged) {
                failures = 0;
                scores(index, fold) = CVUtils::glm_val_loss(alpha, validation,
                    validation_y, validation_w, family, link);
            } else {
                alpha = before;
                std::fill(active.begin(), active.end(), true);
                skip = CVUtils::cv_path_tail_after_failure(failures, position, count) > 0;
            }
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("summary") = CVUtils::aggregate_cv_results(scores, lambdas, count, n_folds, skipped),
        Rcpp::Named("cv_fold_id") = fold_ids, Rcpp::Named("fold_scores") = scores,
        Rcpp::Named("arm_weights") = weights);
}
