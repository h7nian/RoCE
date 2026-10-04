// Diagnostic only: compile against a frozen RoCE source snapshot. This follows
// its density-CV path and records each attempted/skipped lambda without editing
// the installed estimator. The exhaustive option audits its failure heuristic.
// attempted marks a visited path candidate; certified_nonconvergence identifies
// visited candidates handled without calling the optimizer. Tail skips are not
// marked attempted. The original behavior is unchanged when the flag is false.
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp14)]]
#include "cv_utils.h"
#include <chrono>

// A sufficient certificate that no point can meet the requested L-infinity
// KKT tolerance for this convex, affinely extended exponential-tilting loss.
// This is diagnostic and does not alter the path or its stopping decisions.
double recession_margin(const MatrixXd& x, const VectorXd& h, double scale,
                        const VectorXd& linear, const VectorXd& gamma,
                        double lambda, double radius, double tolerance) {
    if (!std::isfinite(radius) || radius <= 0 ||
        radius > std::min(std::log(NumericalConstants::WEIGHT_MAX),
                         -std::log(NumericalConstants::WEIGHT_MIN)) ||
        !gamma.allFinite() || h.minCoeff() < 0) return NA_REAL;
    const double norm = gamma.cwiseAbs().maxCoeff();
    if (norm == 0) return NA_REAL;
    const VectorXd direction = gamma / norm;
    long double value = 0, absolute_sum = 0, l1 = 0;
    for (int j = 0; j < direction.size(); ++j) {
        const long double term = static_cast<long double>(linear(j)) * direction(j);
        value += term;
        absolute_sum += std::abs(term);
        l1 += std::abs(direction(j));
        if (j > 0) {
            const long double penalty = lambda * std::abs(direction(j));
            value += penalty;
            absolute_sum += penalty;
        }
    }
    for (int i = 0; i < x.rows(); ++i) {
        long double projection = 0;
        for (int j = 0; j < x.cols(); ++j)
            projection += static_cast<long double>(x(i, j)) * direction(j);
        const long double slope = projection >= 0 ? -std::exp(-radius) : -std::exp(radius);
        const long double term = scale * h(i) * slope * projection / x.rows();
        value += term;
        absolute_sum += std::abs(term);
    }
    // A negative returned margin certifies a recession slope below the KKT
    // tolerance, with an additional conservative floating-point margin.
    return static_cast<double>(value + tolerance * l1 +
        1e-10L * std::max(1.0L, absolute_sum));
}

// [[Rcpp::export]]
Rcpp::List profile_density_grid_cpp(
        const MatrixXd& Z, const VectorXd& A, const VectorXd& linear,
        const VectorXd& h_all, const VectorXd& lambdas,
        int arm, int n_folds, double radius, double tolerance, int max_iter,
        bool exhaustive = false, bool center_design = false,
        Rcpp::Nullable<Rcpp::NumericVector> cv_fold_id = R_NilValue,
        bool use_kkt_certificate = false) {
    auto selected = CVUtils::filter_treated(A, arm);
    auto folds = CVUtils::create_fold_splits(selected.size(), n_folds, cv_fold_id);
    MatrixXd design = prepend_intercept(subset_rows(Z, selected));
    VectorXd h = CVUtils::slice_elements(h_all, selected);
    const double scale = static_cast<double>(selected.size()) / Z.rows();
    const int count = lambdas.size(), p = design.cols();
    auto order = CVUtils::sort_lambda_descending(lambdas);
    MatrixXd scores = MatrixXd::Constant(count, n_folds, INFINITY);
    const int total = count * n_folds;
    Rcpp::IntegerVector fold_column(total), index_column(total), iterations(total);
    Rcpp::NumericVector penalty(total), seconds(total), kkt(total, NA_REAL),
        coefficient_max(total, NA_REAL), intercept(total, NA_REAL),
        slope_max(total, NA_REAL), certificate(total, NA_REAL);
    Rcpp::LogicalVector attempted(total, false), converged(total, false);
    Rcpp::LogicalVector certified_nonconvergence(total, false);
    Rcpp::NumericVector fold_ids(selected.size());
    const double cv_tol = std::max(tolerance, NumericalConstants::CV_TOL_FLOOR);
    const int cv_max_iter = std::min(max_iter, NumericalConstants::CV_MAX_ITER);
    for (int fold = 0; fold < n_folds; ++fold) {
        MatrixXd train = CVUtils::slice_rows(design, folds.train[fold]);
        MatrixXd validation = CVUtils::slice_rows(design, folds.val[fold]);
        VectorXd train_h = CVUtils::slice_elements(h, folds.train[fold]);
        VectorXd validation_h = CVUtils::slice_elements(h, folds.val[fold]);
        MatrixXd fit_design = train;
        VectorXd fit_linear = linear;
        VectorXd center = VectorXd::Zero(p);
        double fit_tolerance = cv_tol;
        if (center_design) {
            for (int j = 1; j < p; ++j) {
                center(j) = train.col(j).mean();
                fit_design.col(j).array() -= center(j);
                fit_linear(j) -= center(j) * linear(0);
            }
            fit_tolerance /= 1.0 + center.cwiseAbs().maxCoeff();
        }
        for (int i : folds.val[fold]) fold_ids[i] = fold + 1;
        VectorXd gamma = VectorXd::Zero(p);
        std::vector<bool> active(p, true);
        int failures = 0;
        int certified_failures = 0;
        double certified_bound = 0.0;
        bool skip = false;
        for (int position = 0; position < count; ++position) {
            const int index = order[position], row = fold * count + index;
            fold_column[row] = fold + 1;
            index_column[row] = index + 1;
            penalty[row] = lambdas[index];
            if (skip) continue;
            attempted[row] = true;
            VectorXd before = gamma;
            auto start = std::chrono::steady_clock::now();
            const int certified_before = certified_failures;
            auto fit = CVUtils::density_ratio_cv_fit(gamma, active, fit_design,
                train_h, scale, fit_linear, lambdas[index], radius, fit_tolerance, cv_max_iter,
                use_kkt_certificate, certified_bound, certified_failures);
            certified_nonconvergence[row] = certified_failures > certified_before;
            seconds[row] = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
            iterations[row] = fit.iterations;
            VectorXd original_gamma = gamma;
            original_gamma(0) -= center.dot(gamma);
            kkt[row] = CVUtils::density_ratio_score_residual(original_gamma, train,
                train_h, scale, linear, lambdas[index], radius);
            converged[row] = fit.converged && kkt[row] <= cv_tol;
            coefficient_max[row] = original_gamma.cwiseAbs().maxCoeff();
            intercept[row] = original_gamma(0);
            slope_max[row] = original_gamma.tail(p - 1).cwiseAbs().maxCoeff();
            if (converged[row]) {
                failures = 0;
                scores(index, fold) = CVUtils::density_ratio_val_loss(original_gamma, linear,
                    validation, validation_h, scale, radius);
            } else {
                certificate[row] = recession_margin(train, train_h, scale, linear,
                    original_gamma, lambdas[index], radius, cv_tol);
                gamma = before;
                std::fill(active.begin(), active.end(), true);
                int skipped = CVUtils::cv_path_tail_after_failure(failures, position, count);
                if (!exhaustive && skipped > 0) skip = true;
            }
        }
    }
    int skipped = 0;
    for (int i = 0; i < total; ++i) if (!attempted[i]) skipped++;
    Rcpp::List summary;
    std::string error;
    try { summary = CVUtils::aggregate_cv_results(scores, lambdas, count, n_folds, skipped); }
    catch (const std::exception& exception) { error = exception.what(); }
    return Rcpp::List::create(
        Rcpp::Named("summary") = summary, Rcpp::Named("error") = error,
        Rcpp::Named("cv_fold_id") = fold_ids, Rcpp::Named("fold_scores") = scores,
        Rcpp::Named("profile") = Rcpp::DataFrame::create(
            Rcpp::Named("fold") = fold_column, Rcpp::Named("lambda_index") = index_column,
            Rcpp::Named("lambda") = penalty, Rcpp::Named("attempted") = attempted,
            Rcpp::Named("converged") = converged, Rcpp::Named("seconds") = seconds,
            Rcpp::Named("certified_nonconvergence") = certified_nonconvergence,
            Rcpp::Named("iterations") = iterations, Rcpp::Named("kkt") = kkt,
            Rcpp::Named("max_coefficient") = coefficient_max,
            Rcpp::Named("intercept") = intercept, Rcpp::Named("max_slope") = slope_max,
            Rcpp::Named("recession_kkt_margin") = certificate));
}
