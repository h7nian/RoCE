#ifndef UTILS_HPP
#define UTILS_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <limits>
#include "numerical_constants.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// GLM Family enumeration
enum class GLMFamily {
    GAUSSIAN,    // Identity link, squared error loss
    BINOMIAL     // Logit link, logistic regression
};

// GLM Link function enumeration  
enum class LinkFunction {
    IDENTITY,    // η = μ
    LOGIT        // η = log(μ/(1-μ))
};

// GLM utility functions
class GLMUtils {
public:
    // Response function h(η) - maps linear predictor to mean
    static double response_function(double eta, LinkFunction link);
    
    // Canonical potential function H(η) such that h(η) = ∇H(η)
    static double potential_function(double eta, LinkFunction link);
    
    // Variance function V(μ) for exponential family
    static double variance_function(double mu, GLMFamily family);
    
    // Derivative of response function h'(η)
    static double response_derivative(double eta, LinkFunction link);
    
    // IRLS working weight: h'(η)² / V(h(η)), always positive
    // For canonical links (identity, logit, log): equals h'(η)
    // For non-canonical links (inverse): ensures positive-definite Hessian
    static double irls_working_weight(double eta, LinkFunction link, GLMFamily family);
    
    // Get canonical link for GLM family
    static LinkFunction get_canonical_link(GLMFamily family);
    
    // Loss function for GLM families (negative log-likelihood)
    static double calculate_loss(double y, double mu, GLMFamily family);
    
    // GLM negative log-likelihood (without constant terms)
    static double negative_log_likelihood(double y, double eta, GLMFamily family, LinkFunction link);

    // NLL gradient residual: -(y-μ)h'(η)/V(μ), the per-observation gradient of
    // the NLL w.r.t. the linear predictor η.  For canonical links (identity,
    // logit, log) this simplifies to (μ-y).  For the inverse link (Gamma) it
    // equals (y-μ) because h'(η)=-1/η² and V(μ)=1/η², giving a sign flip.
    static double nll_gradient_residual(double y, double mu, LinkFunction link);
};

// Helper function implementations

// Logistic function
inline double logistic_cpp(double x) {
    // Numerical stability handling
    if (x > NumericalConstants::ETA_CLIP_MAX) return 1.0;
    if (x < NumericalConstants::ETA_CLIP_MIN) return 0.0;
    return 1.0 / (1.0 + std::exp(-x));
}

// Vectorized logistic function
inline VectorXd logistic_vec_cpp(const VectorXd& x) {
    return x.unaryExpr(&logistic_cpp);
}

// Soft threshold function
inline double soft_threshold_cpp(double x, double lambda) {
    if (x > lambda) return x - lambda;
    if (x < -lambda) return x + lambda;
    return 0.0;
}

// Stable compensated summation accumulator (Neumaier variant) using long double.
// Useful for long reduction loops (gradient/Hessian/CV score accumulation).
struct StableAccumulator {
    long double sum{0.0L};
    long double c{0.0L};

    inline void add(double x) {
        long double t = sum + static_cast<long double>(x);
        if (std::fabs(sum) >= std::fabs(static_cast<long double>(x))) {
            c += (sum - t) + static_cast<long double>(x);
        } else {
            c += (static_cast<long double>(x) - t) + sum;
        }
        sum = t;
    }

    inline double value() const {
        long double v = sum + c;
        if (v > static_cast<long double>(std::numeric_limits<double>::max())) {
            return std::numeric_limits<double>::max();
        }
        if (v < static_cast<long double>(-std::numeric_limits<double>::max())) {
            return -std::numeric_limits<double>::max();
        }
        return static_cast<double>(v);
    }
};

// Proximal gradient update with L1 penalty (soft-thresholding).
// Applies soft_threshold for non-intercept coordinates (j > 0);
// intercept (j == 0) is updated without penalty (standard glmnet/RCAL practice).
// Returns the updated parameter value.
inline double proximal_update(double current_val, double grad, double step_size,
                              double lambda, int j) {
    double update_val = current_val - step_size * grad;
    if (j == 0) {
        return update_val;  // No penalty on intercept
    }
    return soft_threshold_cpp(update_val, lambda * step_size);
}

// Sanitize a parameter value: replace NaN/Inf/out-of-range with 0.
// Prevents numerical instability from propagating through coordinate descent.
inline double sanitize_param(double val, double max_abs = NumericalConstants::PARAM_MAX) {
    if (std::isnan(val) || std::isinf(val) || std::abs(val) > max_abs) {
        return 0.0;
    }
    return val;
}

// Prepend an intercept column (column of ones) to a matrix.
// Converts (n x p) matrix to (n x (p+1)) with ones in column 0.
inline MatrixXd prepend_intercept(const MatrixXd& X) {
    int n = X.rows();
    int p = X.cols();
    MatrixXd X_int(n, p + 1);
    X_int.col(0) = VectorXd::Ones(n);
    X_int.rightCols(p) = X;
    return X_int;
}

// Filter indices where A(i) == A_val (treated or control unit selection).
// Returns a vector of row indices for subsetting.
inline std::vector<int> filter_treated_indices(const VectorXd& A, int A_val) {
    std::vector<int> idx;
    idx.reserve(A.size());
    for (int i = 0; i < A.size(); i++) {
        if (static_cast<int>(A(i)) == A_val) {
            idx.push_back(i);
        }
    }
    return idx;
}

// Subset rows of a matrix by index vector.
inline MatrixXd subset_rows(const MatrixXd& X, const std::vector<int>& idx) {
    MatrixXd result(idx.size(), X.cols());
    for (size_t i = 0; i < idx.size(); i++) {
        result.row(i) = X.row(idx[i]);
    }
    return result;
}

// Subset elements of a vector by index vector.
inline VectorXd subset_elements(const VectorXd& v, const std::vector<int>& idx) {
    VectorXd result(idx.size());
    for (size_t i = 0; i < idx.size(); i++) {
        result(i) = v(idx[i]);
    }
    return result;
}

// Check convergence using absolute + relative infinity-norm criterion:
//   ||Δ||_∞ <= atol + rtol * max(1, ||new||_∞)
// This is more scale-robust than a pure absolute threshold.
inline double convergence_threshold_cpp(const VectorXd& new_params,
                                        double tol) {
    const double scale_inf =
        std::max(1.0, new_params.lpNorm<Eigen::Infinity>());
    constexpr double relative_tolerance = 1e-8;
    return tol + relative_tolerance * scale_inf;
}

inline bool check_convergence_cpp(const VectorXd& old_params, const VectorXd& new_params, double tol) {
    // Absolute + relative tolerance criterion for scale-robust stopping.
    const double delta_inf = (new_params - old_params).lpNorm<Eigen::Infinity>();
    return delta_inf <= convergence_threshold_cpp(new_params, tol);
}

// Create cross-validation folds using R's RNG for reproducibility
// This function uses Rcpp::runif() which respects set.seed() in R
inline VectorXd create_folds_cpp(int n, int n_folds) {
    VectorXd folds(n);
    for (int i = 0; i < n; i++) {
        folds(i) = i % n_folds;
    }
    
    // Random shuffle using R's RNG (respects set.seed() in R)
    // Use Fisher-Yates shuffle with R's runif for reproducibility
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    
    // Use Rcpp::runif which uses R's RNG state
    Rcpp::NumericVector rand_vals = Rcpp::runif(n);
    
    // Fisher-Yates shuffle using R's random values
    for (int i = n - 1; i > 0; i--) {
        // Use the random value to determine swap index
        int j = static_cast<int>(rand_vals[i] * (i + 1));
        j = std::min(j, i); // Ensure j <= i
        std::swap(indices[i], indices[j]);
    }
    
    VectorXd shuffled_folds(n);
    for (int i = 0; i < n; i++) {
        shuffled_folds(i) = folds(indices[i]);
    }
    
    return shuffled_folds;
}

// Calculate aggregated variance with finite sample corrections and cross-site covariances
// Unified function that handles both with and without penalty term
// Defined in weight_optimization.cpp (not inline — exported to R via Rcpp)
double calculate_aggregated_variance_cpp(const VectorXd& eta, const VectorXd& V_t, 
                                   const VectorXd& V_s, const VectorXd& n_s, 
                                   double V_ot, double n_t, const VectorXd& C_ot,
                                   const MatrixXd& C_cross,
                                   double lambda = 0.0, double mu_ot = 0.0);

// C++ version of matrix safe clipping (internal C++ helper, not exported to R)
inline MatrixXd safe_clip_matrix_cpp(const MatrixXd& X, double min_val, double max_val) {
    MatrixXd result = X;
    
    for (int i = 0; i < result.rows(); i++) {
        for (int j = 0; j < result.cols(); j++) {
            result(i, j) = std::max(min_val, std::min(max_val, result(i, j)));
        }
    }
    
    return result;
}

// C++ version of vector safe clipping (internal C++ helper, not exported to R)
inline VectorXd safe_clip_vector_cpp(const VectorXd& x, double min_val, double max_val) {
    VectorXd result = x;
    
    for (int i = 0; i < result.size(); i++) {
        result(i) = std::max(min_val, std::min(max_val, result(i)));
    }
    
    return result;
}

// C++ version of safe variance calculation (internal C++ helper, not exported to R)
inline double safe_var_cpp(const VectorXd& x, double min_var) {
    if (x.size() <= 1) return min_var;
    
    double mean_x = x.mean();
    // Use (n-1) denominator (Bessel's correction) to match R's var()
    double var_x = (x.array() - mean_x).square().sum() / static_cast<double>(x.size() - 1);
    
    if (std::isnan(var_x) || var_x <= 0) return min_var;
    
    return std::max(var_x, min_var);
}

// C++ version of matrix multiplication with numerical stability (internal C++ helper, not exported to R)
inline MatrixXd stable_matrix_multiply_cpp(const MatrixXd& A, const MatrixXd& B, double reg_factor) {
    
    if (A.cols() != B.rows()) {
        throw std::runtime_error("Matrix dimensions incompatible for multiplication");
    }
    
    MatrixXd result = A * B;
    
    if (result.rows() == result.cols()) {
        result.diagonal() += reg_factor * VectorXd::Ones(result.rows());
    }
    
    return result;
}

// C++ version of matrix inversion with numerical stability (internal C++ helper, not exported to R)
inline MatrixXd stable_matrix_inverse_cpp(const MatrixXd& A, double reg_factor) {
    
    if (A.rows() != A.cols()) {
        throw std::runtime_error("Matrix must be square for inversion");
    }
    
    MatrixXd A_reg = A;
    A_reg.diagonal() += reg_factor * VectorXd::Ones(A.rows());
    
    Eigen::LDLT<MatrixXd> ldlt(A_reg);
    
    if (ldlt.info() != Eigen::Success) {
        throw std::runtime_error("Matrix decomposition failed");
    }
    
    return ldlt.solve(MatrixXd::Identity(A.rows(), A.cols()));
}

// C++ version of solving linear system with numerical stability (internal C++ helper, not exported to R)
inline VectorXd stable_solve_cpp(const MatrixXd& A, const VectorXd& b, double reg_factor) {
    
    if (A.rows() != A.cols()) {
        throw std::runtime_error("Matrix must be square for solving");
    }
    
    if (A.rows() != b.size()) {
        throw std::runtime_error("Matrix and vector dimensions incompatible");
    }
    
    MatrixXd A_reg = A;
    A_reg.diagonal() += reg_factor * VectorXd::Ones(A.rows());
    
    Eigen::LDLT<MatrixXd> ldlt(A_reg);
    
    if (ldlt.info() != Eigen::Success) {
        throw std::runtime_error("Matrix decomposition failed");
    }
    
    return ldlt.solve(b);
}

// GLMUtils class implementation
inline double GLMUtils::response_function(double eta, LinkFunction link) {
    switch (link) {
        case LinkFunction::IDENTITY:
            return eta;
        case LinkFunction::LOGIT:
            // h(η) = exp(η)/(1 + exp(η)) = 1/(1 + exp(-η))
            eta = std::max(NumericalConstants::ETA_CLIP_MIN, std::min(NumericalConstants::ETA_CLIP_MAX, eta)); // Numerical stability
            return 1.0 / (1.0 + std::exp(-eta));
        default:
            return eta;
    }
}

inline double GLMUtils::potential_function(double eta, LinkFunction link) {
    // H(η) such that h(η) = ∇H(η)
    switch (link) {
        case LinkFunction::IDENTITY:
            // H(η) = η²/2
            return 0.5 * eta * eta;
        case LinkFunction::LOGIT:
            // H(η) = log(1 + exp(η))
            eta = std::max(NumericalConstants::ETA_CLIP_MIN, std::min(NumericalConstants::ETA_CLIP_MAX, eta)); // Numerical stability
            if (eta > NumericalConstants::LOG_SUM_EXP_THRESHOLD) {
                return eta; // log(1 + exp(η)) ≈ η for large η
            } else {
                return std::log(1.0 + std::exp(eta));
            }
        default:
            return 0.5 * eta * eta;
    }
}

inline double GLMUtils::variance_function(double mu, GLMFamily family) {
    switch (family) {
        case GLMFamily::GAUSSIAN:
            return 1.0; // V(μ) = 1
        case GLMFamily::BINOMIAL:
            // V(μ) = μ(1-μ)
            mu = std::max(NumericalConstants::HESSIAN_FLOOR, std::min(1.0 - NumericalConstants::HESSIAN_FLOOR, mu)); // Ensure μ ∈ (0,1)
            return mu * (1.0 - mu);
        default:
            return 1.0;
    }
}

inline double GLMUtils::response_derivative(double eta, LinkFunction link) {
    switch (link) {
        case LinkFunction::IDENTITY:
            return 1.0;
        case LinkFunction::LOGIT:
            // h'(η) = h(η)(1 - h(η))
            {
                double h = response_function(eta, link);
                return h * (1.0 - h);
            }
        default:
            return 1.0;
    }
}

// IRLS working weight = h'(η)² / V(μ), always positive.
// For canonical links: h'(η) = V(μ), so working weight = h'(η). No change.
// For inverse link: h'(η) = -1/η², V(μ) = μ² = 1/η²,
//   so working weight = (1/η⁴)/(1/η²) = 1/η² > 0.
// This ensures the Hessian diagonal is always positive for all link functions.
inline double GLMUtils::irls_working_weight(double eta, LinkFunction link, GLMFamily family) {
    (void)family;
    return response_derivative(eta, link);
}

inline LinkFunction GLMUtils::get_canonical_link(GLMFamily family) {
    switch (family) {
        case GLMFamily::GAUSSIAN:
            return LinkFunction::IDENTITY;
        case GLMFamily::BINOMIAL:
            return LinkFunction::LOGIT;
        default:
            return LinkFunction::IDENTITY;
    }
}

inline double GLMUtils::negative_log_likelihood(double y, double eta, GLMFamily family, LinkFunction link) {
    double mu = response_function(eta, link);
    
    switch (family) {
        case GLMFamily::GAUSSIAN:
            // -log p(y|μ) = (y-μ)²/2 (ignoring constant terms)
            return 0.5 * (y - mu) * (y - mu);
        case GLMFamily::BINOMIAL:
            // -log p(y|μ) = -[y*log(μ) + (1-y)*log(1-μ)]
            mu = std::max(NumericalConstants::PROBABILITY_FLOOR,
                          std::min(1.0 - NumericalConstants::PROBABILITY_FLOOR, mu));
            return -(y * std::log(mu) + (1.0 - y) * std::log(1.0 - mu));
        default:
            return 0.5 * (y - mu) * (y - mu);
    }
}

inline double GLMUtils::nll_gradient_residual(double y, double mu, LinkFunction link) {
    (void)link;
    return (mu - y);
}

inline double GLMUtils::calculate_loss(double y, double mu, GLMFamily family) {
    // Compute loss directly on response scale (mu), avoiding double-application
    // of response_function that would occur if delegated to negative_log_likelihood.
    switch (family) {
        case GLMFamily::GAUSSIAN:
            return 0.5 * (y - mu) * (y - mu);
        case GLMFamily::BINOMIAL:
            mu = std::max(NumericalConstants::PROBABILITY_FLOOR,
                          std::min(1.0 - NumericalConstants::PROBABILITY_FLOOR, mu));
            return -(y * std::log(mu) + (1.0 - y) * std::log(1.0 - mu));
        default:
            return 0.5 * (y - mu) * (y - mu);
    }
}

#endif // UTILS_HPP
