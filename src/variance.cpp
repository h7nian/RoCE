// ============================================================================
// variance.cpp — Variance, covariance, and correction term calculations
// ============================================================================
// Functions:
//   - calculate_correction_term_cpp:       DR correction δ̂_{t,s_j}
//   - calculate_source_variance_cpp:       Source variance V̂_s
//   - calculate_target_variance_cpp:       Target variance V̂_{t,s_j}
//   - calculate_covariance_term_cpp:       Target-source covariance Ĉ_{ot,s_j}
// ============================================================================

#include "optimization.hpp"
#include "cv_utils.hpp"
#include <limits>

// Forward declaration for predict_glm_cpp (defined in outcome_model.cpp)
VectorXd predict_glm_cpp(const MatrixXd& W_outcome, const VectorXd& beta,
                        int family_int, int link_int);

namespace {

struct ClipDiagnostics {
    int n_obs = 0;
    int weight_min_clipped = 0;
    int weight_max_clipped = 0;
    int ratio_min_clipped = 0;
    int ratio_max_clipped = 0;
    double max_abs_logit = 0.0;
    double max_raw_weight = 0.0;
};

VectorXd inference_density_ratios(const VectorXd& density_logits,
                                  double M_tau,
                                  ClipDiagnostics& diag) {
    int n = density_logits.size();
    diag.n_obs = n;
    VectorXd density_ratios(n);

    for (int i = 0; i < n; i++) {
        double raw_logit = density_logits(i);
        diag.max_abs_logit = std::max(diag.max_abs_logit, std::abs(raw_logit));

        double clipped_logit = truncation_function(raw_logit, M_tau);
        double raw_weight = std::exp(-clipped_logit);
        if (std::isfinite(raw_weight)) {
            diag.max_raw_weight = std::max(diag.max_raw_weight, raw_weight);
        } else {
            diag.max_raw_weight = std::numeric_limits<double>::infinity();
        }

        if (!std::isfinite(raw_weight) || raw_weight > NumericalConstants::WEIGHT_MAX) {
            diag.weight_max_clipped++;
            density_ratios(i) = NumericalConstants::WEIGHT_MAX;
        } else if (raw_weight < NumericalConstants::WEIGHT_MIN) {
            diag.weight_min_clipped++;
            density_ratios(i) = NumericalConstants::WEIGHT_MIN;
        } else {
            density_ratios(i) = raw_weight;
        }
    }

    return density_ratios;
}

double clipped_ratio_term(double ratio_term, ClipDiagnostics& diag) {
    if (ratio_term < NumericalConstants::RATIO_CLIP_MIN) {
        diag.ratio_min_clipped++;
        return NumericalConstants::RATIO_CLIP_MIN;
    }
    if (ratio_term > NumericalConstants::RATIO_CLIP_MAX) {
        diag.ratio_max_clipped++;
        return NumericalConstants::RATIO_CLIP_MAX;
    }
    return ratio_term;
}

List make_clip_diagnostics(const ClipDiagnostics& diag) {
    bool any_clipped = (diag.weight_min_clipped + diag.weight_max_clipped +
                        diag.ratio_min_clipped + diag.ratio_max_clipped) > 0;
    return List::create(
        Named("n_obs") = diag.n_obs,
        Named("weight_min_clipped") = diag.weight_min_clipped,
        Named("weight_max_clipped") = diag.weight_max_clipped,
        Named("ratio_min_clipped") = diag.ratio_min_clipped,
        Named("ratio_max_clipped") = diag.ratio_max_clipped,
        Named("max_abs_logit") = diag.max_abs_logit,
        Named("max_raw_weight") = diag.max_raw_weight,
        Named("any_clipped") = any_clipped
    );
}

} // namespace

// C++ version of correction term calculation
// Supports separate Z_site (for density ratio) and W_outcome (for predictions)
// M_tau ensures consistent truncation between correction and variance (Problem 6 fix)
// family_int/link_int control the GLM response function for outcome predictions
// [[Rcpp::export]]
List calculate_correction_term_cpp(const MatrixXd& Z_site, const VectorXd& A_source,
                                  const VectorXd& Y_source, const VectorXd& gamma_s,
                                  const VectorXd& alpha_ts,
                                  const MatrixXd& W_outcome,
                                  double M_tau,
                                  int family_int, int link_int,
                                  int A_val) {
    
    if (A_val != 0 && A_val != 1) {
        throw std::runtime_error("calculate_correction_term_cpp: A_val must be 0 or 1.");
    }

    int n = Z_site.rows();
    
    if (W_outcome.rows() != n || W_outcome.cols() <= 0) {
        throw std::runtime_error("calculate_correction_term_cpp: W_outcome must be non-empty and match Z_site rows.");
    }
    
    // Add intercept for density ratio calculation (using Z_site)
    MatrixXd Z_site_int = prepend_intercept(Z_site);
    
    // Calculate density ratios using Z_site and gamma_s
    VectorXd density_logits = Z_site_int * gamma_s;
    
    ClipDiagnostics clip_diag;
    VectorXd density_ratios = inference_density_ratios(density_logits, M_tau, clip_diag);
    
    // Calculate outcome predictions using W_outcome
    VectorXd outcome_preds = predict_glm_cpp(W_outcome, alpha_ts, family_int, link_int);
    
    // Calculate correction components with numerical stability
    // Following eq:if_source_component / eq:dr_estimator in main.tex:
    // w_a(X) = I(A=a) / exp(φ(X)^T γ).
    // Since density_ratios = exp(-g), the A_i == A_val branch multiplies by exp(-g).
    VectorXd correction_components(A_source.size());
    for (int i = 0; i < A_source.size(); i++) {
        if (A_source(i) == A_val) {
            // density_ratios = exp(-g), and this branch enforces A_i == A_val,
            // so ratio term is exactly 1/exp(g) for included units.
            double ratio_term = clipped_ratio_term(density_ratios(i), clip_diag);
            correction_components(i) = ratio_term * (Y_source(i) - outcome_preds(i));
        } else {
            correction_components(i) = 0.0;
        }
    }
    
    // Calculate mean: δ̂_{t,s_j}^1 = Ẽ_{s_j}[w(X)(Y - ψ(X; α̂))]
    double delta_ts = correction_components.mean();
    
    return List::create(Named("delta_ts") = delta_ts,
                       Named("correction_components") = correction_components,
                       Named("clip_diagnostics") = make_clip_diagnostics(clip_diag));
}

// C++ version of source variance calculation
// Supports separate Z_site (for density ratio) and W_outcome (for predictions)
// M_tau parameter ensures consistent truncation with training (Problem 6 fix)
// family_int/link_int control the GLM response function for outcome predictions
// [[Rcpp::export]]
List calculate_source_variance_cpp(const MatrixXd& Z_site, const VectorXd& A_source,
                                  const VectorXd& Y_source, const VectorXd& gamma_s,
                                  const VectorXd& alpha_ts, double delta_ts,
                                  const MatrixXd& W_outcome,
                                  double M_tau,
                                  int family_int, int link_int,
                                  int A_val) {
    
    if (A_val != 0 && A_val != 1) {
        throw std::runtime_error("calculate_source_variance_cpp: A_val must be 0 or 1.");
    }

    int n = Z_site.rows();
    
    if (W_outcome.rows() != n || W_outcome.cols() <= 0) {
        throw std::runtime_error("calculate_source_variance_cpp: W_outcome must be non-empty and match Z_site rows.");
    }
    
    // Add intercept for density ratio calculation (using Z_site)
    MatrixXd Z_site_int = prepend_intercept(Z_site);
    
    // Calculate density ratios using Z_site and gamma_s
    VectorXd density_logits = Z_site_int * gamma_s;
    
    ClipDiagnostics clip_diag;
    VectorXd density_ratios = inference_density_ratios(density_logits, M_tau, clip_diag);
    
    // Calculate outcome predictions using W_outcome
    VectorXd outcome_preds = predict_glm_cpp(W_outcome, alpha_ts, family_int, link_int);
    
    // Calculate influence function components with numerical stability
    // Following eq:if_source_component / eq:dr_estimator in main.tex:
    // ξ̂_{t,s_j,i} = I(A_i=a) exp(-g_i) * (Y_i - ψ(X_i)) - δ̂.
    VectorXd xi_components(A_source.size());
    for (int i = 0; i < A_source.size(); i++) {
        if (A_source(i) == A_val) {
            // density_ratios = exp(-g), and this branch enforces A_i == A_val,
            // so ratio term is exactly 1/exp(g) for included units.
            double ratio_term = clipped_ratio_term(density_ratios(i), clip_diag);
            xi_components(i) = ratio_term * (Y_source(i) - outcome_preds(i)) - delta_ts;
        } else {
            xi_components(i) = -delta_ts;
        }
    }
    
    // Calculate variance: V̂_s = (1/N_s) Σ(ξ_{t,s_j,i} - ξ̄)^2 following eq:variance_components in main.tex (V_{s_j} component)
    double xi_mean = xi_components.mean();
    double V_s = (xi_components.array() - xi_mean).square().mean();
    
    // Ensure positive variance with more conservative bound
    V_s = std::max(V_s, NumericalConstants::VAR_MIN);
    
    return List::create(Named("V_s") = V_s,
                       Named("xi_components") = xi_components,
                       Named("clip_diagnostics") = make_clip_diagnostics(clip_diag));
}

// C++ version of target variance calculation
// family_int/link_int control the GLM response function for outcome predictions
// [[Rcpp::export]]
List calculate_target_variance_cpp(const MatrixXd& target_W_outcome, const VectorXd& alpha_ts,
                                  double M_ts,
                                  int family_int, int link_int) {
    
    int n = target_W_outcome.rows();
    
    if (n == 0) {
        return List::create(Named("V_t") = NumericalConstants::HESSIAN_FLOOR,
                           Named("zeta_components") = VectorXd::Zero(0));
    }
    
    // Calculate outcome predictions
    VectorXd outcome_preds = predict_glm_cpp(target_W_outcome, alpha_ts, family_int, link_int);
    
    // Calculate influence function components: ζ_{t,s_j,i} = h(φ(X_i); β̂) - Ẽ_t[h(φ(X); β̂)]
    VectorXd zeta_components = outcome_preds.array() - M_ts;
    
    // Calculate variance: V̂_{t,s_j} = (1/N_t) Σ(ζ_{t,s_j,i})^2 following eq:variance_components in main.tex (V_{t,s_j} component)
    // Note: mean(zeta_components) = 0 by construction, so no need to subtract mean again
    double V_t = zeta_components.array().square().mean();
    
    // Ensure positive variance with more conservative bound
    V_t = std::max(V_t, NumericalConstants::VAR_MIN);
    
    return List::create(Named("V_t") = V_t,
                       Named("zeta_components") = zeta_components);
}

// C++ version of covariance term calculation
// [[Rcpp::export]]
double calculate_covariance_term_cpp(const VectorXd& varphi_ot, const VectorXd& zeta_components) {
    
    int n = varphi_ot.size();
    
    if (n == 0 || zeta_components.size() == 0 || n != zeta_components.size()) {
        return 0.0;
    }
    
    // Calculate means
    double varphi_mean = varphi_ot.mean();
    // Note: zeta_components are already centered (mean = 0 by construction)
    double zeta_mean = zeta_components.mean(); // Should be ~0
    
    // Calculate covariance: Ĉ_{ot,s_j} = (1/N_t) Σ(φ_{ot,i} - φ̄)(ζ_{t,s_j,i} - ζ̄)
    double C_ot = ((varphi_ot.array() - varphi_mean) * (zeta_components.array() - zeta_mean)).mean();
    
    return C_ot;
}
