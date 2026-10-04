# Numeric model state for controlled comparisons; timing attributes are omitted.
rhc_stage_signatures <- function(fit) {
  coefficients <- function(value) if (is.null(value)) NULL else as.numeric(value)
  initial <- function(values) lapply(values, coefficients)
  rule <- function(value) attr(value, "lambda_rule")
  target <- messages <- initial_weights <- final_weights <- outcomes <- rules <- list()
  for (arm in c("mu1", "mu0")) for (fold in seq_len(fit$n_folds)) {
    record <- fit$arm_results[[arm]]$fold_results[[fold]]
    key <- paste(arm, fold, sep = ":")
    anchor <- record$target_only
    target[[key]] <- list(estimate = anchor$estimate, scores = as.numeric(anchor$varphi_ot),
      outcome_predictions = coefficients(anchor$m_pred),
      arm_propensity = coefficients(anchor$arm_propensity), propensity = coefficients(anchor$prop_scores),
      calibrated_outcome = coefficients(anchor$calibration$outcome),
      calibrated_weight = coefficients(anchor$calibration$weight))
    for (site in names(record$source_results)) {
      source <- record$source_results[[site]]
      label <- paste(key, site, sep = ":")
      messages[[label]] <- initial(source$per_k2_alpha)
      initial_weights[[label]] <- initial(source$per_k2_gamma)
      final_weights[[label]] <- coefficients(source$gamma_s)
      outcomes[[label]] <- coefficients(source$alpha_ts)
      rules[[label]] <- list(initial_weight = vapply(source$per_k2_gamma, rule, character(1L)),
                            weight = rule(source$gamma_s), outcome = rule(source$alpha_ts))
    }
  }
  list(target = target, initial_outcome_messages = messages, initial_weight = initial_weights,
       weight = final_weights, outcome = outcomes, selected_rules = rules)
}
