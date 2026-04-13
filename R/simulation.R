# simulation.R - Simulation experiment functions for FACE-C
#
# This file contains reusable simulation functions extracted from main.R.
# The orchestration script (main.R) calls these functions with experiment-
# specific configuration and checkpoint infrastructure.
#
# Contents:
#   1. run_single_simulation - Run one Monte Carlo replicate
#   2. summarize_results - Aggregate simulation results into summary statistics
#   3. run_simulation_study - Full study orchestration with checkpointing

#' Run a single Monte Carlo simulation replicate
#'
#' Generates data, runs the requested estimators, and returns a data frame
#' of per-method results (estimate, SE, bias, coverage, CI width).
#'
#' @param sim_id Integer simulation replicate ID (also used as RNG seed)
#' @param n_total Integer total sample size across all sites
#' @param K Integer number of source sites
#' @param p Integer number of covariates
#' @param config configuration ("C1", "C2", "C3", or "C4")
#'   - C1: Both models correctly specified
#'   - C2: Outcome model misspecified only
#'   - C3: Site/propensity model misspecified only
#'   - C4: Both models misspecified (tests double robustness boundary)
#' @param methods vector of methods to run
#' @param verbose Logical. Print progress messages via \code{log_info}.
#' @param n_cores_internal number of cores for internal parallelization (source sites).
#'        NULL or 1 for sequential, -1 for all cores minus 1.
#'        Note: When running simulations in parallel, set to 1 to avoid nested parallelism.
#' @param nlambda_init Integer. Number of lambda candidates for initial outcome
#'        model CV (glmnet). Lower values (e.g. 20) speed up fitting; default 100.
#' @param estimand_type Type of estimand for true value calculation:
#'        - "superpopulation" (default): fixed superpopulation parameter (same across all simulations)
#'        - "sample": sample-specific true value that varies across simulations
#' @param site_allocation Method for site allocation:
#'        - "model" (default): multinomial logistic model based on covariates (realistic)
#'        - "uniform": uniform random allocation (simpler, for variance validation)
#'        - "balanced": equal sample sizes per site (simplest, for baseline)
#' @param transform_type Type of covariate transformation:
#'        - "mild" (default): gentle transformations that preserve IF orthogonality
#'        - "strong": aggressive nonlinear transformations
#'        - "none": no transformation
#' @param outcome_type "binary" (logistic link) or "continuous" (identity link)
#' @param heterogeneity_type Level of outcome heterogeneity across sites:
#'        - "none": homogeneous outcome models
#'        - "mild": ~40%% change in 2 non-zero coefficients
#'        - "strong": ~80%% change in 2 non-zero coefficients
#'        - "partial": only first half of source sites have mild heterogeneity
#' @param shift_strength Numeric multiplier for covariate shift intensity (default 1.0)
#' @param n_folds Integer. Number of cross-fitting folds (>= 3, default 10)
#' @param use_lambda_cache Logical. If TRUE, enable lambda caching
#'   within cross-fitting nuisance-model loops.
#' @param estimate_ate Logical. If TRUE, also estimate target-only ATE
#' @return data frame with results
#' @export
run_single_simulation <- function(sim_id, n_total = 1000, K = 3, p = 4,
                                 config = "C1",
                                 methods = c("two_round_crossfit", "one_round_crossfit",
                                           "target_only", "sample_size", "inverse_variance",
                                           "federated_dr", "pooled_dr", "tilted_aipw",
                                           "oracle_dr"),
                                 verbose = TRUE, n_cores_internal = NULL,
                                 nlambda_init = LAMBDA_GRID_SIZE_STANDARD,
                                 estimand_type = "superpopulation",
                                 site_allocation = "model",
                                 transform_type = "mild",
                                 outcome_type = "binary",
                                 heterogeneity_type = "none",
                                 shift_strength = 1.0,
                                 n_folds = N_FOLDS_DEFAULT,
                                 use_lambda_cache = TRUE,
                                 estimate_ate = FALSE,
                                 # FACE paper DGP parameters
                                 dgp_type = "facec",
                                 ate_deviation    = 0.0,
                                 n_deviated_sites = 0L) {
  
  # Track timing for each stage
  sim_start_time <- Sys.time()
  stage_times <- list()
  
  log_info(verbose, "\n[%s] Simulation %d (config=%s, dgp=%s, heterogeneity=%s)\n",
           format(Sys.time(), "%H:%M:%S"), sim_id, config, dgp_type, heterogeneity_type)
  
  # Generate data
  data_gen_start <- Sys.time()
  set.seed(sim_id)
  data <- generate_simulation_data(n_total, K, p, config,
                                   estimand_type    = estimand_type,
                                   site_allocation  = site_allocation,
                                   transform_type   = transform_type,
                                   outcome_type     = outcome_type,
                                   heterogeneity_type = heterogeneity_type,
                                   shift_strength   = shift_strength,
                                   dgp_type         = dgp_type,
                                   ate_deviation    = ate_deviation,
                                   n_deviated_sites = as.integer(n_deviated_sites))
  data_split <- split_data_by_site(data)
  stage_times$data_gen <- as.numeric(difftime(Sys.time(), data_gen_start, units = "secs"))
  
  log_info(verbose, "    Data generated in %.2fs (n_t=%d, n_s=%s)\n", 
           stage_times$data_gen, data_split$t$n,
           paste(sapply(setdiff(names(data_split), "t"), function(s) data_split[[s]]$n), collapse=","))
  
  # Get true potential outcome from generated data
  # The interpretation depends on estimand_type:
  #   - "sample": sample-specific E_n[Y(1)] (varies across simulations)
  #   - "superpopulation": fixed E[Y(1)] (same for all simulations)
  true_potential_outcome <- data$mu1_true
  
  # Find target sample indices (for logging)
  target_idx <- which(data$R == "t")

  if (estimand_type == "superpopulation") {
    log_info(verbose, "  True potential outcome (superpop): %.4f, Sample-specific: %.4f\n", 
             true_potential_outcome, data$mu1_realized)
  } else {
    log_info(verbose, "  True potential outcome (sample): %.4f\n", true_potential_outcome)
  }
  log_info(verbose, "  Target sample size: %d, proportion: %.3f\n", 
           length(target_idx), length(target_idx) / n_total)
  
  # Initialize results as a list (avoid O(n²) rbind-in-loop)
  results_list <- list()
  crossfit_mu1_results <- list()

  # Precompute fold partition once and reuse across one/two-round and ATE reruns
  precomputed_folds <- build_crossfit_folds(data_split, n_folds)
  target_only_ps_cache <- new.env(hash = TRUE, parent = emptyenv())
  
  # Map outcome_type to GLM family for cross-fitting.
  # Use data$outcome_type rather than the argument so that the face DGP
  # (which always produces continuous outcomes) is handled correctly even when
  # the caller passes outcome_type = "binary" (the default).
  family <- switch(data$outcome_type,
    "binary"     = "binomial",
    "continuous" = "gaussian",
    tolower(data$outcome_type)
  )
  if (!family %in% VALID_GLM_FAMILIES) {
    stop(sprintf("Unsupported GLM family '%s'. Supported families: %s",
                 family, paste(VALID_GLM_FAMILIES, collapse = ", ")))
  }
  
  # Run one-round cross-fitting algorithm (if requested)
  if ("one_round_crossfit" %in% methods) {
    one_round_start <- Sys.time()
    log_info(verbose, "    Running one-round cross-fitting (nlambda=%d)...\n", nlambda_init)
    tryCatch({
      one_round_cf_res <- run_crossfit(data_split, n_folds = n_folds,
                  communication_mode = "one_round",
                  lambda_selection = "cv",
                  verbose = FALSE, n_cores = n_cores_internal,
                  nlambda_init = nlambda_init,
                  family = family,
                  use_lambda_cache = use_lambda_cache,
                  precomputed_folds = precomputed_folds,
                  target_only_ps_cache = target_only_ps_cache)
      stage_times$one_round <- as.numeric(difftime(Sys.time(), one_round_start, units = "secs"))
      
      log_info(verbose, "    [OK] One-round done in %.2fs: est=%.4f, bias=%.4f, SE=%.4f\n", 
               stage_times$one_round,
               one_round_cf_res$estimate, one_round_cf_res$estimate - true_potential_outcome,
               one_round_cf_res$se)
      crossfit_mu1_results[["one_round_crossfit"]] <- one_round_cf_res
      
      results_list[[length(results_list) + 1]] <- data.frame(
        sim_id = sim_id,
        method = "one_round_crossfit",
        estimate = one_round_cf_res$estimate,
        se = one_round_cf_res$se,
        bias = one_round_cf_res$estimate - true_potential_outcome,
        coverage = (true_potential_outcome >= one_round_cf_res$ci_lower) &
                  (true_potential_outcome <= one_round_cf_res$ci_upper),
        ci_width = one_round_cf_res$ci_upper - one_round_cf_res$ci_lower,
        n_total = n_total,
        K = K,
        p = p,
        config = config,
        heterogeneity_type = heterogeneity_type,
        estimand_type = estimand_type,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      stage_times$one_round <- as.numeric(difftime(Sys.time(), one_round_start, units = "secs"))
      log_info(verbose, "    [!!] One-round failed (%.2fs): %s\n", stage_times$one_round, e$message)
    })
  }
  
  # Run two-round cross-fitting algorithm
  # Note: n_folds >= 3 is required for proper two-level cross-fitting calibration
  if ("two_round_crossfit" %in% methods) {
    two_round_start <- Sys.time()
    log_info(verbose, "    Running two-round cross-fitting (nlambda=%d)...\n", nlambda_init)
    tryCatch({
      two_round_cf_res <- run_crossfit(data_split, n_folds = n_folds,
                  communication_mode = "two_round",
                  lambda_selection = "cv",
                  verbose = FALSE, n_cores = n_cores_internal,
                  nlambda_init = nlambda_init,
                  family = family,
                  use_lambda_cache = use_lambda_cache,
                  precomputed_folds = precomputed_folds,
                  target_only_ps_cache = target_only_ps_cache)
      stage_times$two_round <- as.numeric(difftime(Sys.time(), two_round_start, units = "secs"))
      
      log_info(verbose, "    [OK] Two-round done in %.2fs: est=%.4f, bias=%.4f, SE=%.4f\n", 
               stage_times$two_round,
               two_round_cf_res$estimate, two_round_cf_res$estimate - true_potential_outcome,
               two_round_cf_res$se)
      crossfit_mu1_results[["two_round_crossfit"]] <- two_round_cf_res
      
      results_list[[length(results_list) + 1]] <- data.frame(
        sim_id = sim_id,
        method = "two_round_crossfit",
        estimate = two_round_cf_res$estimate,
        se = two_round_cf_res$se,
        bias = two_round_cf_res$estimate - true_potential_outcome,
        coverage = (true_potential_outcome >= two_round_cf_res$ci_lower) &
                  (true_potential_outcome <= two_round_cf_res$ci_upper),
        ci_width = two_round_cf_res$ci_upper - two_round_cf_res$ci_lower,
        n_total = n_total,
        K = K,
        p = p,
        config = config,
        heterogeneity_type = heterogeneity_type,
        estimand_type = estimand_type,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      stage_times$two_round <- as.numeric(difftime(Sys.time(), two_round_start, units = "secs"))
      log_info(verbose, "    [!!] Two-round failed (%.2fs): %s\n", stage_times$two_round, e$message)
    })
  }
  
  # Run comparison methods
  comp_start <- Sys.time()
  log_info(verbose, "    Running comparison methods...\n")
  comparison_results <- tryCatch({
    run_all_comparisons(data_split, family = family)
  }, error = function(e) {
    log_info(verbose, "    [!!] All comparisons failed: %s\n", e$message)
    list()  # Return empty list so individual method checks below gracefully skip
  })
  stage_times$comparison <- as.numeric(difftime(Sys.time(), comp_start, units = "secs"))
  
  for (method_name in c("target_only", "sample_size", "inverse_variance", "federated_dr", "pooled_dr", "tilted_aipw")) {
    if (method_name %in% methods) {
      tryCatch({
        method_res <- comparison_results[[method_name]]
        
        # Calculate confidence interval
        ci_lower <- method_res$estimate - Z_ALPHA_05 * method_res$se
        ci_upper <- method_res$estimate + Z_ALPHA_05 * method_res$se
        
        results_list[[length(results_list) + 1]] <- data.frame(
          sim_id = sim_id,
          method = method_name,
          estimate = method_res$estimate,
          se = method_res$se,
          bias = method_res$estimate - true_potential_outcome,
          coverage = (true_potential_outcome >= ci_lower) & (true_potential_outcome <= ci_upper),
          ci_width = ci_upper - ci_lower,
          n_total = n_total,
          K = K,
          p = p,
          config = config,
          heterogeneity_type = heterogeneity_type,
          estimand_type = estimand_type,
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        log_info(verbose, "    [!!] %s failed: %s\n", method_name, e$message)
      })
    }
  }
  
  # Run oracle DR estimator (uses known true parameters)
  # Not available for the FACE paper DGP (gamma_params and alpha1_true are NULL).
  if ("oracle_dr" %in% methods) {
    oracle_start <- Sys.time()
    if (is.null(data$gamma_params) || is.null(data$alpha1_true)) {
      log_info(verbose, "    [--] Oracle DR skipped: true parameters not available for dgp_type='%s'\n",
               data$dgp_type %||% dgp_type)
    } else {
    log_info(verbose, "    Running oracle DR estimator...\n")
    tryCatch({
      target_propensity_true <- NULL
      if (!is.null(data$p_treat_true)) {
        target_propensity_true <- data$p_treat_true[target_idx]
      }
      oracle_res <- estimate_oracle_dr(data_split, data$alpha1_true, data$gamma_params,
                                       outcome_type = data$outcome_type,
                                       target_propensity_true = target_propensity_true)
      stage_times$oracle <- as.numeric(difftime(Sys.time(), oracle_start, units = "secs"))
      
      ci_lower <- oracle_res$estimate - Z_ALPHA_05 * oracle_res$se
      ci_upper <- oracle_res$estimate + Z_ALPHA_05 * oracle_res$se
      
      log_info(verbose, "    [OK] Oracle DR done in %.2fs: est=%.4f, bias=%.4f\n",
               stage_times$oracle, oracle_res$estimate,
               oracle_res$estimate - true_potential_outcome)
      
      results_list[[length(results_list) + 1]] <- data.frame(
        sim_id = sim_id,
        method = "oracle_dr",
        estimate = oracle_res$estimate,
        se = oracle_res$se,
        bias = oracle_res$estimate - true_potential_outcome,
        coverage = (true_potential_outcome >= ci_lower) & (true_potential_outcome <= ci_upper),
        ci_width = ci_upper - ci_lower,
        n_total = n_total,
        K = K,
        p = p,
        config = config,
        heterogeneity_type = heterogeneity_type,
        estimand_type = estimand_type,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      log_info(verbose, "    [!!] Oracle DR failed: %s\n", e$message)
    })
    }  # end else (gamma_params guard)
  }
  
  # ATE estimation: rerun key methods with A_val=0 and compute difference
  if (isTRUE(estimate_ate)) {
    log_info(verbose, "    Running ATE estimation (mu0 arm)...\n")
    ate_start <- Sys.time()
    
    # Get mu0 true value
    mu0_true <- data$mu0_true
    
    # Cross-fitted ATE: rerun cross-fitting algorithms with A_val=0 for mu0
    ate_true <- data$mu1_true - mu0_true
    for (ate_method in c("two_round_crossfit", "one_round_crossfit")) {
      if (ate_method %in% methods) {
        tryCatch({
          comm_mode <- if (ate_method == "two_round_crossfit") "two_round" else "one_round"
          mu0_res <- run_crossfit(data_split, n_folds = n_folds,
                                  communication_mode = comm_mode,
                                  lambda_selection = "cv",
                                  verbose = FALSE, n_cores = n_cores_internal,
                                  nlambda_init = nlambda_init, family = family, A_val = 0L,
                                  use_lambda_cache = use_lambda_cache,
                                  precomputed_folds = precomputed_folds,
                                  target_only_ps_cache = target_only_ps_cache)
          
          # Find the corresponding mu1 estimate from results_list
          mu1_method_label <- ate_method
          mu1_row <- which(sapply(results_list, function(x) x$method[1] == mu1_method_label))
          if (length(mu1_row) > 0) {
            mu1_est <- results_list[[mu1_row[1]]]$estimate[1]
            mu1_se <- results_list[[mu1_row[1]]]$se[1]
            
            ate_est <- mu1_est - mu0_res$estimate
            mu1_cf <- crossfit_mu1_results[[ate_method]]
            mu1_phi <- if (!is.null(mu1_cf) && !is.null(mu1_cf$all_phi_agg)) mu1_cf$all_phi_agg else NULL
            mu0_phi <- if (!is.null(mu0_res$all_phi_agg)) mu0_res$all_phi_agg else NULL

            if (!is.null(mu1_phi) && !is.null(mu0_phi) && length(mu1_phi) == length(mu0_phi)) {
              # Joint IF variance: Var(ATE) = E[(Phi_mu1 - Phi_mu0)^2] / N
              ate_var <- mean((mu1_phi - mu0_phi)^2) / length(mu1_phi)
            } else {
              # Fallback: conservative independence approximation
              ate_var <- mu1_se^2 + mu0_res$se^2
            }
            ate_se <- sqrt(max(ate_var, 0))
            ci_lower <- ate_est - Z_ALPHA_05 * ate_se
            ci_upper <- ate_est + Z_ALPHA_05 * ate_se
            
            results_list[[length(results_list) + 1]] <- data.frame(
              sim_id = sim_id,
              method = paste0(ate_method, "_ate"),
              estimate = ate_est,
              se = ate_se,
              bias = ate_est - ate_true,
              coverage = (ate_true >= ci_lower) & (ate_true <= ci_upper),
              ci_width = ci_upper - ci_lower,
              n_total = n_total,
              K = K,
              p = p,
              config = config,
              heterogeneity_type = heterogeneity_type,
              estimand_type = estimand_type,
              stringsAsFactors = FALSE
            )
            log_info(verbose, "    [OK] %s ATE: est=%.4f, true=%.4f, bias=%.4f\n",
                     ate_method, ate_est, ate_true, ate_est - ate_true)
          }
        }, error = function(e) {
          log_info(verbose, "    [!!] ATE %s failed: %s\n", ate_method, e$message)
        })
      }
    }
    
    # Target-only ATE: joint IF variance (not assuming independence)
    tryCatch({
      target_data <- data_split[["t"]]
      mu0_cf <- estimate_target_only_crossfit(target_data, n_folds = n_folds,
                                               family = family,
                                               A_val = 0)
      # Find the mu1 target-only estimate from results_list
      mu1_target_row <- which(sapply(results_list, function(x) x$method[1] == "target_only"))
      if (length(mu1_target_row) > 0) {
        mu1_est <- results_list[[mu1_target_row[1]]]$estimate[1]
        
        ate_est <- mu1_est - mu0_cf$estimate
        ate_true <- data$mu1_true - mu0_true
        
        # Use joint IF variance: Var(ATE) = E[(psi1_i - psi0_i)^2] / n
        # varphi_ot from each arm is already centered (mean zero).
        # comparison_results[["target_only"]]$varphi_ot holds mu1 IFs.
        mu1_psi <- comparison_results[["target_only"]]$varphi_ot
        mu0_psi <- mu0_cf$varphi_ot
        
        if (!is.null(mu1_psi) && !is.null(mu0_psi) &&
            length(mu1_psi) == length(mu0_psi)) {
          # Joint IF: accounts for covariance between mu1 and mu0 estimates
          ate_var <- mean((mu1_psi - mu0_psi)^2) / length(mu1_psi)
        } else {
          # Fallback: conservative (independence) approximation
          mu1_se <- results_list[[mu1_target_row[1]]]$se[1]
          ate_var <- mu1_se^2 + mu0_cf$variance
        }
        ate_se <- sqrt(ate_var)
        ci_lower <- ate_est - Z_ALPHA_05 * ate_se
        ci_upper <- ate_est + Z_ALPHA_05 * ate_se
        
        results_list[[length(results_list) + 1]] <- data.frame(
          sim_id = sim_id,
          method = "target_only_ate",
          estimate = ate_est,
          se = ate_se,
          bias = ate_est - ate_true,
          coverage = (ate_true >= ci_lower) & (ate_true <= ci_upper),
          ci_width = ci_upper - ci_lower,
          n_total = n_total,
          K = K,
          p = p,
          config = config,
          heterogeneity_type = heterogeneity_type,
          estimand_type = estimand_type,
          stringsAsFactors = FALSE
        )
        log_info(verbose, "    [OK] Target-only ATE: est=%.4f, true=%.4f, bias=%.4f\n",
                 ate_est, ate_true, ate_est - ate_true)
      }
    }, error = function(e) {
      log_info(verbose, "    [!!] ATE estimation failed: %s\n", e$message)
    })
    
    stage_times$ate <- as.numeric(difftime(Sys.time(), ate_start, units = "secs"))
  }
  
  # Print total simulation time summary
  sim_total_time <- as.numeric(difftime(Sys.time(), sim_start_time, units = "secs"))
  log_info(verbose, "    [OK] Comparison methods done in %.2fs\n", stage_times$comparison)
  log_info(verbose, "    [TIME] Sim %d total: %.2fs (data=%.1f%%, 1rnd=%.1f%%, 2rnd=%.1f%%, comp=%.1f%%)\n",
           sim_id, sim_total_time,
           100 * (stage_times$data_gen %||% 0) / sim_total_time,
           100 * (stage_times$one_round %||% 0) / sim_total_time,
           100 * (stage_times$two_round %||% 0) / sim_total_time,
           100 * (stage_times$comparison %||% 0) / sim_total_time)
  
  # Combine results from list into a data.frame (single rbind at the end)
  results <- if (length(results_list) > 0) do.call(rbind, results_list) else data.frame()

  # Stamp dgp_type onto every row so downstream summaries can group by it
  if (nrow(results) > 0) results$dgp_type <- dgp_type

  return(results)
}

#' Summarize simulation results into aggregate statistics
#'
#' Computes mean bias, standard deviation of bias, RMSE, mean SE,
#' coverage probability, and CI width statistics, grouped by method
#' and configuration.
#'
#' @param results Data frame of per-simulation results (output of
#'   \code{run_single_simulation}). Must contain columns: method, config,
#'   heterogeneity_type, n_total, K, bias, se, coverage, ci_width.
#' @return Data frame of summary statistics with one row per
#'   method-config-heterogeneity-n_total-K combination.
#' @export
summarize_results <- function(results) {
  # Handle empty or invalid results
  if (is.null(results) || nrow(results) == 0) {
    cat("[WARN] No results to summarize (empty data frame)\n")
    return(data.frame())
  }
  
  # Check required columns exist
  required_cols <- c("method", "config", "heterogeneity_type", "n_total", "K", "bias", "se", "coverage", "ci_width")
  missing_cols <- setdiff(required_cols, names(results))
  if (length(missing_cols) > 0) {
    cat(sprintf("[WARN] Missing columns in results: %s\n", paste(missing_cols, collapse = ", ")))
    return(data.frame())
  }
  
  # Determine grouping columns (include estimand_type and dgp_type if present)
  group_cols <- c("method", "config", "heterogeneity_type", "n_total", "K")
  if ("estimand_type" %in% names(results)) {
    group_cols <- c(group_cols, "estimand_type")
  }
  if ("dgp_type" %in% names(results)) {
    group_cols <- c(group_cols, "dgp_type")
  }
  group_formula <- as.formula(paste("~ ", paste(group_cols, collapse = " + ")))
  
  # Calculate summary statistics by method and configuration
  # For bias and se, we need mean, sd, and rmse
  bias_stats <- aggregate(
    update(group_formula, bias ~ .),
    data = results,
    FUN = function(x) c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      rmse = sqrt(mean(x^2, na.rm = TRUE))  # RMSE for bias
    )
  )
  
  # For se (IF-based), we only need mean
  se_stats <- aggregate(
    update(group_formula, se ~ .),
    data = results,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  
  # For coverage (IF-based), we only need mean (proportion of TRUE values)
  coverage_stats <- aggregate(
    update(group_formula, coverage ~ .),
    data = results,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  
  # For CI width (IF-based), we need mean and sd
  ci_width_stats <- aggregate(
    update(group_formula, ci_width ~ .),
    data = results,
    FUN = function(x) c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE)
    )
  )
  
  # Build base summary data frame
  summary_df <- data.frame(
    method = bias_stats$method,
    config = bias_stats$config,
    heterogeneity_type = bias_stats$heterogeneity_type,
    n_total = bias_stats$n_total,
    K = bias_stats$K,
    bias_mean = bias_stats$bias[, "mean"],
    bias_sd = bias_stats$bias[, "sd"],
    rmse = bias_stats$bias[, "rmse"],
    se_mean = se_stats$se,
    coverage = coverage_stats$coverage,
    ci_width_mean = ci_width_stats$ci_width[, "mean"],
    ci_width_sd = ci_width_stats$ci_width[, "sd"]
  )
  
  # Add estimand_type column if present in results
  if ("estimand_type" %in% names(bias_stats)) {
    summary_df$estimand_type <- bias_stats$estimand_type
  }
  # Add dgp_type column if present in results
  if ("dgp_type" %in% names(bias_stats)) {
    summary_df$dgp_type <- bias_stats$dgp_type
  }

  return(summary_df)
}


# =============================================================================
# SIMULATION STUDY ORCHESTRATION
# =============================================================================

#' Run a full simulation study with checkpoint/restart support
#'
#' Iterates over a parameter grid (n, K, p, config), running
#' \code{n_sims} Monte Carlo replicates per setting. Supports
#' SLURM preemption-safe checkpointing and nested parallelism.
#'
#' @param n_sims Integer number of MC replicates per setting.
#' @param n_total_vec Integer vector of sample sizes.
#' @param K_vec Integer vector of source-site counts.
#' @param p_vec Integer vector of covariate dimensions.
#' @param configs Character vector of configurations (e.g., "C1").
#' @param n_cores Integer number of available cores.
#' @param checkpoint_config List returned by \code{init_checkpoint_config()},
#'   or NULL to disable checkpointing. Fields: \code{file},
#'   \code{setting_interval}, \code{sim_interval}, \code{preempt_signal},
#'   \code{saved_signal}.
#' @param nlambda_init Integer lambda grid size for initial outcome CV.
#' @param nested_parallel Logical. Enable sim + source-site parallelism.
#' @param estimand_type "superpopulation" or "sample".
#' @param site_allocation "model", "uniform", or "balanced".
#' @param transform_type "strong", "mild", or "none".
#' @param outcome_type "binary" or "continuous".
#' @param heterogeneity_type "none", "mild", "strong", or "partial".
#' @param shift_strength Numeric covariate shift multiplier.
#' @param n_folds Integer cross-fitting fold count.
#' @param use_lambda_cache Logical. If TRUE, enable lambda caching
#'   within cross-fitting nuisance-model loops.
#' @param verbose_every Integer. In sequential mode, run detailed simulation
#'   logging every \\code{verbose_every} simulations (default 10).
#' @param parallel_strategy Character. Nested parallel allocation strategy:
#'   \code{"outer_priority"} (default), \code{"balanced"}, or
#'   \code{"outer_only"}.
#' @param estimate_ate Logical. Also estimate ATE via A_val=0?
#' @param dgp_type "facec" or "face".
#' @param ate_deviation Numeric ATE deviation (FACE DGP only).
#' @param n_deviated_sites Integer deviated sites (FACE DGP only).
#' @return Data frame of simulation results.
#' @export
run_simulation_study <- function(n_sims = 500,
                                n_total_vec = c(500, 1000, 2000),
                                K_vec = c(2, 3, 5),
                                p_vec = c(4, 8),
                                configs = c("C1", "C2", "C3", "C4"),
                                n_cores = 1,
                                checkpoint_config = NULL,
                                nlambda_init = LAMBDA_GRID_SIZE_STANDARD,
                                nested_parallel = FALSE,
                                estimand_type = "superpopulation",
                                site_allocation = "model",
                                transform_type = "mild",
                                outcome_type = "binary",
                                heterogeneity_type = "none",
                                shift_strength = 1.0,
                                n_folds = N_FOLDS_DEFAULT,
                                use_lambda_cache = TRUE,
                                verbose_every = 10L,
                                parallel_strategy = c("outer_priority", "balanced", "outer_only"),
                                estimate_ate = FALSE,
                                dgp_type = "facec",
                                ate_deviation    = 0.0,
                                n_deviated_sites = 0L) {

  parallel_strategy <- match.arg(parallel_strategy)

  # Resolve checkpoint config (NULL disables checkpointing)
  ckpt_file         <- checkpoint_config$file
  ckpt_interval     <- checkpoint_config$setting_interval %||% 1L
  sim_ckpt_interval <- checkpoint_config$sim_interval %||% n_sims
  preempt_file      <- checkpoint_config$preempt_signal
  saved_file        <- checkpoint_config$saved_signal
  use_checkpoint    <- !is.null(ckpt_file)

  # Create parameter grid
  param_grid <- expand.grid(
    n_total = n_total_vec,
    K = K_vec,
    p = p_vec,
    config = configs,
    stringsAsFactors = FALSE
  )
  total_settings <- nrow(param_grid)

  cat(sprintf("\n[FACE-C] Simulation study: %d settings, %d sims each, %d cores\n",
              total_settings, n_sims, n_cores))
  cat(sprintf("  n_total: %s | K: %s | p: %s | configs: %s\n",
              paste(n_total_vec, collapse = ","),
              paste(K_vec, collapse = ","),
              paste(p_vec, collapse = ","),
              paste(configs, collapse = ",")))
  cat(sprintf("  estimand: %s | outcome: %s | heterogeneity: %s\n",
              estimand_type, outcome_type, heterogeneity_type))
  cat(sprintf("  checkpoint: %s (setting every %d, sim every %d)\n",
              if (use_checkpoint) "ON" else "OFF", ckpt_interval, sim_ckpt_interval))

  # Try to load checkpoint
  checkpoint <- if (use_checkpoint) load_checkpoint(ckpt_file) else NULL

  # Variables for simulation-level resume
  resume_sim_results <- NULL
  resume_sim_idx <- 0

  if (!is.null(checkpoint)) {
    all_results <- checkpoint$all_results

    if (!is.null(checkpoint$sim_results) && !is.null(checkpoint$current_sim_idx)) {
      start_idx <- checkpoint$current_setting_idx
      resume_sim_results <- checkpoint$sim_results
      resume_sim_idx <- checkpoint$current_sim_idx
      cat(sprintf("  Resuming setting %d from simulation %d\n",
                  start_idx, resume_sim_idx + 1))
    } else {
      start_idx <- checkpoint$current_setting_idx + 1
      cat(sprintf("  Resuming from setting %d/%d\n", start_idx, total_settings))
    }

    if (checkpoint$current_setting_idx >= total_settings &&
        is.null(checkpoint$sim_results)) {
      cat("  Checkpoint indicates all settings complete. Finalizing.\n")
      final_results <- do.call(rbind, all_results)
      if (use_checkpoint) cleanup_checkpoint(ckpt_file, preempt_file, saved_file)
      return(final_results)
    }
  } else {
    all_results <- list()
    start_idx <- 1
    cat("  Starting fresh (no checkpoint found)\n")
  }

  # Set up parallel processing
  cl <- NULL
  n_cores_internal_parallel <- 1
  n_cores_outer <- n_cores

  if (n_cores > 1) {
    if (nested_parallel) {
      max_k <- max(K_vec)

      if (parallel_strategy == "outer_only" || max_k <= 1) {
        outer_cores <- n_cores
        inner_cores <- 1L
      } else if (parallel_strategy == "outer_priority") {
        # High-core default: prioritize simulation-level throughput.
        # Keep inner parallelism minimal (<=2) to avoid nested overhead.
        inner_cores <- if (n_cores >= 16 && max_k >= 2) min(2L, max_k) else 1L
        outer_cores <- max(1L, floor(n_cores / inner_cores))
      } else {
        # balanced
        inner_cores <- min(max_k, n_cores)
        outer_cores <- max(2L, floor(n_cores / inner_cores))
        inner_cores <- max(1L, floor(n_cores / outer_cores))
        inner_cores <- min(inner_cores, max_k)
      }

      if (outer_cores < 2L) {
        outer_cores <- n_cores
        inner_cores <- 1L
      }

      n_cores_outer <- outer_cores
      n_cores_internal_parallel <- inner_cores

      cat(sprintf("  Nested parallel (%s): %d outer x %d inner (of %d total)\n",
                  parallel_strategy, outer_cores, inner_cores, n_cores))
    } else {
      n_cores_outer <- n_cores
      cat(sprintf("  Using %d cores for simulation-level parallelization\n", n_cores))
    }

    cl <- parallel::makeCluster(n_cores_outer)
    doParallel::registerDoParallel(cl)

    parallel::clusterEvalQ(cl, {
      if (requireNamespace("FACEC", quietly = TRUE)) {
        library(FACEC)
      } else if (requireNamespace("devtools", quietly = TRUE)) {
        devtools::load_all(".")
      } else {
        stop("Workers require either installed FACEC package or devtools.")
      }
    })

    parallel::clusterExport(cl, c(
      "run_single_simulation",
      "n_cores_internal_parallel", "nlambda_init", "estimand_type",
      "site_allocation", "transform_type", "outcome_type",
      "heterogeneity_type", "shift_strength", "n_folds",
      "estimate_ate", "dgp_type", "ate_deviation", "n_deviated_sites"
    ), envir = environment())
  }

  # Main loop over parameter settings
  for (i in start_idx:total_settings) {

    # Check for preemption signal
    if (use_checkpoint && check_preempt_signal(preempt_file)) {
      cat(sprintf("\n[WARN] Preemption signal at setting %d/%d. Saving checkpoint...\n",
                  i, total_settings))
      state <- list(
        all_results = all_results,
        current_setting_idx = i - 1,
        total_settings = total_settings,
        timestamp = Sys.time()
      )
      checkpoint_saved <- save_checkpoint(state, ckpt_file)
      if (checkpoint_saved) signal_checkpoint_saved(saved_file)
      if (!is.null(cl)) parallel::stopCluster(cl)
      cat("  Checkpoint saved. Exiting for requeue.\n")
      quit(save = "no", status = 0)
    }

    params <- param_grid[i, ]
    cat(sprintf("\n[%s] Setting %d/%d: n=%d, K=%d, p=%d, config=%s\n",
                format(Sys.time(), "%H:%M:%S"), i, total_settings,
                params$n_total, params$K, params$p, params$config))

    methods <- c("two_round_crossfit", "one_round_crossfit",
                 "target_only", "sample_size", "inverse_variance",
                 "federated_dr", "pooled_dr", "tilted_aipw",
                 "oracle_dr")

    # Resume data if available
    if (i == start_idx && !is.null(resume_sim_results)) {
      sim_results <- resume_sim_results
      sim_start_idx <- resume_sim_idx + 1
      cat(sprintf("  Resuming with %d completed simulations\n", resume_sim_idx))
      resume_sim_results <- NULL
      resume_sim_idx <- 0
    } else {
      sim_results <- list()
      sim_start_idx <- 1
    }

    setting_start_time <- Sys.time()
    completed_sims_timing <- numeric(0)

    tryCatch({
      if (n_cores_outer > 1) {
        # --- Parallel execution (batched) ---
        batch_size <- n_cores_outer * 2
        first_batch <- ceiling(sim_start_idx / batch_size)
        n_batches <- ceiling(n_sims / batch_size)

        for (batch in first_batch:n_batches) {
          batch_start_time <- Sys.time()
          start_sim <- max((batch - 1) * batch_size + 1, sim_start_idx)
          end_sim <- min(batch * batch_size, n_sims)
          batch_ids <- start_sim:end_sim

          elapsed_secs <- as.numeric(difftime(Sys.time(), setting_start_time, units = "secs"))
          sims_done <- start_sim - sim_start_idx
          eta_str <- if (sims_done > 0) {
            sprintf("ETA: %.1f min",
                    elapsed_secs / sims_done * (n_sims - start_sim + 1) / 60)
          } else {
            "ETA: calculating..."
          }

          pct <- end_sim / n_sims * 100
          bar_width <- 20
          filled <- floor(bar_width * pct / 100)
          bar <- paste0("[", strrep("\u2588", filled),
                        strrep("\u2591", bar_width - filled), "]")
          cat(sprintf("  [%s] %s %.0f%% | Sims %d-%d/%d | Batch %d/%d | %s\n",
                      format(Sys.time(), "%H:%M:%S"), bar, pct,
                      start_sim, end_sim, n_sims, batch, n_batches, eta_str))

          batch_results <- parallel::parLapply(cl, batch_ids, function(sim_id) {
            run_single_simulation(
              sim_id, params$n_total, params$K, params$p,
              params$config, methods,
              verbose = FALSE,
              n_cores_internal = n_cores_internal_parallel,
              nlambda_init = nlambda_init,
              estimand_type = estimand_type,
              site_allocation = site_allocation,
              transform_type = transform_type,
              outcome_type = outcome_type,
              heterogeneity_type = heterogeneity_type,
              shift_strength = shift_strength,
              n_folds = n_folds,
              use_lambda_cache = use_lambda_cache,
              estimate_ate = estimate_ate,
              dgp_type = dgp_type,
              ate_deviation = ate_deviation,
              n_deviated_sites = n_deviated_sites)
          })

          batch_elapsed <- as.numeric(
            difftime(Sys.time(), batch_start_time, units = "secs"))
          n_success <- sum(vapply(batch_results,
                                 function(r) nrow(r) > 0, logical(1)))
          cat(sprintf("    Batch done: %.1fs (%.2fs/sim) | %d/%d succeeded\n",
                      batch_elapsed, batch_elapsed / length(batch_ids),
                      n_success, length(batch_ids)))

          for (j in seq_along(batch_results)) {
            sim_results[[batch_ids[j]]] <- batch_results[[j]]
          }

          # Sim-level checkpoint
          if (use_checkpoint &&
              end_sim %% sim_ckpt_interval == 0 && end_sim < n_sims) {
            state <- list(
              all_results = all_results, current_setting_idx = i,
              total_settings = total_settings, sim_results = sim_results,
              current_sim_idx = end_sim, n_sims = n_sims,
              timestamp = Sys.time()
            )
            save_checkpoint(state, ckpt_file, sim_level = TRUE)
          }
        }
      } else {
        # --- Sequential execution ---
        n_cores_internal <- min(params$K,
                               max(1, parallel::detectCores() - 1))
        if (n_cores_internal > 1) {
          cat(sprintf("  Using %d internal cores for source-site parallelization\n",
                      n_cores_internal))
        }

        for (sim_id in sim_start_idx:n_sims) {
          sim_start <- Sys.time()

          sims_done <- sim_id - sim_start_idx
          eta_str <- if (sims_done > 0 && length(completed_sims_timing) > 0) {
            sprintf("ETA: %.1f min",
                    mean(completed_sims_timing) * (n_sims - sim_id + 1) / 60)
          } else {
            "ETA: calculating..."
          }

          if (sim_id %% 10 == 1 || sim_id == sim_start_idx) {
            cat(sprintf("  [%s] Progress: sim %d/%d (%.0f%%) | %s\n",
                        format(Sys.time(), "%H:%M:%S"),
                        sim_id, n_sims, sim_id / n_sims * 100, eta_str))
          }

          sim_results[[sim_id]] <- run_single_simulation(
            sim_id, params$n_total, params$K, params$p,
            params$config, methods,
            verbose = (sim_id == sim_start_idx) || ((sim_id - sim_start_idx) %% max(1L, as.integer(verbose_every)) == 0L),
            n_cores_internal = n_cores_internal,
            nlambda_init = nlambda_init,
            estimand_type = estimand_type,
            site_allocation = site_allocation,
            transform_type = transform_type,
            outcome_type = outcome_type,
            heterogeneity_type = heterogeneity_type,
            shift_strength = shift_strength,
            n_folds = n_folds,
            use_lambda_cache = use_lambda_cache,
            estimate_ate = estimate_ate,
            dgp_type = dgp_type,
            ate_deviation = ate_deviation,
            n_deviated_sites = n_deviated_sites
          )

          sim_elapsed <- as.numeric(difftime(Sys.time(), sim_start, units = "secs"))
          completed_sims_timing <- c(completed_sims_timing, sim_elapsed)
          if (length(completed_sims_timing) > 10) {
            completed_sims_timing <- tail(completed_sims_timing, 10)
          }

          # Sim-level checkpoint
          if (use_checkpoint &&
              sim_id %% sim_ckpt_interval == 0 && sim_id < n_sims) {
            state <- list(
              all_results = all_results, current_setting_idx = i,
              total_settings = total_settings, sim_results = sim_results,
              current_sim_idx = sim_id, n_sims = n_sims,
              timestamp = Sys.time()
            )
            save_checkpoint(state, ckpt_file, sim_level = TRUE)
          }
        }
      }

      setting_results <- do.call(rbind, sim_results)
      all_results[[i]] <- setting_results

      setting_elapsed <- as.numeric(
        difftime(Sys.time(), setting_start_time, units = "secs"))
      cat(sprintf("  Setting %d/%d done: %d rows, %.1f min (%.2f s/sim)\n",
                  i, total_settings, nrow(setting_results),
                  setting_elapsed / 60, setting_elapsed / n_sims))

    }, error = function(e) {
      cat(sprintf("  [ERROR] Setting %d failed: %s\n", i, e$message))
      all_results[[i]] <<- NULL
    })

    # Setting-level checkpoint
    if (use_checkpoint && (i %% ckpt_interval == 0 || i == total_settings)) {
      state <- list(
        all_results = all_results, current_setting_idx = i,
        total_settings = total_settings, timestamp = Sys.time()
      )
      save_checkpoint(state, ckpt_file)
    }
  }

  if (!is.null(cl)) parallel::stopCluster(cl)

  final_results <- do.call(rbind, all_results)

  if (is.null(final_results) || nrow(final_results) == 0) {
    cat("\n[WARN] No simulation results collected!\n")
    final_results <- data.frame()
  }

  if (use_checkpoint) cleanup_checkpoint(ckpt_file, preempt_file, saved_file)

  total_rows <- if (is.data.frame(final_results)) nrow(final_results) else 0
  cat(sprintf("\nSimulation study completed: %d rows, %d settings\n",
              total_rows, total_settings))

  return(final_results)
}
