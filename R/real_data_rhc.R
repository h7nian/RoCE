# ============================================================================
# real_data_rhc.R — Loader and cohort builder for the Right-Heart
#   Catheterization (RHC) public teaching dataset.
# ============================================================================
# Data source:
#   Connors AF Jr, Speroff T, Dawson NV, et al.
#   "The effectiveness of right heart catheterization in the initial care of
#   critically ill patients." JAMA. 1996;276(11):889-897.
#
#   Teaching mirror curated by Frank Harrell at
#   https://hbiostat.org/data/repo/rhc.csv  (vendored as inst/extdata/rhc.csv).
#
# The dataset is fully public (no DUA / no credentialing); the CSV ships with
# the RoCE package so downstream users do not need to download anything.
# ============================================================================

# ----------------------------------------------------------------------------
# Historical observed-variable profile motivated by Connors et al. (1996)
# and Hirano & Imbens (2001), Health Services and Outcomes Research Methodology.
# This is not an exact reconstruction of their published 72-term design.
# The selected site variable is excluded from the covariate vector.
# ----------------------------------------------------------------------------

.RHC_CONTINUOUS_VARS <- c(
  "age", "edu", "das2d3pc", "surv2md1", "aps1", "scoma1",
  "wtkilo1", "temp1", "meanbp1", "resp1", "hrt1",
  "pafi1", "paco21", "ph1", "wblc1", "hema1",
  "sod1", "pot1", "crea1", "bili1", "alb1", "urin1"
)

.RHC_BINARY_VARS <- c(
  # Prior-comorbidity history
  "cardiohx", "chfhx", "dementhx", "psychhx", "chrpulhx",
  "renalhx", "liverhx", "gibledhx", "malighx", "immunhx",
  "transhx", "amihx",
  # Admission-system diagnosis flags (10 organ systems, "Yes"/"No" strings
  # in the raw CSV; the binary conversion below coerces them to 0/1).
  "resp", "card", "neuro", "gastr", "renal", "meta", "hema",
  "seps", "trauma", "ortho"
)

# cat1 (9 primary-disease levels) is included as a categorical covariate by
# default. When the user selects site_var = "cat1", build_rhc_cohort()
# automatically removes it from the covariate set (via setdiff) so that the
# site label never leaks into the predictors.
.RHC_CATEGORICAL_VARS <- c("sex", "race", "income", "ninsclas", "ca", "dnr1", "cat1")

.RHC_SITE_VAR      <- "cat1"
.RHC_TREATMENT_VAR <- "swang1"


#' Path to the Vendored RHC CSV
#'
#' Returns the absolute path to \code{inst/extdata/rhc.csv}. Useful for
#' diagnostics and reproducibility; end users normally call
#' \code{\link{load_rhc_raw}} instead.
#'
#' @return Character. Absolute path to the CSV (zero-length string if the
#'   file has not yet been downloaded via \code{data/download_rhc.sh}).
#' @export
rhc_csv_path <- function() {
  system.file("extdata", "rhc.csv", package = "RoCE")
}


#' Load the Raw RHC Dataset
#'
#' Reads the vendored \code{inst/extdata/rhc.csv} into a data frame with all
#' 63 original columns intact. No coercion is applied; downstream processing
#' is delegated to \code{\link{build_rhc_cohort}}.
#'
#' @return A data frame with 5,735 rows.
#' @export
load_rhc_raw <- function() {
  path <- rhc_csv_path()
  if (!nzchar(path) || !file.exists(path)) {
    stop(
      "RHC CSV not found. Run `./data/download_rhc.sh` to fetch it from ",
      "https://hbiostat.org/data/repo/rhc.csv into inst/extdata/rhc.csv, ",
      "then reinstall the package."
    )
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}


# Median-impute numeric variables; mode-impute categorical variables. Both are
# computed on the cohort as a whole by default, or on explicit training rows.
# Both helpers fail fast when no valid value is available
# to use as a fill, since silently returning all-NA / unchanged data would
# defer the failure to a far less informative location downstream.
.rhc_impute_continuous <- function(x, var_name = "<unknown>", training_rows = NULL) {
  if (!is.numeric(x)) {
    x_chr <- trimws(as.character(x))
    missing <- is.na(x) | x_chr == ""
    numeric_token <- grepl(
      "^[+-]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?$",
      x_chr
    )
    bad <- which(!missing & !numeric_token)
    if (length(bad) > 0L) {
      stop(sprintf(
        ".rhc_impute_continuous: variable '%s' contains %d non-numeric non-missing value(s); first bad value is '%s'.",
        var_name, length(bad), x_chr[bad[1L]]
      ), call. = FALSE)
    }
    x_num <- rep(NA_real_, length(x_chr))
    x_num[!missing] <- as.numeric(x_chr[!missing])
    x <- x_num
  }
  med <- if (is.null(training_rows)) stats::median(x, na.rm = TRUE) else
    stats::median(x[training_rows], na.rm = TRUE)
  if (is.na(med)) {
    stop(sprintf(
      ".rhc_impute_continuous: variable '%s' has no non-NA values (n=%d); cannot compute a median for imputation.",
      var_name, length(x)
    ), call. = FALSE)
  }
  x[is.na(x)] <- med
  x
}

.rhc_impute_mode <- function(x, var_name = "<unknown>", training_rows = NULL) {
  tab <- if (is.null(training_rows)) sort(table(x, useNA = "no"), decreasing = TRUE) else
    sort(table(x[training_rows], useNA = "no"), decreasing = TRUE)
  if (length(tab) == 0L) {
    stop(sprintf(
      ".rhc_impute_mode: variable '%s' has no non-NA values (n=%d); cannot determine a mode for imputation.",
      var_name, length(x)
    ), call. = FALSE)
  }
  mode_val <- names(tab)[1L]
  x[is.na(x) | x == ""] <- mode_val
  x
}


#' Build the RHC Analysis Cohort
#'
#' Cleans the raw RHC data, defines the treatment and outcome variables, and
#' selects the covariate set of Hirano & Imbens (2001). Returns a wide data
#' frame with one row per patient.
#'
#' @param raw Optional data frame from \code{\link{load_rhc_raw}}. If
#'   \code{NULL} (default), \code{load_rhc_raw()} is invoked.
#' @param outcome Character. \code{"death30"} (30-day all-cause mortality;
#'   binary, the default, aligning with FL-TTE and Connors 1996),
#'   \code{"death180"} (180-day mortality; binary), or \code{"los"}
#'   (hospital length of stay in days; continuous).
#' @param site_var Character. Raw RHC column to use for site partitioning.
#'   Default \code{"cat1"} (primary disease category) reproduces Connors
#'   1996's 9-category clinical split. Alternatives include \code{"ninsclas"}
#'   (insurance class; 6 levels, patients more homogeneous within a site),
#'   \code{"income"} (4 levels), \code{"race"} (3 levels), and \code{"ca"}
#'   (cancer status; 3 levels). The chosen column is automatically removed
#'   from the covariate set so that it is not also used as a predictor.
#' @param site_recode Optional named character vector. Remaps raw values of
#'   the \code{site_var} column before site partitioning. Names are raw
#'   levels, values are the desired output levels; setting a value to
#'   \code{NA} drops all patients with that raw level. Example: to merge
#'   \code{"Medicaid"} and \code{"Medicare & Medicaid"} and drop
#'   \code{"No insurance"} when \code{site_var = "ninsclas"}, pass
#'   \code{c("Medicaid" = "Medicaid_any", "Medicare & Medicaid" =
#'   "Medicaid_any", "No insurance" = NA)}.
#' @param covariate_profile \code{"historical"} preserves the original features.
#'   \code{"log_missing"} replaces urine output, creatinine, bilirubin and white
#'   cell count by \code{log1p} values and adds an observed urine-missingness
#'   indicator. \code{"grouped_log_missing"} additionally groups primary colon
#'   and lung cancer diagnoses and combines trauma/orthopedic diagnosis flags.
#'   These are named sensitivity specifications; no rows or outcomes select
#'   the coding. The grouped profile coarsens clinical information.
#' @return A data frame with columns: \code{A} (0/1 treatment, with 1 = RHC),
#'   \code{Y} (numeric outcome defined by \code{outcome}), \code{site_var}
#'   (the raw value of \code{site_var}, preserved for site partitioning),
#'   plus the covariates listed in \code{.RHC_CONTINUOUS_VARS},
#'   \code{.RHC_BINARY_VARS}, and the one-hot expansion of
#'   \code{.RHC_CATEGORICAL_VARS} (with \code{site_var} excluded).
#' @param preprocessing_rows Optional training-row indices in the retained
#'   cohort, after site and missing-treatment/outcome exclusions. NULL retains
#'   the historical whole-cohort imputation.
#' @param categorical_levels Optional named list of fixed category levels.
#' @export
build_rhc_cohort <- function(raw = NULL,
                             outcome = c("death30", "death180", "los"),
                             site_var = .RHC_SITE_VAR,
                             site_recode = NULL,
                             covariate_profile = c("historical", "log_missing", "grouped_log_missing"),
                             preprocessing_rows = NULL, categorical_levels = NULL) {
  outcome <- match.arg(outcome)
  covariate_profile <- match.arg(covariate_profile)
  if (is.null(raw)) raw <- load_rhc_raw()
  if (!is.character(site_var) || length(site_var) != 1L || !nzchar(site_var)) {
    stop("site_var must be a single non-empty character string.")
  }
  if (!site_var %in% colnames(raw)) {
    stop(sprintf("site_var = '%s' is not a column of the raw RHC data.", site_var))
  }
  if (!is.null(site_recode)) {
    if (!is.character(site_recode) || is.null(names(site_recode)) ||
        any(!nzchar(names(site_recode)))) {
      stop("site_recode must be a named character vector (NA values drop rows).")
    }
    raw_levels_seen <- setdiff(names(site_recode), unique(raw[[site_var]]))
    if (length(raw_levels_seen) > 0L) {
      warning(sprintf(
        "site_recode references level(s) not present in raw$%s: %s",
        site_var, paste(raw_levels_seen, collapse = ", ")
      ))
    }
  }

  required <- c(site_var, .RHC_TREATMENT_VAR,
                .RHC_CONTINUOUS_VARS, .RHC_BINARY_VARS, .RHC_CATEGORICAL_VARS)
  missing_cols <- setdiff(required, colnames(raw))
  if (length(missing_cols) > 0L) {
    stop("RHC raw data is missing expected columns: ",
         paste(missing_cols, collapse = ", "))
  }

  # Apply site_recode before any downstream processing: rows whose raw site
  # label maps to NA are filtered out; remaining rows have their site value
  # replaced by the mapped target. This step keeps the cohort consistent
  # across treatment/outcome/covariate construction below.
  if (!is.null(site_recode)) {
    raw_site <- as.character(raw[[site_var]])
    is_in_recode <- raw_site %in% names(site_recode)
    drop_levels <- names(site_recode)[is.na(site_recode)]
    keep_rows   <- !(raw_site %in% drop_levels)
    raw <- raw[keep_rows, , drop = FALSE]
    raw_site <- raw_site[keep_rows]
    remap_mask <- raw_site %in% setdiff(names(site_recode), drop_levels)
    raw_site[remap_mask] <- unname(site_recode[raw_site[remap_mask]])
    raw[[site_var]] <- raw_site
  }

  # ---- Treatment ----
  A_chr <- toupper(gsub("\\s+", "", raw[[.RHC_TREATMENT_VAR]]))
  A <- ifelse(A_chr %in% c("RHC", "RHCYES"), 1L,
       ifelse(A_chr %in% c("NORHC", "NO", "NORHCNO"), 0L, NA_integer_))

  # ---- Outcome ----
  Y <- switch(
    outcome,
    death30 = {
      if ("dth30" %in% colnames(raw)) {
        y_chr <- toupper(gsub("\\s+", "", as.character(raw$dth30)))
        ifelse(y_chr %in% c("YES", "1", "TRUE"), 1L,
        ifelse(y_chr %in% c("NO",  "0", "FALSE"), 0L, NA_integer_))
      } else {
        stop("Column dth30 not found; cannot build death30 outcome.")
      }
    },
    death180 = {
      if ("death" %in% colnames(raw)) {
        y_chr <- toupper(gsub("\\s+", "", as.character(raw$death)))
        ifelse(y_chr %in% c("YES", "1", "TRUE"), 1L,
        ifelse(y_chr %in% c("NO",  "0", "FALSE"), 0L, NA_integer_))
      } else {
        stop("Column death not found; cannot build death180 outcome.")
      }
    },
    los = {
      # Length of stay in days, computed from admission/discharge dates
      # stored as julian integers in Harrell's mirror.
      date_cols <- c("sadmdte", "dschdte")
      missing_date <- setdiff(date_cols, colnames(raw))
      if (length(missing_date) > 0L) {
        stop("Cannot compute length-of-stay outcome; missing columns: ",
             paste(missing_date, collapse = ", "))
      }
      as.numeric(raw$dschdte) - as.numeric(raw$sadmdte)
    }
  )

  # ---- Keep cases with both A and Y observed ----
  keep <- !is.na(A) & !is.na(Y)
  raw  <- raw[keep, , drop = FALSE]
  A    <- A[keep]
  Y    <- Y[keep]

  .validate_rhc_preprocessing_rows(preprocessing_rows, nrow(raw))

  # ---- Covariate subsets (drop the site variable so it is not a predictor) ----
  continuous_vars  <- setdiff(.RHC_CONTINUOUS_VARS,  site_var)
  binary_vars      <- setdiff(.RHC_BINARY_VARS,      site_var)
  categorical_vars <- setdiff(.RHC_CATEGORICAL_VARS, site_var)

  # ---- Covariates: continuous ----
  X_cont <- setNames(
    lapply(continuous_vars, function(v) .rhc_impute_continuous(raw[[v]], var_name = v,
      training_rows = preprocessing_rows)),
    continuous_vars
  )
  X_cont <- as.data.frame(X_cont, stringsAsFactors = FALSE)
  if (covariate_profile != "historical") {
    for (variable in intersect(c("urin1", "crea1", "bili1", "wblc1"), names(X_cont))) {
      if (any(X_cont[[variable]] < 0)) {
        stop("log1p sensitivity requires nonnegative values for ", variable, call. = FALSE)
      }
      X_cont[[variable]] <- log1p(X_cont[[variable]])
      names(X_cont)[names(X_cont) == variable] <- paste0("log1p_", variable)
    }
  }

  # ---- Covariates: binary indicators (coded 0/1 in source) ----
  X_bin <- setNames(
    lapply(binary_vars, function(v) {
      x <- raw[[v]]
      # Harrell's mirror stores these as 0/1 integers; be defensive anyway.
      if (is.numeric(x)) return(.rhc_impute_continuous(x, var_name = v,
        training_rows = preprocessing_rows))
      ch <- toupper(as.character(x))
      as.integer(ch %in% c("1", "YES", "TRUE"))
    }),
    binary_vars
  )
  X_bin <- as.data.frame(X_bin, stringsAsFactors = FALSE)
  if (covariate_profile != "historical" && "urin1" %in% continuous_vars) {
    X_bin$urin1_missing <- as.integer(is.na(raw$urin1) | trimws(as.character(raw$urin1)) == "")
  }
  if (covariate_profile == "grouped_log_missing") {
    if (!all(c("trauma", "ortho") %in% names(X_bin))) {
      stop("Grouped RHC coding requires both trauma and ortho as covariates.", call. = FALSE)
    }
    X_bin$injury_diagnosis <- as.integer(X_bin$trauma == 1 | X_bin$ortho == 1)
    X_bin[c("trauma", "ortho")] <- NULL
  }

  # ---- Covariates: categorical → one-hot ----
  if (length(categorical_vars) > 0L) {
    cat_frame <- setNames(
      lapply(categorical_vars, function(v) {
        x <- .rhc_impute_mode(as.character(raw[[v]]), var_name = v,
          training_rows = preprocessing_rows)
        if (covariate_profile == "grouped_log_missing" && v == "cat1") {
          x[x %in% c("Colon Cancer", "Lung Cancer")] <- "Solid Cancer"
        }
        if (is.null(categorical_levels)) return(factor(x))
        levels <- categorical_levels[[v]]
        if (is.null(levels) || anyNA(x) || any(!x %in% levels)) {
          stop("Unknown or missing value in the fixed RHC category dictionary: ", v,
               call. = FALSE)
        }
        factor(x, levels = levels)
      }),
      categorical_vars
    )
    cat_frame <- as.data.frame(cat_frame, stringsAsFactors = TRUE)
    # Build a full-rank treatment-coded design.  Using `~ . - 1` here is not
    # equivalent: model.matrix() keeps all levels of the first factor when the
    # formula has no intercept.  The nuisance solvers add their own intercept,
    # so that construction made (for example) sexFemale + sexMale exactly
    # collinear with the fitted intercept.  Build the ordinary intercept model
    # and then remove only its intercept column instead.
    cat_mm <- stats::model.matrix(~ ., data = cat_frame)
    cat_mm <- cat_mm[, colnames(cat_mm) != "(Intercept)", drop = FALSE]
    X_cat <- as.data.frame(cat_mm, stringsAsFactors = FALSE, check.names = TRUE)
  } else {
    X_cat <- data.frame(row.names = seq_along(A))
  }

  cohort <- data.frame(
    A        = A,
    Y        = Y,
    site_var = as.character(raw[[site_var]]),
    X_cont,
    X_bin,
    X_cat,
    stringsAsFactors = FALSE,
    check.names      = TRUE
  )

  attr(cohort, "outcome")         <- outcome
  attr(cohort, "site_var_col")    <- site_var
  attr(cohort, "continuous_vars") <- colnames(X_cont)
  attr(cohort, "binary_vars")     <- colnames(X_bin)
  attr(cohort, "onehot_vars")     <- colnames(X_cat)
  cohort
}


# Consolidate the 9 cat1 categories into the top K by count, merging the
# remainder into an "Other" bucket. Returns a factor with K levels.
.rhc_site_factor <- function(site_var, K) {
  if (K < 2L) stop("K must be >= 2; RoCE requires at least 1 target and 1 source site.")
  counts <- sort(table(site_var), decreasing = TRUE)
  top    <- names(counts)[seq_len(min(K, length(counts)))]
  if (length(top) < K) {
    stop(sprintf("site_var has only %d distinct categories; cannot build K = %d sites.",
                 length(top), K))
  }
  # If more than K categories, fold the tail into the smallest retained site.
  # Record the raw composition explicitly so a retained display label can
  # never be mistaken for an unmodified raw stratum downstream.
  collapsed <- ifelse(site_var %in% top, site_var, top[length(top)])
  site_factor <- factor(collapsed, levels = top)
  raw_to_retained <- stats::setNames(
    ifelse(names(counts) %in% top, names(counts), top[length(top)]),
    names(counts)
  )
  site_components <- stats::setNames(
    vapply(
      top,
      function(label) {
        paste(names(raw_to_retained)[raw_to_retained == label], collapse = ";")
      },
      character(1L)
    ),
    top
  )
  attr(site_factor, "site_components") <- site_components
  site_factor
}


#' Assemble an RHC \code{data_split} for RoCE
#'
#' Partitions the cohort into a target site and \code{K - 1} source sites
#' using the primary disease category (\code{cat1}) as the site variable.
#' The largest category is used as the target unless overridden.
#'
#' @param cohort Output of \code{\link{build_rhc_cohort}}. If \code{NULL},
#'   the cohort is built with default options.
#' @param K Integer. Total number of sites including the target
#'   (\code{K} is the number of \emph{source} sites plus 1). Defaults to 5.
#' @param target_site Optional character naming the cat1 level to use as the
#'   target. If \code{NULL}, the largest retained category is selected.
#' @param phi Optional basis-expansion function applied to the numeric
#'   covariate matrix before constructing \code{W_outcome} / \code{Z_site}.
#'   Must accept a numeric matrix and return a numeric matrix with at least
#'   as many columns. Defaults to \code{\link[base]{identity}} so that no
#'   dimension inflation occurs; a reference B-spline implementation is
#'   available in \code{\link{phi_rhc_bspline}} for opt-in use.
#' @param standardize_features Logical. If \code{TRUE} (default), center and
#'   scale every working-feature column after applying \code{phi}.  The pooled
#'   first and second moments used here are federated-computable summaries; raw
#'   covariates retained in \code{X} and \code{X_dagger} are not modified.
#' @param seed Integer. Random seed passed to R's RNG for any tie-breaking
#'   steps; retained for reproducibility (currently unused). Default 42.
#' @param preprocessing_rows Optional retained-cohort training rows used to
#'   compute feature centers and scales. The same transform is applied to all rows.
#' @return A \code{data_split} list with elements named \code{"t"} (target)
#'   and \code{"s1"}, ..., \code{"s<K-1>"} (source sites). Each element is a
#'   list with fields \code{n}, \code{X}, \code{X_dagger}, \code{A},
#'   \code{Y}, \code{Z_site}, and \code{W_outcome}, matching the output of
#'   \code{\link{split_data_by_site}} used elsewhere in RoCE.
#' @seealso \code{\link{run_crossfit}}, \code{\link{split_data_by_site}},
#'   \code{\link{phi_rhc_bspline}}.
#' @export
build_rhc_data_split <- function(cohort      = NULL,
                                 K           = 5L,
                                 target_site = NULL,
                                 phi         = base::identity,
                                 standardize_features = TRUE,
                                 seed        = 42L,
                                 preprocessing_rows = NULL) {
  if (is.null(cohort)) cohort <- build_rhc_cohort()
  if (!is.function(phi)) {
    stop("phi must be a function mapping a numeric matrix to a numeric matrix.")
  }
  if (length(standardize_features) != 1L || is.na(standardize_features) ||
      !is.logical(standardize_features)) {
    stop("standardize_features must be TRUE or FALSE.")
  }
  K <- as.integer(K)
  set.seed(seed)

  site_factor <- .rhc_site_factor(cohort$site_var, K = K)
  retained <- levels(site_factor)
  retained_components <- attr(site_factor, "site_components")

  if (is.null(target_site)) {
    target_site <- retained[1L]  # largest
  } else if (!target_site %in% retained) {
    stop(sprintf("target_site = '%s' is not among the retained sites: %s",
                 target_site, paste(retained, collapse = ", ")))
  }

  covariate_cols <- setdiff(colnames(cohort), c("A", "Y", "site_var"))
  X_full <- as.matrix(cohort[, covariate_cols, drop = FALSE])
  storage.mode(X_full) <- "double"

  X_phi_raw <- phi(X_full)
  if (!is.matrix(X_phi_raw) || nrow(X_phi_raw) != nrow(X_full)) {
    stop("phi must return a numeric matrix with the same number of rows as its input.")
  }
  storage.mode(X_phi_raw) <- "double"
  if (any(!is.finite(X_phi_raw))) {
    stop("phi must return only finite numeric working features.")
  }

  .validate_rhc_preprocessing_rows(preprocessing_rows, nrow(X_phi_raw))

  feature_center <- if (isTRUE(standardize_features)) {
    if (is.null(preprocessing_rows)) colMeans(X_phi_raw) else
      colMeans(X_phi_raw[preprocessing_rows, , drop = FALSE])
  } else {
    rep(0, ncol(X_phi_raw))
  }
  centered_features <- sweep(X_phi_raw, 2L, feature_center, FUN = "-")
  feature_scale <- if (isTRUE(standardize_features)) {
    if (is.null(preprocessing_rows)) sqrt(colMeans(centered_features^2)) else
      sqrt(colMeans(centered_features[preprocessing_rows, , drop = FALSE]^2))
  } else {
    rep(1, ncol(X_phi_raw))
  }
  # A constant feature carries no scale information.  Keeping it unchanged
  # preserves a user-supplied phi() interface while avoiding division by zero;
  # the default RHC design is full rank and has no constant working columns.
  invalid_scale <- !is.finite(feature_scale) |
    feature_scale <= sqrt(.Machine$double.eps)
  feature_scale[invalid_scale] <- 1
  X_phi <- sweep(centered_features, 2L, feature_scale, FUN = "/")
  colnames(X_phi) <- colnames(X_phi_raw)

  # Map target/source labels to the RoCE convention ("t", "s1", ...).
  source_sites <- setdiff(retained, target_site)
  site_labels  <- character(length(site_factor))
  site_labels[site_factor == target_site] <- "t"
  for (j in seq_along(source_sites)) {
    site_labels[site_factor == source_sites[j]] <- paste0("s", j)
  }

  build_one <- function(lbl) {
    idx <- which(site_labels == lbl)
    list(
      n         = length(idx),
      X         = X_full[idx, , drop = FALSE],
      X_dagger  = X_full[idx, , drop = FALSE],
      A         = as.integer(cohort$A[idx]),
      Y         = cohort$Y[idx],
      Z_site    = X_phi[idx, , drop = FALSE],
      W_outcome = X_phi[idx, , drop = FALSE]
    )
  }

  data_split <- list()
  data_split[["t"]] <- build_one("t")
  for (j in seq_along(source_sites)) {
    lbl <- paste0("s", j)
    data_split[[lbl]] <- build_one(lbl)
  }

  attr(data_split, "outcome")         <- attr(cohort, "outcome")
  attr(data_split, "site_mapping")    <- c(target = target_site,
                                           stats::setNames(source_sites,
                                                           paste0("s", seq_along(source_sites))))
  attr(data_split, "site_components") <- c(
    target = unname(retained_components[target_site]),
    stats::setNames(
      unname(retained_components[source_sites]),
      paste0("s", seq_along(source_sites))
    )
  )
  attr(data_split, "phi_applied")     <- !identical(phi, base::identity)
  attr(data_split, "features_standardized") <- isTRUE(standardize_features)
  attr(data_split, "feature_center")  <- feature_center
  attr(data_split, "feature_scale")   <- feature_scale
  attr(data_split, "covariate_cols")  <- covariate_cols
  data_split
}

.rhc_method_row <- function(method_label, estimate, se) {
  if (is.null(estimate) || !is.finite(estimate) || !is.finite(se)) {
    return(data.frame(
      method = method_label,
      estimate = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    method = method_label,
    estimate = estimate,
    se = se,
    ci_lower = estimate - Z_ALPHA_05 * se,
    ci_upper = estimate + Z_ALPHA_05 * se,
    stringsAsFactors = FALSE
  )
}


#' Run the Full RHC Real-Data Experiment
#'
#' Builds the RHC \code{data_split}, runs arm-specific RoCE (one-round
#' cross-fitting)
#' together with the comparison methods implemented in
#' \code{\link{run_all_comparisons}}, and returns a uniform result list
#' ready for downstream tabulation / plotting. This is the single entry
#' point invoked by \code{realdata.R} and is exported so that interactive
#' users can reproduce the paper's Section 6 figures without dealing with
#' SLURM.
#'
#' @param K Integer. Number of sites (target + sources). Default 5.
#' @param target_site Character. See \code{\link{build_rhc_data_split}};
#'   defaults to the largest \code{cat1} category (ARF).
#' @param outcome Character. Passed to \code{\link{build_rhc_cohort}};
#'   one of \code{"death30"} (default), \code{"death180"}, or \code{"los"}.
#' @param site_var Character. Raw RHC column to use for site partitioning;
#'   passed to \code{\link{build_rhc_cohort}} and
#'   \code{\link{build_rhc_data_split}}. Defaults to
#'   \code{"cat1"} (primary disease category). Use \code{"ninsclas"} for a
#'   more homogeneous insurance-level split that typically yields tighter
#'   covariate overlap across sites.
#' @param site_recode Optional named character vector passed to
#'   \code{\link{build_rhc_cohort}} to remap / drop raw site levels before
#'   partitioning (see that function's documentation).
#' @param family GLM family. Inferred from \code{outcome} when \code{NULL}
#'   (default): \code{death30}/\code{death180} \eqn{\Rightarrow}
#'   \code{"binomial"}, \code{los} \eqn{\Rightarrow} \code{"gaussian"}.
#' @param n_folds Integer. Number of outer cross-fitting folds; default
#'   \code{N_FOLDS_DEFAULT}.
#' @param A_val Integer (0/1). Treatment value for the potential-outcome
#'   contrast. Default 1 (RHC exposure).
#' @param phi Basis-expansion function applied to the covariate matrix
#'   inside \code{\link{build_rhc_data_split}}. Default \code{identity}.
#' @param seed Integer. Passed to \code{build_rhc_data_split}.
#' @param verbose Logical. Propagated to \code{\link{run_crossfit}}.
#' @param n_cores Number of cores for parallel processing of source sites
#'   inside \code{\link{run_crossfit}}. \code{NULL} or \code{1L} keeps the
#'   default sequential path; \code{-1L} uses all but one available core.
#'   Passing \code{-1L} is recommended on MSI compute nodes to avoid the
#'   hours-long wall time of single-threaded cross-fitting on K source
#'   sites.
#' @param M_tau_inference Numeric. Truncation cap on the linear predictor
#'   \eqn{\phi(\boldsymbol{X})^\top\boldsymbol{\gamma}} inside the
#'   inference-time IF, which in turn bounds the calibration weight
#'   \eqn{\exp(-\phi^\top\boldsymbol{\gamma})} to \eqn{[e^{-M}, e^M]}. The
#'   default \code{M_TAU_INFERENCE_RHC = 3} bounds the weight by
#'   \eqn{\approx 20}, which comfortably
#'   covers realistic RHC density ratios while preventing the catastrophic
#'   pairwise estimates (e.g. \eqn{\widehat{\mu}^a_{t,s}} at order
#'   \eqn{10^2} or beyond) that arise in small sources with limited
#'   covariate overlap. The TATE manuscript analysis explicitly uses
#'   \code{M_TAU_DEFAULT = 5} for both fitting and inference.
#' @return A list with elements:
#' \describe{
#'   \item{methods}{Data frame with one row per method, columns
#'     \code{method}, \code{estimate}, \code{se}, \code{ci_lower},
#'     \code{ci_upper}. Methods and labels are unified with the simulation
#'     figure: \code{"Target-only"}, \code{"SS"}, \code{"IVW"},
#'     \code{"Federated-DR"}, \code{"Pooled-DR"}, and the one-round
#'     \code{"RoCE"} estimator.}
#'   \item{pairwise}{Data frame of one-round source-assisted pairwise
#'     estimates with columns \code{source}, \code{estimate}, \code{se}.}
#'   \item{weights}{Named numeric vector of RoCE aggregation weights
#'     \eqn{\widehat{\boldsymbol{\eta}}} (length K-1; names match the
#'     \code{s1, s2, ...} source labels of \code{pairwise}).}
#'   \item{target_only}{Target-only cross-fitted point estimate.}
#'   \item{metadata}{Named list with run configuration
#'     (\code{K}, \code{target_site}, \code{outcome}, \code{family},
#'     \code{n_folds}, site sample sizes, and RoCE pairwise gap
#'     penalties \eqn{d_j^2}).}
#' }
#' @seealso \code{\link{build_rhc_data_split}},
#'   \code{\link{run_crossfit}}, \code{\link{run_all_comparisons}},
#'   \code{\link{plot_forest_methods}}.
#' @export
run_rhc_experiment <- function(K               = 5L,
                               target_site     = NULL,
                               outcome         = c("death30", "death180", "los"),
                               family          = NULL,
                               n_folds         = N_FOLDS_DEFAULT,
                               A_val           = 1L,
                               phi             = base::identity,
                               seed            = 42L,
                               verbose         = TRUE,
                               n_cores         = NULL,
                               M_tau_inference = M_TAU_INFERENCE_RHC,
                               site_var        = .RHC_SITE_VAR,
                               site_recode     = NULL) {
  outcome <- match.arg(outcome)
  if (is.null(family)) {
    family <- if (outcome %in% c("death30", "death180")) "binomial" else "gaussian"
  }
  family  <- match.arg(family, choices = c("binomial", "gaussian"))
  K       <- as.integer(K)
  n_folds <- as.integer(n_folds)
  A_val   <- as.integer(A_val)

  cohort <- build_rhc_cohort(outcome = outcome, site_var = site_var,
                             site_recode = site_recode)
  data_split <- build_rhc_data_split(cohort = cohort, K = K,
                                     target_site = target_site,
                                     phi = phi, seed = seed)

  if (verbose) {
    cat(sprintf(
      "[run_rhc_experiment] outcome=%s family=%s K=%d n_folds=%d site_var=%s target=%s\n",
      outcome, family, K, n_folds, site_var,
      attr(data_split, "site_mapping")["target"]
    ))
  }

  # ---- RoCE: two-round and one-round cross-fitted variants -------------
  # Both communication modes share the same penalized aggregation but differ
  # in how source-site initial nuisance estimators are obtained (two-round:
  # source-specific; one-round: target-only plug-in). The main manuscript
  # treats them as separate methods, so both are reported in the results
  # table. Source-site parallelism (n_cores) is opt-in: the default NULL/1
  # preserves the reproducible sequential path, while -1 on compute nodes
  # unlocks the K-1 parallel source-site fits that otherwise serialise the
  # cross-fitting loop.
  roce_one_round <- run_crossfit(
    data_split         = data_split,
    n_folds            = n_folds,
    communication_mode = "one_round",
    family             = family,
    A_val              = A_val,
    n_cores            = n_cores,
    M_tau_inference    = M_tau_inference,
    verbose            = verbose
  )
  # The one-round variant is the reported RoCE estimator (matching the
  # simulation figure); it also supplies the pairwise source-assisted estimates,
  # aggregation weights, and target-only reference used below.
  roce_res <- roce_one_round

  # ---- Baselines --------------------------------------------------------
  comparison_res <- run_all_comparisons(
    data_split   = data_split,
    family       = family,
    A_val        = A_val,
    use_crossfit = TRUE,
    n_folds      = n_folds,
    include_tilted = FALSE,
    n_cores      = n_cores
  )

  # Assemble uniform "methods" data frame expected by plot_forest_methods().
  pull <- function(method_name, label = method_name) {
    res <- comparison_res[[method_name]]
    if (is.null(res)) return(.rhc_method_row(label, NA_real_, NA_real_))
    .rhc_method_row(
      label,
      estimate = res$estimate,
      se = if (!is.null(res$se)) res$se else sqrt(res$variance)
    )
  }

  # Method set and display labels are unified with the simulation figure: the
  # one-round RoCE plus the penalized baselines, with publication labels. The
  # unregularized tilted-AIPW baseline and the two-round variant are omitted so
  # the comparison forest plot reports a single, consistent set of estimators.
  methods_df <- rbind(
    pull("target_only",      "Target-only"),
    pull("sample_size",      "SS"),
    pull("inverse_variance", "IVW"),
    pull("federated_dr",     "Federated-DR"),
    pull("pooled_dr",        "Pooled-DR"),
    .rhc_method_row(
      "RoCE",
      estimate = roce_one_round$estimate,
      se = roce_one_round$se
    )
  )
  methods_df$method <- factor(
    methods_df$method,
    levels = c("Target-only", "SS", "IVW", "Federated-DR", "Pooled-DR", "RoCE"))

  # ---- Pairwise source-assisted estimates ------------------------------
  # Pairwise estimates come from the reported one-round RoCE run.
  source_labels <- names(roce_res$source_estimates)
  pairwise_df <- data.frame(
    source   = source_labels,
    estimate = as.numeric(roce_res$source_estimates),
    se       = NA_real_,
    stringsAsFactors = FALSE
  )

  # Attach source labels to the aggregation-weight vector so downstream
  # callers (realdata.R, plot_aggregation_weights) can rely on names()
  # without hard-coding the label convention.
  weights_named <- stats::setNames(as.numeric(roce_res$weights), source_labels)

  # Squared target-source gaps (penalty scale) for optional plot annotation.
  d_sq <- (roce_res$target_only$estimate - pairwise_df$estimate)^2

  metadata <- list(
    K              = length(source_labels),
    n_sites        = length(data_split),
    site_var       = site_var,
    target_site    = attr(data_split, "site_mapping")["target"],
    source_sites   = attr(data_split, "site_mapping")[paste0("s", seq_along(source_labels))],
    site_components = attr(data_split, "site_components"),
    outcome        = outcome,
    family         = family,
    n_folds        = n_folds,
    A_val          = A_val,
    site_n         = vapply(data_split, function(s) s$n, integer(1L)),
    pairwise_d_sq  = d_sq
  )

  list(
    methods     = methods_df,
    pairwise    = pairwise_df,
    weights     = weights_named,
    target_only = roce_res$target_only$estimate,
    metadata    = metadata
  )
}

# Build only controls supported by the selected target nuisance program.
.rhc_calibration_control <- function(target_program, source_program, radius, source_lambda_rules = NULL) {
  target_program <- match.arg(target_program, c("hou_calibrated", "lasso"))
  source_program <- match.arg(source_program, c("calibrated", "standard"))
  control <- list(recipe = "score_derivative", target_propensity_initialization = "logistic",
                  source_nuisance_method = source_program)
  if (target_program == "hou_calibrated") control$target_radius <- radius
  if (!is.null(source_lambda_rules)) control$source_lambda_rules <- source_lambda_rules
  .validate_calibration_control(control, target_program)
  control
}

#' RoCE TATE analysis of the RHC data
#'
#' Runs the treated and control nuisance fits on a shared fold partition and
#' applies the selected TATE aggregation rule. The historical common-weight
#' default is retained; joint TATE uses two arm-specific weight vectors.
#' The existing arm-specific \code{run_rhc_experiment()} interface is
#' unchanged and remains available for secondary potential-outcome summaries.
#'
#' @inheritParams run_rhc_experiment
#' @param K Integer. Total number of retained insurance strata, including the
#'   target. The support-screened manuscript driver passes \code{K = 4}, which
#'   produces three source sites; the returned metadata therefore records
#'   \code{K = 3} and \code{n_sites = 4}. The function default remains five
#'   strata so earlier all-site sensitivity analyses stay reproducible.
#' @param aggregation_lambda Positive truncated-Wald multiplier. Its reciprocal
#'   is the source penalty-activation cutoff. The locked manuscript analysis
#'   uses \code{0.5} (cutoff \code{2}); historical defaults are unchanged.
#' @param comparison_methods Character vector of TATE comparison methods.
#'   Defaults to the four baselines reported in the manuscript. Supply
#'   \code{character(0)} to run only target-only and RoCE.
#' @param variance_method Comparison-method variance estimator, passed to
#'   \code{\link{run_all_comparisons_tate}}.
#' @param n_bootstrap Number of multiplier-bootstrap draws for comparison-method
#'   standard errors when \code{variance_method = "bootstrap"}.
#' @param parallel_arms Logical. Fit the two RoCE treatment-arm nuisance
#'   pipelines concurrently. When enabled, \code{n_cores} is interpreted per
#'   arm and the caller should allocate approximately \code{2 * n_cores} CPUs.
#' @param nlambda_init Number of candidate penalties in each initial nuisance
#'   path. The manuscript uses \code{LAMBDA_GRID_SIZE_STANDARD = 100}.
#' @param nuisance_lambda_rule Nuisance-model CV selection rule passed to
#'   \code{\link{run_tate_crossfit}}. The primary analysis uses \code{"min"};
#'   \code{"1se"} is available as a more strongly regularized stability
#'   sensitivity analysis.
#' @param M_tau Numeric fitting-stage truncation radius. The manuscript uses
#'   a source radius of \code{3} and a separately specified target radius of \code{5}.
#' @param aggregation_mode,target_nuisance_method,source_validation_method,crossfit_layers,calibration_control,calibration_layout,nuisance_solver,nuisance_tol,nuisance_cv_certificate,checkpoint_dir
#'   Explicit fitting controls forwarded to \code{\link{run_tate_crossfit}}.
#'   Defaults preserve the historical RHC analysis. Current calibrated analyses
#'   should select their layer, anchor, calibration recipe and radii explicitly.
#' @param fold_seed Optional positive integer for RoCE's outer fold assignment.
#'   \code{NULL} preserves the sample-size-based historical assignment.
#'   The ordinary target reference uses these same outer folds; other baselines
#'   retain their own documented fitting procedures on the identical cohort.
#' @param covariate_profile Named covariate sensitivity passed to
#'   \code{\link{build_rhc_cohort}}. The historical default is unchanged.
#' @param preprocessing \code{"cohort"} preserves historical preprocessing.
#'   \code{"outer_fold"} learns the shared transform outside each RoCE outer
#'   evaluation fold. Initial/calibration CV still reuse that outer-training
#'   transform. SS/IVW exclude their own site-validation rows when learning
#'   preprocessing; the two DR baselines retain their full-sample fitting protocol.
#' @return A list containing the TATE method comparison, pairwise
#'   source-assisted TATE estimates and arm-labeled Wald diagnostics,
#'   aggregation weights, ordinary and calibrated target references,
#'   the full TATE fit, and run metadata.
#' @export
run_rhc_tate_experiment <- function(
    K = 5L,
    target_site = NULL,
    outcome = c("death30", "death180", "los"),
    family = NULL,
    n_folds = N_FOLDS_DEFAULT,
    phi = base::identity,
    seed = 42L,
    verbose = TRUE,
    n_cores = NULL,
    M_tau_inference = M_TAU_INFERENCE_DEFAULT,
    M_tau = M_TAU_DEFAULT,
    site_var = .RHC_SITE_VAR,
    site_recode = NULL,
    aggregation_lambda = AGG_WALD_LAMBDA,
    comparison_methods = c(
      "sample_size", "inverse_variance", "federated_dr", "pooled_dr"
    ),
    variance_method = c("bootstrap", "analytic"),
    n_bootstrap = BOOTSTRAP_REPLICATES_DEFAULT,
    parallel_arms = FALSE,
    nlambda_init = LAMBDA_GRID_SIZE_STANDARD,
    nuisance_lambda_rule = c("min", "1se"),
    aggregation_mode = c("common_tate", "separate_arms", "joint_tate"),
    target_nuisance_method = c("lasso", "hou_calibrated"),
    source_validation_method = c("initial", "calibrated"),
    crossfit_layers = NULL,
    calibration_control = NULL,
    calibration_layout = c("block", "compact"),
    nuisance_solver = NULL,
    nuisance_tol = TOL_DEFAULT,
    nuisance_cv_certificate = NULL,
    checkpoint_dir = NULL,
    fold_seed = NULL,
    covariate_profile = c("historical", "log_missing", "grouped_log_missing"),
    preprocessing = c("cohort", "outer_fold")) {
  outcome <- match.arg(outcome)
  covariate_profile <- match.arg(covariate_profile)
  preprocessing <- match.arg(preprocessing)
  variance_method <- match.arg(variance_method)
  aggregation_mode <- match.arg(aggregation_mode)
  target_nuisance_method <- match.arg(target_nuisance_method)
  calibration_layout <- match.arg(calibration_layout)
  if (!is.null(fold_seed) && (length(fold_seed) != 1L ||
      !is.numeric(fold_seed) || !is.finite(fold_seed) || fold_seed < 1 ||
      fold_seed > .Machine$integer.max - K || fold_seed != floor(fold_seed))) {
    stop("fold_seed must be NULL or one positive integer.", call. = FALSE)
  }
  n_bootstrap <- .validate_bootstrap_replicates(
    n_bootstrap, "run_rhc_tate_experiment"
  )
  nuisance_lambda_rule <- .match_nuisance_lambda_rule(
    nuisance_lambda_rule, "run_rhc_tate_experiment"
  )
  if (is.null(family)) {
    family <- if (outcome %in% c("death30", "death180")) "binomial" else "gaussian"
  }
  family <- match.arg(family, choices = c("binomial", "gaussian"))
  K <- as.integer(K)
  n_folds <- as.integer(n_folds)

  cohort <- build_rhc_cohort(
    outcome = outcome, site_var = site_var, site_recode = site_recode,
    covariate_profile = covariate_profile
  )
  data_split <- build_rhc_data_split(
    cohort = cohort, K = K, target_site = target_site, phi = phi, seed = seed
  )
  folds <- build_crossfit_folds(data_split, n_folds)
  if (!is.null(fold_seed)) {
    folds$target_folds <- partition_into_folds(data_split$t, n_folds, seed = fold_seed)
    source_sites <- setdiff(names(data_split), "t")
    folds$source_folds <- stats::setNames(lapply(seq_along(source_sites), function(j) {
      partition_into_folds(data_split[[source_sites[j]]], n_folds, seed = fold_seed + j)
    }), source_sites)
  }
  propensity_cache <- new.env(hash = TRUE, parent = emptyenv())

  comparison_data <- data_split
  if (preprocessing == "outer_fold") {
    prepared <- .prepare_rhc_outer_preprocessing(cohort, data_split, folds,
      cohort_arguments = list(outcome = outcome, site_var = site_var,
        site_recode = site_recode, covariate_profile = covariate_profile),
      split_arguments = list(K = K, target_site = target_site, phi = phi, seed = seed))
    folds <- prepared$folds
    comparison_data <- prepared$baseline_data
  }

  tate_fit <- run_tate_crossfit(
    data_split = data_split,
    n_folds = n_folds,
    communication_mode = "one_round",
    family = family,
    n_cores = n_cores,
    M_tau = M_tau,
    M_tau_inference = M_tau_inference,
    lambda_selection = aggregation_lambda,
    nlambda_init = nlambda_init,
    nuisance_lambda_rule = nuisance_lambda_rule,
    precomputed_folds = folds,
    target_only_ps_cache = propensity_cache,
    parallel_arms = parallel_arms,
    aggregation_mode = aggregation_mode,
    target_nuisance_method = target_nuisance_method,
    source_validation_method = source_validation_method,
    crossfit_layers = crossfit_layers,
    calibration_control = calibration_control,
    calibration_layout = calibration_layout,
    nuisance_solver = nuisance_solver,
    nuisance_tol = nuisance_tol,
    nuisance_cv_certificate = nuisance_cv_certificate,
    checkpoint_dir = checkpoint_dir,
    verbose = verbose
  )

  ordinary_target <- if (target_nuisance_method == "hou_calibrated") {
    .target_tate_reference(folds$target_folds, family, nuisance_lambda_rule)
  } else {
    tate_fit$target_only
  }
  target_anchor <- tate_fit$target_only

  comparison_labels <- c(
    sample_size = "SS",
    inverse_variance = "IVW",
    federated_dr = "Federated-DR",
    pooled_dr = "Pooled-DR"
  )
  comparison_methods <- unique(as.character(comparison_methods))
  unknown_comparisons <- setdiff(
    comparison_methods, names(comparison_labels)
  )
  if (length(unknown_comparisons) > 0L) {
    stop(sprintf(
      "run_rhc_tate_experiment: unknown comparison method(s): %s.",
      paste(unknown_comparisons, collapse = ", ")
    ), call. = FALSE)
  }

  comparison_results <- if (length(comparison_methods) > 0L) {
    run_all_comparisons_tate(
      data_split = comparison_data,
      family = family,
      n_folds = n_folds,
      variance_method = variance_method,
      n_bootstrap = n_bootstrap,
      include_tilted = FALSE,
      methods = comparison_methods,
      n_cores = n_cores
    )
  } else {
    list()
  }
  comparison_rows <- lapply(comparison_methods, function(method_name) {
    method_result <- comparison_results[[method_name]]
    .rhc_method_row(
      comparison_labels[[method_name]],
      method_result$estimate,
      method_result$se
    )
  })
  methods <- do.call(
    rbind,
    c(
      list(.rhc_method_row(
        "Target-only", ordinary_target$estimate, ordinary_target$se
      )),
      if (target_nuisance_method == "hou_calibrated") list(.rhc_method_row(
        "Calibrated target-only", target_anchor$estimate, target_anchor$se
      )),
      comparison_rows,
      list(.rhc_method_row("RoCE", tate_fit$estimate, tate_fit$se))
    )
  )
  method_levels <- c(
    "Target-only",
    if (target_nuisance_method == "hou_calibrated") "Calibrated target-only",
    unname(comparison_labels[comparison_methods]),
    "RoCE"
  )
  methods$method <- factor(methods$method, levels = method_levels)

  source_labels <- names(tate_fit$source_estimates)
  fold_wald <- tate_fit$fold_wald_statistics
  fold_penalty <- tate_fit$fold_penalty_coefficients
  arm_specific <- aggregation_mode != "common_tate"
  weights_by_arm <- if (arm_specific) tate_fit$weights_by_arm else
    list(mu1 = tate_fit$weights, mu0 = tate_fit$weights)
  pairwise <- data.frame(source = source_labels,
    estimate = as.numeric(tate_fit$source_estimates), stringsAsFactors = FALSE)
  for (arm in if (arm_specific) c("mu1", "mu0") else "common") {
    suffix <- if (arm_specific) paste0("_", arm) else ""
    columns <- if (arm_specific) paste0(arm, ":", source_labels) else source_labels
    weights <- if (arm_specific) weights_by_arm[[arm]] else tate_fit$weights
    stopifnot(length(weights) == length(source_labels),
              all(columns %in% colnames(fold_wald)),
              all(columns %in% colnames(fold_penalty)))
    pairwise[[paste0("weight", suffix)]] <- as.numeric(weights[source_labels])
    pairwise[[paste0("mean_wald_statistic", suffix)]] <- colMeans(fold_wald[, columns, drop = FALSE])
    pairwise[[paste0("max_wald_statistic", suffix)]] <- apply(fold_wald[, columns, drop = FALSE], 2L, max)
    pairwise[[paste0("penalty_activation_fraction", suffix)]] <- colMeans(fold_penalty[, columns, drop = FALSE] > 0)
  }
  pairwise_gap <- pairwise$estimate - target_anchor$estimate
  recode_specification <- if (is.null(site_recode)) {
    "none"
  } else {
    paste(
      paste0(
        names(site_recode), "->",
        ifelse(is.na(site_recode), "<excluded>", unname(site_recode))
      ),
      collapse = ";"
    )
  }

  list(
    methods = methods,
    pairwise = pairwise,
    weights = tate_fit$weights,
    weights_by_arm = weights_by_arm,
    target_only = ordinary_target,
    target_anchor = target_anchor,
    tate_fit = tate_fit,
    comparison_results = comparison_results,
    metadata = list(
      estimand = "TATE",
      preprocessing = preprocessing,
      data_sha256 = digest::digest(data_split, algo = "sha256"),
      covariate_profile = covariate_profile,
      fold_sha256 = digest::digest(folds, algo = "sha256"),
      K = length(source_labels),
      n_sites = length(data_split),
      site_var = site_var,
      site_recode_specification = recode_specification,
      excluded_site_levels = if (is.null(site_recode)) {
        "none"
      } else {
        paste(names(site_recode)[is.na(site_recode)], collapse = ";")
      },
      target_site = attr(data_split, "site_mapping")["target"],
      source_sites = attr(data_split, "site_mapping")[
        paste0("s", seq_along(source_labels))
      ],
      site_components = attr(data_split, "site_components"),
      outcome = outcome,
      family = family,
      n_folds = n_folds,
      seed = seed,
      fold_seed = fold_seed,
      aggregation_mode = aggregation_mode,
      crossfit_layers = tate_fit$crossfit_levels,
      target_nuisance_method = target_nuisance_method,
      source_validation_method = tate_fit$source_validation_method,
      calibration_control = tate_fit$calibration_control,
      nuisance_solver = tate_fit$nuisance_solver,
      nuisance_tol = nuisance_tol,
      nuisance_cv_certificate = tate_fit$nuisance_cv_certificate,
      nlambda_init = nlambda_init,
      nuisance_lambda_rule = nuisance_lambda_rule,
      M_tau = M_tau,
      M_tau_inference = M_tau_inference,
      aggregation_lambda = aggregation_lambda,
      aggregation_cutoff = 1 / aggregation_lambda,
      target_anchor_weight = if (!arm_specific) 1 - sum(tate_fit$weights) else NULL,
      target_anchor_weight_by_arm = vapply(weights_by_arm, function(weight) 1 - sum(weight), numeric(1L)),
      comparison_methods = comparison_methods,
      variance_method = variance_method,
      n_bootstrap = n_bootstrap,
      parallel_arms = parallel_arms,
      source_cores_per_arm = n_cores,
      site_n = vapply(data_split, function(site) site$n, integer(1L)),
      pairwise_gap = pairwise_gap,
      pairwise_gap_sq = pairwise_gap^2
    )
  )
}


#' B-Spline Basis Expansion for RHC Covariates (opt-in)
#'
#' Applies a cubic B-spline basis (with \code{knots_per_var} internal knots)
#' to each continuous column of \code{X} and leaves binary / one-hot columns
#' unchanged. Provided as a reference \code{phi} for
#' \code{\link{build_rhc_data_split}} when the user wishes to test RoCE
#' in the high-dimensional regime; \emph{not} applied by default.
#'
#' @param X Numeric matrix of raw covariates.
#' @param knots_per_var Integer. Internal knot count for each continuous
#'   column (default 5).
#' @param degree Integer. Spline degree (default 3, cubic).
#' @return A numeric matrix whose columns are the spline basis of each
#'   continuous column concatenated with the original binary / one-hot
#'   columns.
#' @export
phi_rhc_bspline <- function(X, knots_per_var = 5L, degree = 3L) {
  if (!is.matrix(X)) stop("X must be a numeric matrix.")
  requireNamespace("splines", quietly = TRUE)
  knots_per_var <- as.integer(knots_per_var)
  degree        <- as.integer(degree)

  # A column is treated as continuous if it is not binary-valued.
  col_is_binary <- apply(X, 2L, function(col) {
    vals <- unique(col[!is.na(col)])
    length(vals) <= 2L && all(vals %in% c(0, 1))
  })

  expanded <- vector("list", ncol(X))
  for (j in seq_len(ncol(X))) {
    col <- X[, j]
    if (col_is_binary[j]) {
      expanded[[j]] <- matrix(col, ncol = 1L, dimnames = list(NULL, colnames(X)[j]))
    } else {
      probs <- seq(0, 1, length.out = knots_per_var + 2L)[-c(1L, knots_per_var + 2L)]
      knots <- stats::quantile(col, probs = probs, na.rm = TRUE, names = FALSE)
      bs    <- splines::bs(col, knots = knots, degree = degree, intercept = FALSE)
      colnames(bs) <- paste0(colnames(X)[j], "_bs", seq_len(ncol(bs)))
      expanded[[j]] <- bs
    }
  }
  do.call(cbind, expanded)
}
