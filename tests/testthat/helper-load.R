# helper-load.R - Auto-loaded by testthat before any test files
#
# Loads the RoCE code so all exported functions, constants, and compiled C++
# code are available in the test environment. Internal (non-exported) objects
# are attached for convenience. Shared synthetic-data factories used across
# multiple test files are defined here as well, so they live on the same
# search path as the RoCE internals they call.
#
# By default, prefer the *local* source tree via devtools::load_all() when the
# repo is present, even if an installed RoCE package exists. This avoids the
# common pitfall of running tests against an outdated installed version.
# To force testing the installed package, set ROCE_TEST_INSTALLED=1.

force_installed <- Sys.getenv("ROCE_TEST_INSTALLED", "0") %in% c("1", "TRUE", "true", "True")

is_roce_repo <- function(path) {
  desc <- file.path(path, "DESCRIPTION")
  if (!file.exists(desc)) return(FALSE)
  hdr <- readLines(desc, n = 20L)
  any(grepl("^Package:\\s*RoCE\\s*$", hdr))
}

repo_root <- if (is_roce_repo(getwd())) {
  normalizePath(getwd(), mustWork = FALSE)
} else {
  normalizePath(file.path(dirname(getwd()), ".."), mustWork = FALSE)
}

if (force_installed) {
  if (!requireNamespace("RoCE", quietly = TRUE)) {
    stop("ROCE_TEST_INSTALLED=1 requires an installed RoCE package; source compilation is disabled.")
  }
  library(RoCE)
} else if (is_roce_repo(repo_root) && requireNamespace("devtools", quietly = TRUE)) {
  # IMPORTANT: helpers=FALSE prevents devtools from sourcing testthat helper
  # files (including this one), which would otherwise recurse indefinitely.
  devtools::load_all(repo_root, helpers = FALSE)
} else if (requireNamespace("RoCE", quietly = TRUE)) {
  library(RoCE)
} else if (requireNamespace("devtools", quietly = TRUE) && is_roce_repo(repo_root)) {
  devtools::load_all(repo_root, helpers = FALSE)
} else {
  stop("Tests require either the RoCE source tree + devtools, or an installed RoCE package.")
}

# Attach non-exported package internals in an isolated search-path environment.
# This keeps test convenience (unqualified access to internal helpers/C++ bindings)
# while avoiding .GlobalEnv pollution and mask/conflict noise.
pkg_env <- asNamespace("RoCE")
exported <- getNamespaceExports("RoCE")
all_namespace_names <- ls(pkg_env, all.names = TRUE)
reserved_namespace_names <- c(
  ".__NAMESPACE__.", ".__S3MethodsTable__.", ".packageName"
)
internal_names <- setdiff(
  all_namespace_names,
  c(exported, reserved_namespace_names)
)

if ("ROCE_test_internals" %in% search()) {
  detach("ROCE_test_internals", character.only = TRUE)
}

test_internal_env <- new.env(parent = emptyenv())
for (nm in internal_names) {
  assign(nm, get(nm, envir = pkg_env), envir = test_internal_env)
}

# ============================================================================
# Shared synthetic-data factories (used across multiple test files)
# ----------------------------------------------------------------------------
# These live in helper-load.R rather than a separate helper file because
# testthat's helper-file discovery order and target environment vary by
# version; defining them in the same file that drives the package load and
# the internals attach guarantees they end up on the same search path as
# RoCE internals like MIN_TREATED_FOR_MODEL.
#
# IMPORTANT: do not call generate_simulation_data with p < 4. The function
# transform_covariates() in R/data_generation.R indexes X[, 4] unconditionally
# for transform_type "mild" and "strong". The factory below enforces p >= 4.
# ============================================================================

assign("make_small_data_split", function(n_per_site = 40L, K = 1L, p = 4L,
                                         seed = 42L, n_folds = 3L,
                                         A_val = 1L, max_tries = 40L,
                                         outcome_type = "continuous") {
  if (p < 4L) {
    stop("make_small_data_split: p must be >= 4 (transform_covariates indexes X[, 4]).",
         call. = FALSE)
  }
  min_treated <- max(n_folds * 3L, MIN_TREATED_FOR_MODEL)

  for (attempt in seq_len(max_tries)) {
    set.seed(seed + attempt - 1L)
    n_total <- n_per_site * (K + 1L)
    data <- generate_simulation_data(
      n_total, K, p,
      config       = "C1",
      outcome_type = outcome_type,
      dgp_type     = "roce",   # estimator smoke factory: pin to the RoCE DGP
      warn_ignored = FALSE
    )
    data_split <- split_data_by_site(data)

    treated_ok <- vapply(names(data_split), function(site) {
      sum(data_split[[site]]$A == A_val) >= min_treated
    }, logical(1))

    if (all(treated_ok)) {
      attr(data_split, "raw_data") <- data
      return(data_split)
    }
  }

  stop("make_small_data_split: failed to find a balanced data split after ",
       max_tries, " attempts.", call. = FALSE)
}, envir = test_internal_env)

# Shared smoke cache for run_crossfit smoke results. Lives inside
# test_internal_env so that subsequent helper-side functions (get_smoke_*)
# can find it via the attached search-path entry.
assign(".shared_smoke_cache", new.env(parent = emptyenv()),
       envir = test_internal_env)

assign("get_smoke_data_split", function() {
  cache <- get(".shared_smoke_cache", envir = test_internal_env)
  if (!exists("data_split", envir = cache, inherits = FALSE)) {
    # Two-level cross-fitting excludes two outer folds when constructing an
    # initial nuisance fit. With only 40 observations per site, the remaining
    # treatment-arm inner-CV training samples can have fewer rows than the
    # intercept-augmented five-parameter model. Use a still-small but properly
    # identified integration fixture; production simulations use 1000/site.
    cache$data_split <- get("make_small_data_split", envir = test_internal_env)(
      n_per_site = 120L, K = 1L, p = 4L, seed = 2026L, n_folds = 3L,
      outcome_type = "continuous"
    )
  }
  cache$data_split
}, envir = test_internal_env)

assign("get_smoke_result", function(mode = c("two_round", "one_round")) {
  mode <- match.arg(mode)
  cache <- get(".shared_smoke_cache", envir = test_internal_env)
  key <- paste0("res_", mode)
  if (!exists(key, envir = cache, inherits = FALSE)) {
    data_split <- get("get_smoke_data_split", envir = test_internal_env)()
    cache[[key]] <- run_crossfit(
      data_split,
      n_folds            = 3L,
      communication_mode = mode,
      lambda_selection   = 0.1,
      verbose            = FALSE,
      nlambda_init       = 3L,
      use_lambda_cache   = TRUE,
      family             = "gaussian"
    )
  }
  cache[[key]]
}, envir = test_internal_env)

# Repository root for tests that source repository-only scripts (scripts/slurm).
# It lives on the attached search-path environment because the helper's own
# environment is not visible to the test files after load_all().
assign("repo_root", repo_root, envir = test_internal_env)

attach(test_internal_env, name = "ROCE_test_internals", warn.conflicts = FALSE)
