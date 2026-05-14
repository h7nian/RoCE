# helper-load.R - Auto-loaded by testthat before any test files
#
# Loads the FACEC code so all exported functions, constants, and compiled C++
# code are available in the test environment. Internal (non-exported) objects
# are attached for convenience. Shared synthetic-data factories used across
# multiple test files are defined here as well, so they live on the same
# search path as the FACEC internals they call.
#
# By default, prefer the *local* source tree via devtools::load_all() when the
# repo is present, even if an installed FACEC package exists. This avoids the
# common pitfall of running tests against an outdated installed version.
# To force testing the installed package, set FACEC_TEST_INSTALLED=1.

force_installed <- Sys.getenv("FACEC_TEST_INSTALLED", "0") %in% c("1", "TRUE", "true", "True")

is_facec_repo <- function(path) {
  desc <- file.path(path, "DESCRIPTION")
  if (!file.exists(desc)) return(FALSE)
  hdr <- tryCatch(readLines(desc, n = 20L, warn = FALSE), error = function(e) character(0))
  any(grepl("^Package:\\s*FACEC\\s*$", hdr))
}

repo_root <- if (is_facec_repo(getwd())) {
  normalizePath(getwd(), mustWork = FALSE)
} else {
  normalizePath(file.path(dirname(getwd()), ".."), mustWork = FALSE)
}

if (!force_installed && is_facec_repo(repo_root) && requireNamespace("devtools", quietly = TRUE)) {
  # IMPORTANT: helpers=FALSE prevents devtools from sourcing testthat helper
  # files (including this one), which would otherwise recurse indefinitely.
  devtools::load_all(repo_root, helpers = FALSE)
} else if (requireNamespace("FACEC", quietly = TRUE)) {
  library(FACEC)
} else if (requireNamespace("devtools", quietly = TRUE) && is_facec_repo(repo_root)) {
  devtools::load_all(repo_root, helpers = FALSE)
} else {
  stop("Tests require either the FACEC source tree + devtools, or an installed FACEC package.")
}

# Attach non-exported package internals in an isolated search-path environment.
# This keeps test convenience (unqualified access to internal helpers/C++ bindings)
# while avoiding .GlobalEnv pollution and mask/conflict noise.
pkg_env <- asNamespace("FACEC")
exported <- getNamespaceExports("FACEC")
all_public_names <- ls(pkg_env, all.names = FALSE)
internal_names <- setdiff(all_public_names, exported)

if ("FACEC_test_internals" %in% search()) {
  detach("FACEC_test_internals", character.only = TRUE)
}

test_internal_env <- new.env(parent = emptyenv())
for (nm in internal_names) {
  obj <- tryCatch(get(nm, envir = pkg_env), error = function(e) NULL)
  if (!is.null(obj)) {
    assign(nm, obj, envir = test_internal_env)
  }
}

# ============================================================================
# Shared synthetic-data factories (used across multiple test files)
# ----------------------------------------------------------------------------
# These live in helper-load.R rather than a separate helper file because
# testthat's helper-file discovery order and target environment vary by
# version; defining them in the same file that drives the package load and
# the internals attach guarantees they end up on the same search path as
# FACEC internals like MIN_TREATED_FOR_MODEL.
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
    cache$data_split <- get("make_small_data_split", envir = test_internal_env)(
      n_per_site = 40L, K = 1L, p = 4L, seed = 2026L, n_folds = 3L,
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

attach(test_internal_env, name = "FACEC_test_internals", warn.conflicts = FALSE)
