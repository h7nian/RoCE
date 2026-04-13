# helper-load.R - Auto-loaded by testthat before any test files
#
# Loads the FACEC code so all exported functions, constants, and compiled C++
# code are available in the test environment. Internal (non-exported) objects
# are attached for convenience.
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

attach(test_internal_env, name = "FACEC_test_internals", warn.conflicts = FALSE)
