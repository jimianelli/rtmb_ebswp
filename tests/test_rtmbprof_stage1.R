if (!requireNamespace("RTMBprof", quietly = TRUE)) {
  stop("The RTMBprof package is required for this test.", call. = FALSE)
}

repo_root <- normalizePath(".", mustWork = TRUE)
rtmb_env <- new.env(parent = globalenv())
rtmb_env$rm <- function(...) invisible(NULL)
rtmb_env$source <- function(file, ...) {
  base::source(file, local = parent.frame(), ...)
}
base::source(file.path(repo_root, "R", "config.R"), local = rtmb_env)
base::source(file.path(repo_root, "R", "profile_components.R"))

obj <- rtmb_env$obj
objective <- obj$fn(obj$par)
report <- obj$report(obj$par)
components <- pollock_profile_components(report)

stopifnot(
  "avgsel_like" %in% names(report),
  all(c("Bzero", "q_bts", "q_ats") %in% names(report)),
  all(is.finite(components)),
  abs(objective - sum(components)) <=
    1e-6 + 1e-8 * max(1, abs(objective), abs(sum(components)))
)

cat("RTMBprof Stage 1 component closure test passed.\n")
