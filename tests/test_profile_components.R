if (!requireNamespace("RTMB", quietly = TRUE)) {
  stop("The RTMB package is required for this test.", call. = FALSE)
}

source(file.path("R", "profile_components.R"))

toy_nll <- function(par) {
  part_a <- 0.5 * (par$theta - 2)^2
  part_b <- 0.5 * (par$beta + par$theta - 1)^2
  RTMB::REPORT(part_a)
  RTMB::REPORT(part_b)
  part_a + part_b
}

toy_components <- function(report) {
  c(A = report$part_a, B = report$part_b)
}

obj <- RTMB::MakeADFun(
  toy_nll,
  parameters = list(theta = 0, beta = 0)
)
fit <- nlminb(obj$par, obj$fn, obj$gr)
stopifnot(fit$convergence == 0L)

grid <- c(1.5, 2, 2.5)
slice <- profile_components_slice(
  obj, "theta", grid,
  component_fun = toy_components,
  par = fit$par
)
stopifnot(nrow(slice) == length(grid))
stopifnot(all(abs(slice$objective - rowSums(slice[c("A", "B", "Other")])) < 1e-10))
stopifnot(all(abs(slice$Other) < 1e-10))

profile <- profile_components_reopt(
  obj, "theta", grid,
  component_fun = toy_components,
  par = fit$par,
  trace = FALSE
)
stopifnot(all(profile$convergence == 0L))
stopifnot(all(profile$max_gradient < 1e-6))
stopifnot(all(abs(profile$objective - 0.5 * (grid - 2)^2) < 1e-8))
stopifnot(all(abs(profile$objective - rowSums(profile[c("A", "B", "Other")])) < 1e-10))

duplicate_par <- c(x = 1, x = 2, y = 3)
stopifnot(profile_resolve_parameter(duplicate_par, "x[2]")$index == 2L)
ambiguous <- try(profile_resolve_parameter(duplicate_par, "x"), silent = TRUE)
stopifnot(inherits(ambiguous, "try-error"))

long <- profile_components_long(profile)
stopifnot(all(c("component", "nll", "delta") %in% names(long)))
stopifnot(nrow(long) == nrow(profile) * 4L)
stopifnot("Objective" %in% long$component)

if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("ggthemes", quietly = TRUE)) {
  figure <- plot_profile_components(profile)
  stopifnot(inherits(figure, "ggplot"))
  stopifnot(inherits(figure$theme, "theme"))
}

cat("Likelihood-profile helper tests passed.\n")
