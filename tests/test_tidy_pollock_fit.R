source(file.path("R", "tidy-pollock-fit.R"))

fixture <- list(
  label = "fixture",
  objective = 10,
  evaluated_total = 10,
  convergence = 0L,
  max_gradient = 1e-6,
  seconds = 1,
  fixed_parameters = list(log_q = 0),
  parameter_table = data.frame(
    term = "log_q", estimate = 0, std.error = 0.1, gradient = 1e-6,
    estimation_type = "fixed_effect", module = NA_character_, year = NA_integer_,
    age = NA_integer_, engine = "custom RTMB", model = "fixture"
  ),
  n_parameters = 1L,
  report = list(
    N = matrix(c(10, 5, 12, 6), 2, 2, dimnames = list(2000:2001, 1:2)),
    SSB = c(8, 9), obs_catch = c(2, 3), pred_catch = c(2.1, 2.9),
    oac_fsh = matrix(c(.4, .6), 1), phat_fsh = matrix(c(.5, .5), 1),
    oac_bts = matrix(numeric(), 0, 2), phat_bts = matrix(numeric(), 0, 2),
    oac_ats = matrix(numeric(), 0, 2), phat_ats = matrix(numeric(), 0, 2),
    yrs_fsh_data = 2001L, yrs_bts_data = integer(), yrs_ats_data = integer(),
    sam_fsh = 100, sam_bts = numeric(), sam_ats = numeric(), tot_like = 10
  )
)

fit <- as_pollock_fit_rtmb(fixture)
stopifnot(inherits(fit, "pollock_fit"))
stopifnot(nrow(generics::tidy(fit)) == 1L)
stopifnot(nrow(generics::tidy(fit, parameters = "derived_quantity")) == 4L)
stopifnot(nrow(generics::glance(fit)) == 1L)
stopifnot(all(c(".truth", ".pred", "data_type", "fleet") %in% names(generics::augment(fit))))
stopifnot(sum(generics::augment(fit)$data_type == "catch") == 2L)
stopifnot(sum(generics::augment(fit)$data_type == "age_composition") == 2L)

cat("tidy pollock interface tests passed\n")
