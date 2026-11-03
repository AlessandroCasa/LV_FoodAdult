# Scenario definitions ----------------------------------------------------
# Choose scenario_id in run_simulations.R; this file only defines the cases.

scenarios <- data.frame(
  ID = 1:19,
  name = c("S1_n100_p50", "S1_n100_p100", "S1_n100_p150",
           "S1_n250_p50", "S1_n250_p100", "S1_n250_p150",
           "S1_n500_p50", "S1_n500_p100", "S1_n500_p150",
           "S2_hat_weak", "S2_hat_strong", "S2_spline_weak", "S2_spline_strong",
           "S3_weak_low", "S3_weak_medium", "S3_strong_low", "S3_strong_medium",
           "AR1", "Student_t5"),
  n = c(rep(c(100, 250, 500), each = 3), rep(250, 10)),
  p = c(rep(c(50, 100, 150), 3), rep(100, 10)),
  shape = c(rep("hat", 11), "spline", "spline", rep("hat", 6)),
  strength = c(rep("standard", 9), "weak", "strong", "weak", "strong",
               "weak", "weak", "strong", "strong", "standard", "standard"),
  adulteration = c(rep("mixed", 13), "low", "medium", "low", "medium", "mixed", "mixed"),
  error = c(rep("normal", 18), "student"),
  correlation = c(rep(0, 17), 0.9, 0),
  stringsAsFactors = FALSE
)

define_scenario <- function(scenario_id) {
  stopifnot(length(scenario_id) == 1, scenario_id %in% scenarios$ID)
  setting <- scenarios[scenarios$ID == scenario_id, ]
  n <- setting$n
  p <- setting$p

  # Signal ---------------------------------------------------------------
  if (setting$shape == "hat") {
    w <- seq(from = -5, to = 5, length.out = p)
    deltaTrue <- (1 - w^2) * exp(-w^2 / 2)
    deltaTrue <- ifelse(deltaTrue > 0.6, 0.6, deltaTrue)
    deltaTrue <- ifelse(deltaTrue < -0.4, -0.4, deltaTrue)
    deltaTrue <- ifelse(deltaTrue > -0.1 & deltaTrue < 0, 0, deltaTrue)
    multiplier <- 4
    if (setting$strength == "weak") multiplier <- 2
    if (setting$strength == "strong") multiplier <- 6
    deltaTrue <- multiplier * deltaTrue
    mu_p <- rep(0, p)
    spline_details <- NULL
  } else {
    # Fixed curves across Monte Carlo replications; same coefficients in
    # the weak and strong cases, with a narrower/wider affected region.
    set.seed(20241100)
    s <- seq(0, 1, length.out = p)
    knots <- seq(0, 1, length.out = 17)[2:16]
    basis <- splines::bs(s, knots = knots, degree = 3, intercept = FALSE,
                        Boundary.knots = c(0, 1))
    beta_coef <- runif(ncol(basis), min = -2, max = 2)
    increments <- runif(ncol(basis), min = 2, max = 4)
    changed <- if (setting$strength == "weak") 9 else 8:10
    shift_coef <- rep(0, ncol(basis))
    shift_coef[changed] <- increments[changed]
    mu_p <- as.vector(basis %*% beta_coef)
    # This is the difference between the adulterant and pure profiles.
    deltaTrue <- as.vector(basis %*% shift_coef)
    spline_details <- list(knots = knots, beta_coef = beta_coef,
                           shift_coef = shift_coef, changed = changed)
  }

  # Adulteration levels and noise ----------------------------------------
  g_levels <- c(0.1, 0.2, 0.3)
  if (setting$adulteration == "low") g_levels <- 0.05
  if (setting$adulteration == "medium") g_levels <- 0.2

  sigma <- 0.1
  SigmaTrue <- sigma^2 * setting$correlation^abs(outer(1:p, 1:p, "-"))
  if (setting$correlation == 0) {
    OmegaTrue <- solve(SigmaTrue)
  } else {
    # Analytic AR(1) inverse: exact zeros and both off-diagonal bands.
    rho <- setting$correlation
    OmegaTrue <- matrix(0, nrow = p, ncol = p)
    diag(OmegaTrue) <- c(1, rep(1 + rho^2, p - 2), 1) /
      (sigma^2 * (1 - rho^2))
    j <- 1:(p - 1)
    OmegaTrue[cbind(j, j + 1)] <- -rho / (sigma^2 * (1 - rho^2))
    OmegaTrue[cbind(j + 1, j)] <- -rho / (sigma^2 * (1 - rho^2))
  }

  list(setting = setting, n = n, p = p, mu_p = mu_p, deltaTrue = deltaTrue,
       g_levels = g_levels, sigma = sigma, SigmaTrue = SigmaTrue,
       OmegaTrue = OmegaTrue, spline_details = spline_details)
}
