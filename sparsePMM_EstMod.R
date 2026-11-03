# Revised PMM: sum Gaussian log-likelihood minus uniform L1 penalties.
# Only this ordered sparse-fused-lasso D is supported, so its df are explicit.
sparse_fused_D <- function(p) rbind(diff(diag(p)), diag(p))
soft_threshold <- function(x, a) sign(x) * pmax(abs(x) - a, 0)
solve_spd <- function(A, b) {
  R <- chol(A)
  backsolve(R, forwardsolve(t(R), b))
}

gaussian_loglik <- function(X, g, delta, Omega) {
  residual <- X - outer(g, delta)
  n <- nrow(X); p <- ncol(X)
  n / 2 * (2 * sum(log(diag(chol(Omega)))) - p * log(2 * pi)) -
    sum((residual %*% Omega) * residual) / 2
}

penalized_objective <- function(X, g, delta, Omega, D, known,
                                lambda_delta, lambda_g, lambda_omega) {
  gaussian_loglik(X, g, delta, Omega) - lambda_delta * sum(abs(D %*% delta)) -
    lambda_g * sum(g[is.na(known)]) - lambda_omega * sum(abs(Omega))
}

membership_update <- function(X, delta, Omega, lambda_g, known, upper_g) {
  denom <- drop(crossprod(delta, Omega %*% delta))
  g <- if (denom <= .Machine$double.eps) rep(0, nrow(X)) else
    pmin(upper_g, pmax(0, (drop(X %*% Omega %*% delta) - lambda_g) / denom))
  g[!is.na(known)] <- known[!is.na(known)]
  g
}

initialize_pmm <- function(X, known, upper_g = 0.5) {
  labelled <- which(!is.na(known))
  if (!any(known[labelled] > 0)) stop("A positive known concentration is required.")
  g <- numeric(nrow(X)); g[labelled] <- known[labelled]
  delta <- drop(crossprod(g, X)) / sum(g^2)
  old <- Inf
  for (it in seq_len(300L)) {
    g <- membership_update(X, delta, diag(ncol(X)), 0, known, upper_g)
    delta <- drop(crossprod(g, X)) / sum(g^2)
    err <- sum((X - outer(g, delta))^2)
    if (is.finite(old) && abs(old - err) <= 1e-8 * max(1, old)) break
    old <- err
  }
  variance <- colMeans((X - outer(g, delta))^2)
  list(g = g, delta = delta, Omega = diag(1 / pmax(variance, 1e-6)))
}

# Exact zero/equality constraints selected by the ADMM auxiliary variable.
fused_structure <- function(z, p, tol = 1e-8) {
  fused <- abs(z[seq_len(p - 1L)]) <= tol
  group <- cumsum(c(TRUE, !fused))
  zero <- abs(z[(p):(2L * p - 1L)]) <= tol
  active <- vapply(split(zero, group), function(x) !any(x), logical(1))
  keep <- which(active)
  basis <- matrix(0, p, length(keep))
  if (length(keep)) for (j in seq_along(keep)) basis[group == keep[j], j] <- 1
  list(basis = basis, df = length(keep))
}

delta_update <- function(X, g, Omega, D, lambda, initial, control) {
  H <- sum(g^2) * Omega
  b <- drop(Omega %*% crossprod(X, g))
  if (lambda == 0) {
    d <- drop(solve_spd(H, b))
    return(list(delta = d, z = drop(D %*% d), df = length(d),
                converged = TRUE, iterations = 0L, polished = TRUE))
  }
  p <- ncol(X); m <- nrow(D)
  # rho is a numerical ADMM parameter, unrelated to glasso's rho_omega.
  rho <- max(1, median(diag(H)) / 4)
  delta <- initial; z <- drop(D %*% delta); u <- numeric(m)
  DtD <- crossprod(D)
  R <- chol(H + rho * DtD)
  converged <- FALSE
  for (it in seq_len(control$admm_maxit)) {
    rhs <- b + rho * drop(crossprod(D, z - u))
    delta <- drop(backsolve(R, forwardsolve(t(R), rhs)))
    Dz <- drop(D %*% delta); previous <- z
    z <- soft_threshold(Dz + u, lambda / rho)
    u <- u + Dz - z
    primal <- sqrt(sum((Dz - z)^2))
    dual <- rho * sqrt(sum(crossprod(D, z - previous)^2))
    ep <- sqrt(m) * control$admm_abs + control$admm_rel *
      max(sqrt(sum(Dz^2)), sqrt(sum(z^2)))
    ed <- sqrt(p) * control$admm_abs + control$admm_rel *
      sqrt(sum((rho * crossprod(D, u))^2))
    if (primal <= ep && dual <= ed) { converged <- TRUE; break }
    if (it %% 25L == 0L) {
      newrho <- if (primal > 10 * dual) rho * 2 else
        if (dual > 10 * primal) rho / 2 else rho
      if (newrho != rho) {
        u <- u * rho / newrho; rho <- newrho
        R <- chol(H + rho * DtD)
      }
    }
  }
  structure <- fused_structure(z, p, control$zero_tol)
  # Polish on the selected zero/fusion constraints, keeping the L1 penalty.
  # This is not an unpenalized post-selection refit.
  B <- structure$basis
  polished <- if (!ncol(B)) numeric(p) else drop(B %*% solve_spd(
    crossprod(B, H %*% B), crossprod(B, b - lambda * drop(crossprod(D, sign(z))))))
  q <- function(d) drop(crossprod(d, H %*% d) / 2 - crossprod(b, d)) +
    lambda * sum(abs(D %*% d))
  dp <- drop(D %*% polished)
  nz <- abs(z) > control$zero_tol
  # A tiny nonzero ADMM entry may become exactly zero when its fused group
  # is anchored at zero. Reaching that boundary is valid; crossing it is not.
  sign_ok <- all(sign(dp[nz]) == sign(z[nz]) | abs(dp[nz]) <= control$zero_tol)
  accepted <- converged && sign_ok &&
    q(polished) <= q(delta) + 1e-6 * max(1, abs(q(delta))) &&
    max(abs(polished - delta)) <= 1e-4 * max(1, max(abs(delta)))
  if (accepted) delta <- polished
  reason <- if (!converged) "admm_maxit" else if (!sign_ok) "polish_sign" else
    if (!accepted) "polish_tolerance" else "converged"
  list(delta = delta, z = drop(D %*% delta), df = structure$df,
       converged = converged && accepted, iterations = it, polished = accepted,
       reason = reason)
}

precision_update <- function(residual, rho_omega, control) {
  # Mean residual is already specified as zero by the model. Do NOT use cov():
  # cov() subtracts another mean and divides by n-1.
  S <- crossprod(residual) / nrow(residual)
  fit <- glasso::glasso(s = S, rho = rho_omega, penalize.diagonal = TRUE,
                        thr = control$glasso_thr, maxit = control$glasso_maxit)
  Omega <- (fit$wi + t(fit$wi)) / 2
  chol(Omega)  # Fail explicitly if not positive definite.
  list(Omega = Omega, converged = fit$errflag == 0 &&
         fit$niter < control$glasso_maxit, iterations = fit$niter)
}

parameter_count <- function(g, Omega, known, df_delta, zero_tol = 1e-8) {
  c(delta = df_delta,
    omega = ncol(Omega) + sum(abs(Omega[upper.tri(Omega)]) > zero_tol),
    g = sum(g[is.na(known)] > zero_tol))
}

sparsePMM <- function(X, g_known, lambda_delta, rho_omega, lambda_g,
                      control, upper_g = 0.5, initial = NULL) {
  lambda_delta <- unname(lambda_delta)
  rho_omega <- unname(rho_omega)
  lambda_g <- unname(lambda_g)
  stopifnot(is.matrix(X), all(is.finite(X)), ncol(X) >= 2,
            length(g_known) == nrow(X), lambda_delta >= 0,
            rho_omega > 0, lambda_g >= 0,
            all(g_known[!is.na(g_known)] >= 0),
            all(g_known[!is.na(g_known)] <= upper_g))
  if (is.null(initial)) initial <- initialize_pmm(X, g_known, upper_g)
  g <- initial$g; delta <- initial$delta; Omega <- initial$Omega
  D <- sparse_fused_D(ncol(X)); lambda_omega <- nrow(X) * rho_omega / 2
  objective <- function(g, d, O) penalized_objective(
    X, g, d, O, D, g_known, lambda_delta, lambda_g, lambda_omega)
  old <- objective(g, delta, Omega); trace <- old
  converged <- FALSE; reason <- "outer_maxit"
  admm_iterations <- 0L; max_drop <- 0
  for (it in seq_len(control$maxit)) {
    previous <- list(g = g, delta = delta, Omega = Omega)
    upd <- delta_update(X, g, Omega, D, lambda_delta, delta, control)
    admm_iterations <- admm_iterations + upd$iterations
    if (!upd$converged) { reason <- upd$reason; break }
    delta <- upd$delta
    g <- membership_update(X, delta, Omega, lambda_g, g_known, upper_g)
    oup <- precision_update(X - outer(g, delta), rho_omega, control)
    if (!oup$converged) { reason <- "glasso_not_converged"; break }
    Omega <- oup$Omega
    value <- objective(g, delta, Omega)
    max_drop <- max(max_drop, old - value)
    if (!is.finite(value) || value < old - control$objective_tol * max(1, abs(old))) {
      reason <- "objective_decreased"; break
    }
    trace <- c(trace, value)
    change <- max(max(abs(g - previous$g)),
                  max(abs(delta - previous$delta)) / max(1, max(abs(previous$delta))),
                  norm(Omega - previous$Omega, "F") / max(1, norm(previous$Omega, "F")))
    if (abs(value - old) / max(1, abs(old)) <= control$outer_tol &&
        change <= control$parameter_tol) {
      converged <- TRUE; reason <- "converged"; break
    }
    old <- value
  }
  df_delta <- if (lambda_delta == 0) ncol(X) else
    fused_structure(drop(D %*% delta), ncol(X), control$zero_tol)$df
  counts <- parameter_count(g, Omega, g_known, df_delta, control$zero_tol)
  ll <- gaussian_loglik(X, g, delta, Omega)
  list(membership = g, meanShift = delta, precision = Omega,
       gamma = drop(D %*% delta), df = counts, nParam = sum(counts),
       loglik_unpenalized = ll, objective = objective(g, delta, Omega),
       BIC = 2 * ll - sum(counts) * log(nrow(X)),
       hyperparameters = c(lambda_delta = lambda_delta, rho_omega = rho_omega,
                           lambda_g = lambda_g, lambda_omega = lambda_omega),
       converged = converged, reason = reason, nIter = it,
       admm_iterations = admm_iterations, objective_trace = trace,
       max_objective_drop = max_drop)
}
