LS_partial_est <- function(X, g_init, g_known = NULL, tol = 10^(-3), 
                           mean_pure = NULL) 
{
  # load packages and functions needed 
  require(rio)
  SS_f <- function(x, h, mu, delta) {
    out <- sum((x- mu - h*delta)^2)
    return(out)
  }
  
  # semi-supervision 
  if (!is.null(g_known)) {
    sub_known <- which(!is.na(g_known))
    g_init[sub_known] <- g_known[sub_known]
  } else {
    sub_known <- NULL
  }
  
  # initialize quantities
  N <- nrow(X)
  P <- ncol(X)
  G <- matrix(g_init, N, P, byrow = F)
  # here trying to find clever initial value for mu considering potential (semi)supervision
  if (!is.null(mean_pure)) {
    mu_init <- mean_pure 
  } else {
    if (sum(g_known == 0, na.rm = T) != 0) {
      if (length(which(g_known == 0)) > 1) {
      mu_init <- apply(X[which(g_known == 0), ], 2, mean)
      } else {
        mu_init <- apply(X[which(g_known < quantile(g_known, na.rm = T, probs = 0.2)), ], 2, mean)
      }
    } else {
      mu_init <- apply(X, 2, mean)
    }
  }
  MU_p <- matrix(mu_init, nrow = N, ncol = P, byrow = T)
  
  criterion <- TRUE
  SSold <- Inf 
  tol <- tol
  
  
  while (criterion) {
    delta <- apply(G * (X - MU_p), 2, sum)/apply(G^2, 2, sum)
    
    DELTA <- matrix(delta, nrow = N, ncol = P, byrow = T)
    if (is.null(mean_pure)) {
      mu_p <- apply(X - G * DELTA, 2, mean)
    } else {
      mu_p <- mean_pure
    }
    MU_p <- matrix(mu_p, nrow = N, ncol = P, byrow = T)
    
    g <- g_init 
    ghat <- rep(NA, N)
    for (n in 1:N) {
      if (!is.null(g_known) & n %in% sub_known) {
        ghat[n] <- g[n]
      } else {
        fit <- nlminb(g[n], SS_f, mu = mu_p, delta = delta,
                      x = X[n,], lower = 0, upper = 1)
        ghat[n] <- fit$par 
      }
    }
    g_init <- ghat
    G <- matrix(g_init, N, P, byrow = F)
    
    SSnew <- sum((X - MU_p - G*DELTA)^2)
    criterion <- ((SSold - SSnew) > tol)
    SSold <- SSnew
  }
  
  out <- list(pureMean = mu_p, meanShift = delta, membership = ghat)
  return(out)
  
}