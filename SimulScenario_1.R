require(LaplacesDemon)
require(mvtnorm)
require(glasso)
require(patchwork)
source("LS_partial_est.R")
source("penPMM.R")
plot_matrix <- function(matrix, highlight_sparsity = FALSE) {
  require(ggplot2) 
  require(patchwork)
  
  matrix_long <- reshape2::melt(matrix[nrow(matrix):1, , drop = FALSE])
  
  plot_mat <- ggplot(matrix_long, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill = value)) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0
    ) +
    labs(x = "", y = "", fill = "") +
    theme_void()
  if (highlight_sparsity == FALSE) {
    return(plot_mat)
  }
  
  plot_mat +
    geom_point(
      data = dplyr::filter(matrix_long, value == 0),
      pch = "x",
      stroke = 4
    )
  
}
createD <- function(k, p) {
  require(genlasso)
  
  ordvec <- seq(from = 0, to = k, by = 1)
  mat <- matrix(NA, nrow = k+1, ncol = p)
  for (j in ordvec) {
    if (j == 0) {
      mat[j+1,] <- c(1, rep(0, p-1))
    } else {
      mat[j+1,] <- getDtf(n = p, ord = j-1)[1,]
    }
  }
  out <- rbind(mat, getDtf(n = p, ord = k))
  return(out)
}
gen_sparse_posdef_mat_erdos_renyi <- function(prob_of_connections, V) {
  sample_graph <-
    igraph::as_adjacency_matrix(igraph::sample_gnp(n = V, p = prob_of_connections), sparse = FALSE)
  
  starting_matrix <- matrix(.9, nrow = V, ncol = V)
  diag(starting_matrix) <- 1
  
  mixggm:::icf(
    Sigma = diag(diag(starting_matrix), nrow = V, ncol = V),
    starting_matrix,
    sample_graph,
    10,
    mixggm:::controlICF()$tol,
    mixggm:::controlICF()$itMax,
    FALSE,
    FALSE,
    0
  )$sigma
  
}



#  Setting parameters --------------------------------------------------------

# sample size
n <- 250
# number of variables
p <- 100

# modified mexican hat definition
w <- seq(from = -5, to = 5, length.out = p)
deltaTrue <- (1 - w^2)*exp(-(w^2)/2)
deltaTrue <- ifelse(deltaTrue > 0.6, 0.6, deltaTrue)
deltaTrue <- ifelse(deltaTrue < -0.4, -0.4, deltaTrue)
deltaTrue <- ifelse(deltaTrue > -0.1 & deltaTrue < 0, 0, deltaTrue)
deltaTrue <- deltaTrue*4
mu_p <- rep(0, p)

# range of the membership
range <- c(0.1,0.2,0.3)
# probability of pure food
prob_noAd <- 0.7 

# percentage of supervision
percSup <- 0.2

# generate fixed precision
eps <- 100
OmegaTrue <- diag(eps, nrow = p)

# choose structure for D (generalized lasso)
D <- createD(k = 0, p = p)
D <- D[-1, ]
D <- rbind(D, diag(1, nrow = p))

# number of MonteCarlo samples
B <- 100 

# object to be filled
mae <- acc <- sens <- spec <- mse <- vector(length = B)
accDelta <- sensDelta <- specDelta <- vector(length = B)

for (ind in 1:B) {
  set.seed(ind)
  
  # generate true membership labels 
  g_true <- sample(range, size = n, replace = T)
  g_true <- g_true*rbinom(n = n, size = 1, prob = 1 - prob_noAd)
  
  # generate the data
  meanMat <- outer(g_true, deltaTrue, FUN = "*")
  dat <- rmvnp(n = n, mu = meanMat, Omega = OmegaTrue)
  # center the data wrt pure mean 
  X <- t(apply(dat, 1, function(x) x - mu_p))
  
  # set degree of supervision
  g_known <- rep(NA, length = length(g_true))
  sub1 <- sample(which(g_true == 0), 
                 size = round(percSup*sum(g_true == 0)), replace = F)
  sub2 <- sample(which(g_true != 0), 
                 size = round(percSup*sum(g_true != 0)), replace = F)
  g_known[sub1] <- g_true[sub1]
  g_known[sub2] <- g_true[sub2]
  supervisedSample <- c(sub1, sub2)
  
  # initialization
  g_init <- rdirichlet(n = nrow(X), alpha = c(10,0.5))[,2]
  est <- LS_partial_est(X = dat, g_init = g_init, 
                        g_known = g_known, mean_pure = mu_p)
  
  DeltaGrid <- seq(from = 0, to = 25, length.out = 8)
  MembGrid <- seq(from = 0, to = 100, length.out = 8)
  OmegaGrid <- c(0.001, 0.01, 0.1, 1)
  
  cond <- F
  DeltaOld <- 1
  MembOld <- 1
  OmegaOld <- 1
  
  while(cond == FALSE) {
    
    mod0 <- penPMM(data = X, Dmat = D, penalty_lasso = DeltaOld,
                   penalty_graph = OmegaGrid, penalty_memb = MembOld, 
                   rho = 0.7, 
                   init_memb = est$membership, 
                   init_prec = solve(cov(X) + diag(0.001, p)), 
                   init_shift = est$meanShift, 
                   supervision = g_known, verbose = F, tol_llik = 0.01)
    models_BICS <- sapply(mod0, "[[", "BIC")
    max_bic_model <- which.max(models_BICS)
    opt0 <- mod0[[max_bic_model]]
    OmegaOpt <- as.numeric(opt0$hyperpar["penalty_graph"])
    cat("\n OmegaOpt = ", OmegaOpt)
    
    mod1 <- penPMM(data = X, Dmat = D, penalty_lasso = DeltaGrid,
                   penalty_graph = OmegaOpt, penalty_memb = MembOld, 
                   rho = 0.7, 
                   init_memb = est$membership, 
                   init_prec = solve(cov(X) + diag(0.001, p)), 
                   init_shift = est$meanShift, 
                   supervision = g_known, verbose = F, tol_llik = 0.01)
    models_BICS <- sapply(mod1, "[[", "BIC")
    max_bic_model <- which.max(models_BICS)
    opt1 <- mod1[[max_bic_model]]
    DeltaOpt <- as.numeric(opt1$hyperpar["penalty_lasso"])
    cat("\n DeltaOpt = ", DeltaOpt)
    
    mod2 <- penPMM(data = X, Dmat = D, penalty_lasso = DeltaOpt,
                   penalty_graph = OmegaOpt, penalty_memb = MembGrid, 
                   rho = 0.7, 
                   init_memb = est$membership + runif(n, min = -0.01, max = 0.01), 
                   init_prec = solve(cov(X) + diag(0.001, p)), 
                   init_shift = est$meanShift + runif(p, min = -0.1, max = 0.1), 
                   supervision = g_known, verbose = F, tol_llik = 0.01)
    models_BICS <- sapply(mod2, "[[", "BIC")
    max_bic_model <- which.max(models_BICS)
    opt2 <- mod2[[max_bic_model]]
    MembOpt <- as.numeric(opt2$hyperpar["penalty_memb"])
    cat("\n MembOpt = ", MembOpt)
    
    cond <- (MembOpt == MembOld) & 
      (DeltaOpt == DeltaOld) & (OmegaOpt == OmegaOld)
    
    DeltaOld <- DeltaOpt
    MembOld <- MembOpt
    OmegaOld <- OmegaOpt
    
  }
  
  mod <- opt2
  datC <- X - outer(mod$Membership, as.vector(mod$meanShift), FUN = "*")
  mae[ind] <- mean(abs(g_true - mod$Membership))
  tab <- table(mod$Membership == 0, g_true == 0)
  acc[ind] <- sum(diag(tab))/n
  sens[ind] <- tab[2,2]/(tab[2,2] + tab[1,2])
  spec[ind] <- tab[1,1]/(tab[1,1] + tab[2,1])
  mse[ind] <- mean((mod$meanShift - deltaTrue)^2)
  
  tabDelta <- table(deltaTrue == 0, mod$gammaEst[p:length(mod$gammaEst)] == 0)
  if (ncol(tabDelta) != 1) {
    accDelta[ind] <- sum(diag(tabDelta))/p
    sensDelta[ind] <- tabDelta[2,2]/(tabDelta[2,2] + tabDelta[1,2])
    specDelta[ind] <- tabDelta[1,1]/(tabDelta[1,1] + tabDelta[2,1])
  } else {
    accDelta[ind] <- NA
    sensDelta[ind] <- NA
    specDelta[ind] <- NA
  }
  cat("\n", c(ind, date()), sep = " --- ")
  
}


