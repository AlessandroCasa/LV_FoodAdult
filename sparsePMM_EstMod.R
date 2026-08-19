sparsePMM <- function(data, Dmat,
                      init_memb, init_shift, init_prec,
                      maxit = 1000, maxit_ADMM = 1000,
                      rho, 
                      penalty_lasso, 
                      penalty_graph, 
                      penalty_memb,
                      adaptive = FALSE, 
                      tol_llik = 0.01, 
                      tol_rel = 0.001, tol_abs = 0.01, 
                      supervision, verbose = TRUE) {
  
  iter_Ext <- 1
  llikDiff <- llikOld <- Inf
  llikHist <- llik_noPen <- vector(length = maxit)
  
  nDmat <- nrow(Dmat)
  p <- ncol(data)
  n <- nrow(data)
  
  gEst <- init_memb
  deltaEst <- init_shift
  OmegaEst <- init_prec
  
  while(llikDiff > tol_llik) {
    iter_Ext <- iter_Ext + 1
    
    tDmat <- t(Dmat)
    a <- solve(sum(gEst^2)*OmegaEst + rho*crossprod(Dmat)) 
    OmDat <- data%*%OmegaEst
    sumD <- apply(gEst*OmDat, 2, sum) 
    iter_ADMM <- 0
    primalRes <- dualRes <- Inf
    tol_primal <- tol_dual <- 0
    while ((primalRes > tol_primal | dualRes > tol_dual) & 
           iter_ADMM < maxit_ADMM) {
      
      iter_ADMM <- iter_ADMM + 1  
      
      if (iter_ADMM == 1) {
        gamma <- runif(n = nDmat)
        gammaNew <- vector(length = nDmat)
        u <- rnorm(n = nDmat)
      }
      aa <- sumD + rho*(tDmat%*%gamma) - tDmat%*%u
      deltaEst <- a%*%aa

      Dstar <- Dmat%*%deltaEst
      
      condEl <- u + rho*Dstar
      gammaNew <- ifelse(abs(condEl) < penalty_lasso,
                         0,
                         sign(condEl)*(abs(condEl) - penalty_lasso))
      
      
      gammaNew <- gammaNew/rho
    
      dualRes <- rho*(tDmat%*%diag(-1, nrow = nDmat)%*%(gammaNew - gamma))
      dualRes <- sqrt(sum(dualRes^2))
      gamma <- gammaNew
      
      primalRes <- Dmat%*%deltaEst - gamma

      primalRes <- sqrt(sum(primalRes^2))
      
      unew <- u + rho*(Dstar - gamma)
      u <- unew
      
      selMax <- max(sqrt(sum((Dmat%*%deltaEst)^2)), 
                    sqrt(sum((diag(-1, nrow = nDmat)%*%gamma)^2)),
                    0
      )
      tol_primal <- sqrt(nDmat)*tol_abs + tol_rel*selMax
      tol_dual <- sqrt(p)*tol_abs + 
        tol_rel*sqrt(sum((tDmat%*%u)^2))
      
      if (adaptive == TRUE) {
        rho <- ifelse(primalRes > 5*dualRes,
                      2*rho,
                      ifelse(dualRes > 5*primalRes,
                             rho/2,
                             rho))
      }
    }
    
    if(verbose)  cat("\n ADMM iter = ", iter_ADMM, 
                     "- DONE, rho = ",rho, "\n ------ \n")
    
    denomMemb <- t(deltaEst)%*%OmegaEst%*%deltaEst
    nomMemb <- t(deltaEst)%*%OmegaEst%*%t(data)
    gEstNew <- ifelse(abs(nomMemb) < penalty_memb, 0,
                      sign(nomMemb)*(abs(nomMemb) - penalty_memb))
    gEstNew <- gEstNew/as.numeric(denomMemb)
    if (!is.null(supervision)) {
      gEstNew[!is.na(supervision)] <- supervision[!is.na(supervision)]
    }
    gEst <- ifelse(gEstNew < 0, 0, ifelse(gEstNew > 0.5, 0.5, gEstNew))
    gEst <- as.vector(gEst)
    
    datC <- data - outer(gEst, as.vector(deltaEst), FUN = "*")
    
    Smat <- cov(datC)
    matRho <- 1/abs(solve(cov(datC)))
    grapLas <- glasso(s = Smat, rho = matRho*penalty_graph)

    OmegaEst <- grapLas$wi
    OmegaEst <- (OmegaEst + t(OmegaEst))/2
    
    if(verbose) cat("\n glasso iter = ", iter_Ext, 
                    "- DONE","\n ---------- \n")
    
    llik <- vector(length = n)
    for (ll in 1:n) {
      llik[ll] <- LaplacesDemon::dmvnp(x = data[ll,], 
                                       mu = as.vector(gEst[ll]*deltaEst), 
                                       Omega = as.positive.definite(OmegaEst), 
                                       log = TRUE) 
    }
    llik_noPen[iter_Ext] <- sum(llik)
    llikHist[iter_Ext] <- llik_noPen[iter_Ext] - 
      penalty_lasso*sum(abs(Dmat%*%deltaEst)) - 
      penalty_graph*sum(abs(OmegaEst)) - 
      penalty_memb*sum(abs(gEst))
    llikDiff <- abs(llikHist[iter_Ext - 1] - llikHist[iter_Ext])/abs(llikHist[iter_Ext - 1])
    if (verbose) print(llikDiff)
  }
  
  # counting parameters for BIC 
  gamMat <- cbind(gamma[p:nDmat], c(NA, gamma[1:(p-1)]))
  sumPar <- 0
  for (ii in 2:nrow(gamMat)) {
    if (gamMat[ii,2] == 0 & 
        gamMat[ii-1,1] != 0 & 
        gamMat[ii,1] != 0) {
      sumPar = sumPar + 1
    }
  }
  paramGamma <- p - sum(gamMat[,1] == 0) - sumPar
  nparam <- paramGamma + sum(OmegaEst != 0) + sum(gEst != 0)
  
  OUT <- 
    list(
      meanShift = deltaEst, 
      nParam = nparam, 
      membership = gEst,
      precision = OmegaEst,
      loglik = llikHist[iter_Ext],
      BIC = 2*llik_noPen[iter_Ext] - nparam*log(n),
      nIter = iter_Ext, 
      llik_trace = llikHist[2:iter_Ext],
      hyperparameters = c(penalty_lasso, penalty_graph, penalty_memb),
      gamma = gamma,
      llik_noPen = llik_noPen[iter_Ext]
    )
  
  return(OUT)
}

