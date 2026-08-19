source("sparsePMM_EstMod.R")

penPMM <- function(data, 
                   Dmat, 
                   penalty_lasso, 
                   penalty_graph, 
                   penalty_memb, 
                   rho, 
                   init_memb, init_prec, init_shift, 
                   tol_abs = 1e-03, tol_rel = 1e-03, 
                   tol_llik = 1e-05, maxit = 500,
                   maxit_ADMM = 10000,
                   supervision = NULL,
                   adaptive = FALSE, 
                   verbose = TRUE) {
  
  all_hyperpar <- expand.grid(penalty_lasso = penalty_lasso,
                              penalty_graph = penalty_graph, 
                              penalty_memb = penalty_memb)
  n_mod <- nrow(all_hyperpar)
  mod_container <- vector(mode = "list", length = n_mod)
  
  if (verbose) {
    cat("Fitting:\n")
    utils::flush.console()
    pbar <- utils::txtProgressBar(min = 0, max = n_mod, style = 3)
    on.exit(close(pbar))
    ipbar <- 0
  }
  
  for (model in 1:n_mod) {
    
    mod_container[[model]] <- 
      tryCatch(sparsePMM(
        data = data, Dmat = D, 
        penalty_lasso = all_hyperpar[model,"penalty_lasso"],
        penalty_graph = all_hyperpar[model,"penalty_graph"],
        penalty_memb = all_hyperpar[model,"penalty_memb"],
        init_memb = init_memb, 
        init_shift = init_shift, 
        init_prec = init_prec,
        rho = rho,
        maxit = maxit, maxit_ADMM = maxit_ADMM, 
        adaptive = adaptive, 
        tol_llik = tol_llik, supervision = supervision, 
        # tol_primal = tol_primal, tol_dual = tol_dual
        tol_abs = tol_abs, tol_rel = tol_rel
      ), error = function(e) {
        message("Error encountered: ", e$message)
        list(bic = NA)
      })
    
    if (verbose) {
      ipbar <- ipbar + 1
      utils::setTxtProgressBar(pbar, ipbar)
    }
  }
  return(mod_container)
}
