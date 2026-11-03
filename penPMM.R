# Fit one or more triples of penalties, as in the original penPMM function.
source(file.path(code_dir, "sparsePMM_EstMod.R"))

penPMM <- function(data, penalty_lasso, penalty_graph, penalty_memb,
                   supervision, initial, control, upper_g = 0.5) {
  all_hyperpar <- expand.grid(penalty_lasso = penalty_lasso,
                             penalty_graph = penalty_graph,
                             penalty_memb = penalty_memb)
  mod_container <- vector("list", nrow(all_hyperpar))

  for (model in 1:nrow(all_hyperpar)) {
    mod_container[[model]] <- tryCatch(
      sparsePMM(X = data, g_known = supervision,
                lambda_delta = all_hyperpar$penalty_lasso[model],
                rho_omega = all_hyperpar$penalty_graph[model],
                lambda_g = all_hyperpar$penalty_memb[model],
                control = control, upper_g = upper_g, initial = initial),
      error = function(e) list(BIC = NA_real_, converged = FALSE,
                               reason = conditionMessage(e))
    )
    # Keep the model and its diagnostic, but do not select a failed fit by BIC.
    if (!mod_container[[model]]$converged) mod_container[[model]]$BIC <- NA_real_
  }
  mod_container
}

# Small table of the candidate fits; stored for each Monte Carlo replication.
candidate_table <- function(models, DeltaGrid, OmegaGrid, MembGrid, step) {
  tab <- expand.grid(lambda_delta = DeltaGrid, rho_omega = OmegaGrid, lambda_g = MembGrid)
  tab$step <- step
  tab$BIC <- sapply(models, "[[", "BIC")
  tab$converged <- sapply(models, "[[", "converged")
  tab$reason <- sapply(models, "[[", "reason")
  tab
}

# Fixed 2 x 2 table, even when every prediction is pure/zero.
# Keep the paper's convention: sens = recovery of zeros, spec = recovery of nonzeros.
classification_metrics <- function(predicted_zero, true_zero) {
  tab <- table(factor(predicted_zero, levels = c(FALSE, TRUE)),
               factor(true_zero, levels = c(FALSE, TRUE)))
  sens <- if (sum(tab[,2]) > 0) tab[2,2] / sum(tab[,2]) else NA_real_
  spec <- if (sum(tab[,1]) > 0) tab[1,1] / sum(tab[,1]) else NA_real_
  c(acc = sum(diag(tab)) / sum(tab), sens = sens, spec = spec)
}

auc_score <- function(score, adulterated) {
  n1 <- sum(adulterated)
  n0 <- sum(!adulterated)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(score)[adulterated]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
