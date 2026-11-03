# Simulations: one scenario, B Monte Carlo replications -------------------
# Working directory: the folder containing these five R scripts.

code_dir <- "."  # Set the working directory to this folder.
if (dir.exists(file.path(code_dir, ".Rlibrary"))) {
  .libPaths(c(file.path(code_dir, ".Rlibrary"), .libPaths()))
}
library(glasso)
library(caret)
library(pls)
source(file.path(code_dir, "scenarios.R"))
source(file.path(code_dir, "penPMM.R"))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# Choose the scenario and number of repetitions ---------------------------

scenario_id <- 1
B <- 100
resume <- FALSE

print(scenarios)
scenario <- define_scenario(scenario_id)
n <- scenario$n
p <- scenario$p
mu_p <- scenario$mu_p
deltaTrue <- scenario$deltaTrue
SigmaTrue <- scenario$SigmaTrue
OmegaTrue <- scenario$OmegaTrue
g_levels <- scenario$g_levels

prob_noAd <- 0.7
percSup <- 0.2
upper_g <- 0.5
zero_tol <- 1e-8

# Tuning grids: OmegaGrid contains the scalar passed to glasso.
# In the sum log-likelihood its penalty is lambda_omega = n * rho_omega / 2.
DeltaGrid <- c(0, 0.1, 0.5, 1, 2.5, 5, 10, 25)
OmegaGrid <- c(0.00001, 0.0001, 0.0005, 0.001, 0.005, 0.01)
MembGrid <- c(0, 0.5, 1, 2.5, 5, 10, 25, 50, 100)
max_cycles <- 6

control <- list(maxit = 1000, outer_tol = 1e-6, parameter_tol = 1e-4,
                admm_maxit = 5000, admm_abs = 1e-7, admm_rel = 1e-6,
                glasso_thr = 1e-7, glasso_maxit = 10000,
                zero_tol = zero_tol, objective_tol = 1e-7)

results_dir <- file.path(code_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
results_file <- file.path(results_dir,
                         paste0("res_", scenario_id, "_", scenario$setting$name, ".Rdata"))

# Objects to fill --------------------------------------------------------

mae <- rmse <- bias <- acc <- sens <- spec <- mse <- rep(NA_real_, B)
accDelta <- sensDelta <- specDelta <- mseOmega <- rep(NA_real_, B)
accOmega <- sensOmega <- specOmega <- rep(NA_real_, B)
accPLS <- sensPLS <- specPLS <- aucPMM <- aucPLS <- rep(NA_real_, B)
maePure <- maeAdulterated <- rep(NA_real_, B)
ncompPLS <- rep(NA_integer_, B)

# Columns correspond to replications. Known samples are excluded from metrics.
gTrue <- gKnown <- gHat <- predPLS <- scorePLS <- matrix(NA_real_, n, B)
deltaHat <- matrix(NA_real_, p, B)
models <- tuningHistory <- plsDetails <- vector("list", B)
validPMM <- validPLS <- completed <- rep(FALSE, B)
errorsPMM <- errorsPLS <- rep(NA_character_, B)
timePMM <- timePLS <- rep(NA_real_, B)

# Different seeds for generation, supervision and PLS cross-validation.
# They do not depend on fitting times, convergence or the number of candidates.
seeds <- data.frame(membership = 100000 + 1:B, noise = 200000 + 1:B,
                    supervision = 300000 + 1:B, pls_cv = 400000 + 1:B)
session_info <- sessionInfo()
rng_kind <- RNGkind()
code_used <- lapply(file.path(code_dir, c("run_simulations.R", "scenarios.R",
                                         "penPMM.R", "sparsePMM_EstMod.R")), readLines)

if (resume) {
  # Load previous settings and accumulated vectors from this scenario.
  # Completed repetitions are skipped in the for loop below.
  load(results_file)
} else if (file.exists(results_file)) {
  stop("This result file already exists. Use resume <- TRUE, or change results_file.")
}

# Monte Carlo loop -------------------------------------------------------

for (ind in 1:B) {
  if (completed[ind]) next
  cat("\nScenario", scenario_id, "- repetition", ind, "-", date(), "\n")

  # Generate adulteration levels and spectra -----------------------------
  set.seed(seeds$membership[ind])
  g_true <- g_levels[sample.int(length(g_levels), size = n, replace = TRUE)]
  g_true <- g_true * rbinom(n, size = 1, prob = 1 - prob_noAd)

  set.seed(seeds$noise[ind])
  if (scenario$setting$error == "student") {
    # Independent t errors at every wavelength, each with variance sigma^2.
    noise <- matrix(rt(n * p, df = 5), nrow = n, ncol = p)
    noise <- noise * scenario$sigma * sqrt(3 / 5)
  } else {
    noise <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% chol(SigmaTrue)
  }
  meanMat <- matrix(mu_p, nrow = n, ncol = p, byrow = TRUE) + outer(g_true, deltaTrue)
  dat <- meanMat + noise
  X <- sweep(dat, 2, mu_p, "-")

  # Same labelled samples for PMM and PLS-DA ------------------------------
  set.seed(seeds$supervision[ind])
  nPureKnown <- round(percSup * sum(g_true == 0))
  nAdKnown <- round(percSup * sum(g_true != 0))
  if (min(nPureKnown, nAdKnown) < 2) {
    stop("Too few labelled samples in one class for PLS cross-validation.")
  }
  sub1 <- sample(which(g_true == 0), size = nPureKnown, replace = FALSE)
  sub2 <- sample(which(g_true != 0), size = nAdKnown, replace = FALSE)
  g_known <- rep(NA_real_, n)
  g_known[sub1] <- g_true[sub1]
  g_known[sub2] <- g_true[sub2]
  supervisedSample <- c(sub1, sub2)
  testSample <- which(is.na(g_known))
  gTrue[,ind] <- g_true
  gKnown[,ind] <- g_known

  # PMM: initialization and sequential BIC tuning -------------------------
  start_time <- proc.time()[3]
  history <- NULL
  fitPMM <- tryCatch({
    initial <- initialize_pmm(X, g_known, upper_g)
    DeltaOld <- 1
    OmegaOld <- 0.001
    MembOld <- 1
    cond <- FALSE
    cycle <- 0

    while (!cond && cycle < max_cycles) {
      cycle <- cycle + 1

      mod0 <- penPMM(data = X, penalty_lasso = DeltaOld,
                     penalty_graph = OmegaGrid, penalty_memb = MembOld,
                     supervision = g_known, initial = initial, control = control)
      history <- rbind(history, candidate_table(mod0, DeltaOld, OmegaGrid, MembOld, "Omega"))
      models_BICS <- sapply(mod0, "[[", "BIC")
      if (all(!is.finite(models_BICS))) stop("No converged candidate for Omega.")
      opt0 <- mod0[[which.max(models_BICS)]]
      OmegaOpt <- unname(opt0$hyperparameters["rho_omega"])

      mod1 <- penPMM(data = X, penalty_lasso = DeltaGrid,
                     penalty_graph = OmegaOpt, penalty_memb = MembOld,
                     supervision = g_known, initial = initial, control = control)
      history <- rbind(history, candidate_table(mod1, DeltaGrid, OmegaOpt, MembOld, "Delta"))
      models_BICS <- sapply(mod1, "[[", "BIC")
      if (all(!is.finite(models_BICS))) stop("No converged candidate for delta.")
      opt1 <- mod1[[which.max(models_BICS)]]
      DeltaOpt <- unname(opt1$hyperparameters["lambda_delta"])

      mod2 <- penPMM(data = X, penalty_lasso = DeltaOpt,
                     penalty_graph = OmegaOpt, penalty_memb = MembGrid,
                     supervision = g_known, initial = initial, control = control)
      history <- rbind(history, candidate_table(mod2, DeltaOpt, OmegaOpt, MembGrid, "Membership"))
      models_BICS <- sapply(mod2, "[[", "BIC")
      if (all(!is.finite(models_BICS))) stop("No converged candidate for memberships.")
      opt2 <- mod2[[which.max(models_BICS)]]
      MembOpt <- unname(opt2$hyperparameters["lambda_g"])

      cond <- (DeltaOpt == DeltaOld) && (OmegaOpt == OmegaOld) && (MembOpt == MembOld)
      DeltaOld <- DeltaOpt
      OmegaOld <- OmegaOpt
      MembOld <- MembOpt
      cat("Cycle", cycle, "- Delta:", DeltaOpt, "Omega:", OmegaOpt, "Memb:", MembOpt, "\n")
    }
    list(model = opt2, stable = cond)
  }, error = function(e) e)

  timePMM[ind] <- proc.time()[3] - start_time
  tuningHistory[[ind]] <- history
  if (inherits(fitPMM, "error")) {
    errorsPMM[ind] <- conditionMessage(fitPMM)
  } else {
    mod <- fitPMM$model
    models[[ind]] <- mod
    gHat[,ind] <- mod$membership
    deltaHat[,ind] <- mod$meanShift
    validPMM[ind] <- mod$converged && fitPMM$stable && all(history$converged)
    if (!validPMM[ind]) {
      errorsPMM[ind] <- "Check candidate convergence/tuning stability in tuningHistory."
    } else {
      # Evaluate ONLY the unlabelled samples, directly within this replication.
      m <- classification_metrics(mod$membership[testSample] <= zero_tol,
                                  g_true[testSample] == 0)
      acc[ind] <- m["acc"]
      sens[ind] <- m["sens"]  # P(predicted pure | pure), as in the old paper.
      spec[ind] <- m["spec"]  # P(predicted adulterated | adulterated).
      error_g <- mod$membership[testSample] - g_true[testSample]
      mae[ind] <- mean(abs(error_g))
      rmse[ind] <- sqrt(mean(error_g^2))
      bias[ind] <- mean(error_g)
      maePure[ind] <- mean(abs(error_g[g_true[testSample] == 0]))
      maeAdulterated[ind] <- mean(abs(error_g[g_true[testSample] != 0]))
      aucPMM[ind] <- auc_score(mod$membership[testSample], g_true[testSample] != 0)
      mse[ind] <- mean((mod$meanShift - deltaTrue)^2)
      mDelta <- classification_metrics(abs(mod$meanShift) <= zero_tol,
                                       abs(deltaTrue) <= zero_tol)
      accDelta[ind] <- mDelta["acc"]
      sensDelta[ind] <- mDelta["sens"]
      specDelta[ind] <- mDelta["spec"]
      mseOmega[ind] <- mean((mod$precision - OmegaTrue)^2)
      upper <- upper.tri(OmegaTrue)
      mOmega <- classification_metrics(abs(mod$precision[upper]) <= zero_tol,
                                       abs(OmegaTrue[upper]) <= zero_tol)
      accOmega[ind] <- mOmega["acc"]
      sensOmega[ind] <- mOmega["sens"]
      specOmega[ind] <- mOmega["spec"]
    }
  }

  # PLS-DA, on the same dataset and within the same for loop ---------------
  start_time <- proc.time()[3]
  fitPLS <- tryCatch({
    set.seed(seeds$pls_cv[ind])
    Xtrain <- X[supervisedSample, , drop = FALSE]
    Xtest <- X[testSample, , drop = FALSE]
    colnames(Xtrain) <- colnames(Xtest) <- paste0("V", 1:p)
    classTrain <- factor(ifelse(g_known[supervisedSample] == 0, "Pure", "Adulterated"),
                          levels = c("Pure", "Adulterated"))
    k_folds <- min(5, min(table(classTrain)))
    folds <- caret::createFolds(classTrain, k = k_folds, returnTrain = TRUE)
    # Components must be estimable in every training fold.
    fold_ranks <- sapply(folds, function(i) {
      qr(scale(Xtrain[i, , drop = FALSE], center = TRUE, scale = FALSE))$rank
    })
    ncomp_max <- min(10, p, min(sapply(folds, length)) - 1, min(fold_ranks))
    if (ncomp_max < 1) stop("No estimable PLS component in the training folds.")
    ctrlPLS <- caret::trainControl(method = "cv", number = k_folds, index = folds,
                                   classProbs = TRUE, returnResamp = "all")
    pls_fit <- caret::train(x = Xtrain, y = classTrain, method = "pls",
                            metric = "Accuracy", trControl = ctrlPLS,
                            tuneGrid = data.frame(ncomp = 1:ncomp_max))
    if (any(!is.finite(pls_fit$resample$Accuracy))) stop("Incomplete PLS cross-validation.")
    prediction <- predict(pls_fit, newdata = Xtest)
    probability <- predict(pls_fit, newdata = Xtest, type = "prob")$Adulterated
    if (anyNA(prediction) || any(!is.finite(probability))) stop("Invalid PLS predictions.")
    list(prediction = prediction, probability = probability,
         ncomp = pls_fit$bestTune$ncomp, CVresults = pls_fit$results,
         CVfolds = folds, CVresamples = pls_fit$resample)
  }, error = function(e) e)

  timePLS[ind] <- proc.time()[3] - start_time
  if (inherits(fitPLS, "error")) {
    errorsPLS[ind] <- conditionMessage(fitPLS)
  } else {
    validPLS[ind] <- TRUE
    plsDetails[[ind]] <- fitPLS
    predPLS[testSample,ind] <- as.numeric(fitPLS$prediction == "Adulterated")
    scorePLS[testSample,ind] <- fitPLS$probability
    ncompPLS[ind] <- fitPLS$ncomp
    mPLS <- classification_metrics(fitPLS$prediction == "Pure", g_true[testSample] == 0)
    accPLS[ind] <- mPLS["acc"]
    sensPLS[ind] <- mPLS["sens"]
    specPLS[ind] <- mPLS["spec"]
    aucPLS[ind] <- auc_score(fitPLS$probability, g_true[testSample] != 0)
  }

  # One Rdata file for this scenario, updated after each repetition --------
  completed[ind] <- TRUE
  save(scenario_id, scenario, B, n, p, mu_p, deltaTrue, SigmaTrue, OmegaTrue, g_levels,
       prob_noAd, percSup, upper_g, zero_tol, DeltaGrid, OmegaGrid, MembGrid,
       max_cycles, control, seeds, session_info, rng_kind, code_used,
       mae, rmse, bias, acc, sens, spec, mse, accDelta, sensDelta, specDelta,
       mseOmega, accOmega, sensOmega, specOmega, maePure, maeAdulterated,
       accPLS, sensPLS, specPLS, aucPMM, aucPLS, ncompPLS,
       gTrue, gKnown, gHat, predPLS, scorePLS, deltaHat, models, tuningHistory, plsDetails,
       validPMM, validPLS, errorsPMM, errorsPLS, timePMM, timePLS, completed,
       file = results_file)
  cat("Saved:", results_file, "- PMM:", validPMM[ind], "PLS-DA:", validPLS[ind], "\n")
  if (!validPMM[ind]) cat(errorsPMM[ind], "\n")
  if (!validPLS[ind]) cat(errorsPLS[ind], "\n")
}

cat("\nCompleted:", sum(completed), "of", B,
    "- valid PMM:", sum(validPMM), "- valid PLS-DA:", sum(validPLS), "\n")
# Load the Rdata file in produceSimResulPaper.R to inspect the results.
