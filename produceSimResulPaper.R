# Read the B replications saved for one scenario --------------------------
code_dir <- "."  # Set the working directory to this folder.
source(file.path(code_dir, "scenarios.R"))

scenario_id <- 1
setting <- scenarios[scenarios$ID == scenario_id, ]
results_file <- file.path(code_dir, "results",
                         paste0("res_", scenario_id, "_", setting$name, ".Rdata"))
loaded_objects <- load(results_file)

cat("Completed:", sum(completed), "of", B, "\n")
cat("Valid PMM:", sum(validPMM), "- valid PLS-DA:", sum(validPLS), "\n")
if (any(completed & (!validPMM | !validPLS))) {
  print(data.frame(replication = which(completed & (!validPMM | !validPLS)),
                   PMM = errorsPMM[completed & (!validPMM | !validPLS)],
                   PLS = errorsPLS[completed & (!validPMM | !validPLS)]))
}

# Paired comparison: use exactly the same successful repetitions ----------
ok <- completed & validPMM & validPLS
cat("Paired repetitions:", sum(ok), "\n")
average <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_
}

comparison <- data.frame(
  metric = c("Accuracy", "Pure correctly classified (sens)",
             "Adulterated detected (spec)", "False positive rate", "AUC"),
  PMM = c(average(acc[ok]), average(sens[ok]), average(spec[ok]),
          average(1 - sens[ok]), average(aucPMM[ok])),
  PLSDA = c(average(accPLS[ok]), average(sensPLS[ok]), average(specPLS[ok]),
            average(1 - sensPLS[ok]), average(aucPLS[ok])),
  MCSE_difference = c(mcse(acc[ok] - accPLS[ok]), mcse(sens[ok] - sensPLS[ok]),
                      mcse(spec[ok] - specPLS[ok]), mcse(sensPLS[ok] - sens[ok]),
                      mcse(aucPMM[ok] - aucPLS[ok]))
)
print(comparison)

# Quantification and signal recovery, as in the paper --------------------
okPMM <- completed & validPMM
estimation <- c(MAE_x100 = average(mae[okPMM]) * 100,
                RMSE_x100 = average(rmse[okPMM]) * 100,
                bias_x100 = average(bias[okPMM]) * 100,
                MSE_delta_x100 = average(mse[okPMM]) * 100,
                accuracy_delta = average(accDelta[okPMM]),
                sensitivity_zeros_delta = average(sensDelta[okPMM]),
                specificity_nonzeros_delta = average(specDelta[okPMM]),
                MSE_Omega = average(mseOmega[okPMM]))
print(estimation)
cat("MCSE of MAE x100:", mcse(mae[okPMM]) * 100, "\n")

# Detection probability as a function of the true adulteration level ------
# No supervised samples contribute to these calculations.
by_concentration <- NULL
for (level in sort(unique(as.vector(gTrue)))) {
  detectionPMM <- detectionPLS <- errorPMM <- rep(NA_real_, B)
  for (ind in which(ok)) {
    ix <- which(is.na(gKnown[,ind]) & gTrue[,ind] == level)
    if (length(ix) > 0) {
      detectionPMM[ind] <- mean(gHat[ix,ind] > zero_tol)
      detectionPLS[ind] <- mean(predPLS[ix,ind] == 1)
      errorPMM[ind] <- mean(abs(gHat[ix,ind] - level))
    }
  }
  by_concentration <- rbind(by_concentration, data.frame(
    concentration = level, repetitions = sum(is.finite(detectionPMM)),
    detection_PMM = average(detectionPMM), detection_PLSDA = average(detectionPLS),
    MCSE_difference = mcse(detectionPMM - detectionPLS), MAE_PMM = average(errorPMM)))
}
print(by_concentration)
# At concentration 0, detection means a false positive.
# Results remain in the workspace; this script does not create extra files.
