# A Latent Variable Model for Identifying and Characterizing Food Adulteration

R code for the simulated-data experiments associated with the revised manuscript by Alessandro Casa, Thomas Brendan Murphy and Michael Fop.

The scripts generate synthetic spectra, fit the proposed method, compare classification performance with PLS-DA, and summarize estimation and variable-selection results. No external datasets are required. This repository contains the simulation workflow only.

## Files

| File | Purpose |
| --- | --- |
| [run_simulations.R](run_simulations.R) | Run one scenario with a standard Monte Carlo loop. |
| [scenarios.R](scenarios.R) | Define the 19 configurations and their signal and covariance structures. |
| [sparsePMM_EstMod.R](sparsePMM_EstMod.R) | Initialize and fit the model, including the ADMM and graphical-lasso updates. |
| [penPMM.R](penPMM.R) | Fit candidate penalty combinations and compute classification metrics. |
| [produceSimResulPaper.R](produceSimResulPaper.R) | Summarize the saved results for one scenario. |

The historical abbreviation `PMM` is retained in function and object names. Initialization is now implemented by `initialize_pmm()` inside `sparsePMM_EstMod.R`; the old `LS_partial_est.R` file is not needed.

## Running the simulations

Download or clone the repository and set the R working directory to the folder containing these files. Install the dependencies if needed:

```r
install.packages(c("glasso", "caret", "pls"))
```

Open `run_simulations.R` and choose the configuration and number of repetitions near the top:

```r
scenario_id <- 1
B <- 100
resume <- FALSE
```

The penalty grids and numerical controls are directly below these settings. Run the script:

```r
source("run_simulations.R")
```

Both methods are fitted within the same replication loop. The proposed method uses sequential BIC tuning. PLS-DA is fitted using `caret::train(method = "pls")`, with the number of components selected by cross-validated classification accuracy on the labelled subset.

Set the same `scenario_id` in `produceSimResulPaper.R`, then run:

```r
source("produceSimResulPaper.R")
```

This prints the number of completed and valid repetitions, paired classification comparisons, estimation errors, variable-selection metrics and detection rates by adulteration level.

## Simulation configurations

| ID | Configuration |
| --- | --- |
| 1, 2, 3 | Scenario 1: n = 100; p = 50, 100, 150, respectively. |
| 4, 5, 6 | Scenario 1: n = 250; p = 50, 100, 150, respectively. |
| 7, 8, 9 | Scenario 1: n = 500; p = 50, 100, 150, respectively. |
| 10, 11 | Scenario 2: weak/strong Mexican-hat signal. |
| 12, 13 | Scenario 2: weak/strong spline signal. |
| 14, 15 | Scenario 3: weak signal with 5%/20% adulteration. |
| 16, 17 | Scenario 3: strong signal with 5%/20% adulteration. |
| 18 | Gaussian errors with AR(1) covariance and correlation parameter 0.9. |
| 19 | Independent Student t errors with 5 degrees of freedom. |

IDs 10–19 use n = 250 and p = 100. ID 5 is the independent Gaussian reference for IDs 18 and 19. All settings have marginal error variance 0.01. The spline profiles are fixed across replications.

## Saved results and reproducibility

Each scenario produces one `.Rdata` file under `results/`, updated after every replication. For example:

```text
results/res_1_S1_n100_p50.Rdata
```

The file includes estimates, true parameters, labels, metrics, tuning histories, convergence flags, fixed seeds, numerical controls, session information and copies of the source code used. Separate seeds control data generation, label selection and PLS cross-validation.

An existing result file is not overwritten by a new run. To continue an interrupted experiment with the same settings, set `resume <- TRUE`. The saved settings are restored and completed replications are skipped. For a different grid or configuration, start a new experiment with a different output filename.

The base simulation runs used R 4.4.2, `glasso` 1.11, `caret` 7.0-1 and `pls` 2.9-0. The installation command above installs the available package versions rather than a frozen environment.

This compact release contains the simulation design and base analysis. Some final manuscript results also used expanded tuning grids and numerical recovery steps, which are not included here; exact agreement with every reported table entry is therefore not implied. Grids and controls can be edited directly in the runner.

## Metrics and convergence

Classification and membership estimation errors are evaluated only on unlabelled observations. PLS-DA uses the binary labels of the same supervised samples; the proposed method uses their known adulteration fractions and also the unlabelled spectra during fitting.

Pure/zero is the positive class:

- Sensitivity: proportion of truly pure samples correctly classified as pure.
- Specificity: proportion of truly adulterated samples correctly classified as adulterated.

The same zero/nonzero convention is used for variable selection. `scorePLS` contains class probabilities, not estimates of adulterant concentration. MAE and MSE are saved in their original units and multiplied by 100 by the summary script; MAE times 100 is in percentage points.

The runner retains the following conservative acceptance rule:

```r
validPMM[ind] <- mod$converged && fitPMM$stable && all(history$converged)
```

Thus, a repetition can be excluded even when its selected model converges, if another candidate did not converge. The selected model, when available, is still saved, but its metrics remain `NA`. Inspect `errorsPMM` and `tuningHistory` when assessing excluded repetitions. The rule can be edited before running a new experiment; changing a saved flag alone does not recompute missing metrics. The summary compares classifiers on the same valid repetitions and reports their count.

BIC is maximized and defined as `2 * logLik - df * log(n)`. Precision entries are counted once, fixed memberships are excluded from the free-parameter count, and mean-shift degrees of freedom count nonzero fused blocks. `OmegaGrid` is on the graphical-lasso scale, with likelihood-scale penalty `n * rho_omega / 2`. The precision penalty is uniform, including diagonal entries. Tuning proceeds over Omega, delta and memberships one coordinate grid at a time.
