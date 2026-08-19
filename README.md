# Main R Functions

The repository contains three main functions used to initialize and estimate the model.

## `LS_partial_est.R`

This function implements the initialization procedure described in **Section 4.2 of the paper**. 

### Main arguments

* `X`: an (n \times p) matrix containing the observed spectra, with (n) samples measured over (p) wavelengths.
* `g_init`: a vector of length (n) containing the initial values of the membership parameters (g_i).
* `g_known`: a vector of length (n) used for the semi-supervised version of the model. For a fully unsupervised analysis, this can be initialized as a vector of `NA` values. If the membership parameter (g_i) is known for some observations, the corresponding entries can instead be set to their known values.
* `tol`: convergence tolerance for the iterative procedure.
* `mean_pure`: mean spectrum of the pure samples. In our analyses the spectra were centered beforehand, so this argument was typically left as `NULL`. Alternatively, a known mean spectrum can be supplied or estimated from the data.

The estimates obtained from this function can be used to initialize the membership and shift parameters in the full model.

## `penPMM.R`

This function is a wrapper for fitting the model over a grid of hyperparameter values.

It calls the main estimation routine in `sparsePMM_EstMod.R` for the requested combinations of penalty parameters and can therefore be used for hyperparameter selection or sensitivity analyses.

## `sparsePMM_EstMod.R`

This is the main model-fitting function. It estimates the model for a given set of hyperparameter values and contains the core optimization routines, including the **ADMM algorithm** and **graphical lasso** estimation.

### Main arguments

* `data`: the matrix containing the observed spectra.

* `Dmat`: the matrix (D) defined in the paper. This needs to be constructed before fitting the model using the appropriate finite-difference/design matrix construction for the application.

* `penalty_lasso`: lasso penalty parameter.

* `penalty_graph`: graphical-lasso penalty parameter.

* `penalty_memb`: penalty parameter associated with the membership component.

  The penalty arguments can be supplied as either scalar values or vectors, depending on whether a single model or a grid of models is being fitted.

* `rho`: initial value of the ADMM parameter (\rho).

* `init_memb`: initial values of the membership parameters (g_i). In our analyses, these were generally obtained from `LS_partial_est.R`.

* `init_prec`: initial estimate of the precision matrix. A sample-based estimate can be used as a starting point.

* `init_shift`: initial values of the shift parameter (\delta). In our analyses, these were also obtained from `LS_partial_est.R`.

* `supervision`: vector containing known membership values for the semi-supervised version of the model. It follows the same convention as `g_known` in `LS_partial_est.R`: unknown entries should be set to `NA`, while known (g_i) values should be entered explicitly.

## Typical workflow

A typical analysis proceeds as follows:

1. Preprocess and, if appropriate, center the spectral data.
2. Construct the matrix (D) required by the model.
3. Use `LS_partial_est.R` to obtain initial estimates for the membership parameters and shifts.
4. Choose an initial precision matrix.
5. Fit the model using `sparsePMM_EstMod.R` for a particular combination of penalty parameters.
6. Alternatively, use `penPMM.R` to fit the model over a grid of hyperparameter values.
