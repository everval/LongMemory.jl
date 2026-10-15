# Changelog

All notable changes to this project will be documented in this file.

## [Version 1.0.0] - 2024-06-30

First release of `LongMemory.jl`. This package provides tools for generating, estimating, and forecasting long memory time series models in Julia. The package includes functions for fractional differencing, stochastic error duration, cross-sectional aggregation, and various estimators for the Hurst effect. It also offers plotting capabilities to illustrate long memory dynamics and forecasting. The implementation is designed to be efficient and user-friendly, leveraging the speed of Julia.

## [Version 1.0.1] - 2026-04-14

### New Features

- Added `hitransform` function to `Forecasters.jl`, the harmonic inverse transformation of Hassler and Hosseinkouchack (2020), i.e. the inverse of the harmonically weighted filter `hwfilter`. The filter is applied with the fast Fourier transform (FFT).

- To support the new `hitransform` function, added `gregory_coeffs`, which computes the (absolute values of the) Gregory coefficients by recursion.

### Bug Fixes and Performance Improvements

#### ParametricEstimators.jl
- Fixed inefficient double matrix factorizations in `fi_llk` and `csa_llk` by switching to Cholesky decomposition (significant speedup for MLE optimization)
- Optimized `fi_var_vals` to use in-place multiplication (`.*=` instead of `=`) to reduce allocations
- Optimized `my_toeplitz` construction with single-pass triangle filling and `@inbounds` annotation
- Improved `fi_mle_est` and `csa_mle_est` σ estimation using Cholesky decomposition
- Optimized `har_est` inner loop accumulation with `@view` and in-place operations to eliminate temporary allocations

#### Forecasters.jl
- Fixed `hitransform` to not mutate input array; now returns a new array with results

## [Version 2.0.0] - Unreleased

This is a breaking release. The main changes for users are that plotting requires loading Plots.jl, that functions return vectors instead of `T×1` matrices, and that the data sets are returned as named tuples. See "Breaking changes" below.

### Breaking changes
- **Plotting moved to a package extension.** Plots.jl is no longer installed with LongMemory.jl. The plotting functions (`autocorrelation_plot`, `log_variance_plot`, `rescaled_range_plot`, `dfa_plot`, `periodogram_plot`, `fi_forecast_plot`, `csa_forecast_plot`, `har_forecast_plot`, `arfima_forecast_plot`, `NileDataPlot`, `NHTempDataPlot`, `LMPlot`) are available after `using Plots`; install it once with `import Pkg; Pkg.add("Plots")`. Without Plots, these functions throw an error explaining how to enable them. Loading LongMemory.jl is about 3 times faster (0.5 s instead of 1.65 s) and the installation is much lighter.
- **Vectors instead of `T×1` matrices.** Functions that returned a single column as a `T×1` matrix now return a vector, also when the input is a `T×1` matrix. This affects `fracdiff` (with `d = 0` or `d = 1`), `csa_gen(T, N, p, q)`, `edm_gen`, `sds_gen`, `fi_survival_probs`, `autocovariance`, `autocorrelation`, `rescaled_range`, `fi_forecast`, `csa_forecast` and `arfima_forecast` (without `σ`), and `hitransform`. Indexing with `x[i]` or `x[i, 1]` works as before, but `size(x)` changes and `x == ones(T, 1)` is now `false`; compare with `ones(T)` instead. Outputs with several columns, such as forecasts with confidence bands, are still matrices.
- **Data sets as named tuples.** `NileData()` and `NHTempData()` return named tuples of vectors instead of `DataFrame`s, so CSV.jl and DataFrames.jl are no longer dependencies. Columns are accessed as before, e.g. `NileData().NileMin`, and `DataFrame(NileData())` recovers the `DataFrame`.
- **Removed:** `fi` (use `fi_gen`, which is identical), the keyword `detrend` of `exact_whittle_est` (use `demean`), and the keyword `flag` of `autocorrelation` (use `autocorrelation_plot`). `my_toeplitz` is no longer exported.
- **Forecast confidence bands** of `fi_forecast`, `csa_forecast` and `har_forecast` now widen with the horizon: they are ± 2 standard errors of the `k`-step forecast error, computed from the MA representation of the model. Previously they were ± 2σ at every horizon. `σ` is the standard deviation of the innovations.
- `fi_forecast` and `csa_forecast` demean the series and add the sample mean back to the forecasts (the plotting functions already demeaned).
- Several bug fixes change numerical output, see "Bug Fixes" below.

### New Features

#### Estimation and testing
- `arfima_whittle_est`: parametric Whittle estimation of ARFIMA(p,d,q) models, imposing stationarity and invertibility.
- `feasible_exact_whittle_est`: the feasible exact local Whittle estimator of Shimotsu (2010) for series with unknown mean, valid also for `d > 1/2`.
- `durbin_levinson`: one-step predictions and prediction error variances of a stationary series. `fi_llk`, `csa_llk`, `fi_mle_est` and `csa_mle_est` now use it, so the likelihood costs `O(T²)` instead of `O(T³)`, with the same values (e.g. `fi_mle_est` with `T = 2000` takes 0.06 s instead of 1.2 s).
- `dfa_est`, `dfa_plot`, `dfa_fluctuation`: detrended fluctuation analysis.
- `lo_rs_test`: Lo's (1991) modified rescaled range test, with an Andrews (1991) data-dependent bandwidth.
- `vs_test`: the V/S test of Giraitis et al. (2003).
- `qu_test`: Qu's (2011) test against spurious long memory.
- `lm_change_test(y; flag = true)` also returns the test statistic and a parametric bootstrap p-value.

#### Generation
- `davies_harte`: exact simulation of a stationary Gaussian series with given autocovariances.
- `fi_exact_gen`, `csa_exact_gen`: exact (stationary) simulation of the FI and CSA processes.

#### Forecasting
- `arfima_forecast`, `arfima_forecast_plot`, `arfima_ar_coefs`, `arfima_ma_coefs`: forecasting with ARFIMA(p,d,q) models.
- `fi_ma_coefs`: MA coefficients of the fractionally differenced process.

### Other changes
- `csa_mle_est` uses deterministic starting values (`p = 2`, `q = 2(1-d)` with `d` from `gph_est`) instead of random ones, so the estimates are reproducible.
- `log_variance_est` and `log_variance_plot` throw an error when fewer than 3 points enter the regression.
- All functions accept any `AbstractArray` (ranges, views, ...) instead of only `Array`.
- Dependencies: Plots.jl is now a weak dependency (package extension `LongMemoryPlotsExt`); CSV.jl, DataFrames.jl and Random were removed; DelimitedFiles was added.

### Bug Fixes

Some of these fixes change numerical output.

#### GeneratingFunctions.jl
- `arfima_gen`: the MA filter weighted the `j`-th MA coefficient by `q+2-j` (e.g. MA(1) gave `2e_t + θe_{t-1}`); it now computes `e_t + Σ θ_j e_{t-j}`.
- `arfima_gen`, `arfi_gen`, `fi_gen`, `csa_gen`, `edm_gen`, `sds_gen`: `μ` is now the mean of the generated series, as documented. Previously it was the innovation mean (so the series mean grew with `t` for `d > 0`), an AR intercept, or ignored entirely (`arfima_gen` without AR terms).
- `edm_gen`/`sds_gen`: return `T` observations; the pre-sample was not discarded (length was `T + t`).
- `fracdiff(x, d::Int)`: integer `d` other than 0 and 1 threw `UndefVarError` (e.g. `fi_gen(T, 1)`); it now falls back to the truncated filter. `fracdiff` accepts any `Real` `d`.
- `hwp_gen`: removed a stray random draw.

#### LogPeriodEstimators.jl
- `gph_est`, `whittle_est`, `exact_whittle_est`: use Fourier frequencies `j = 1, ..., round(T^m)` (previously `j = 1, ..., round(T^m) - 1`), consistent with the variance formulas, and cap the band at the Nyquist frequency instead of throwing a `BoundsError`.

#### ClassicEstimators.jl
- `rescaled_range` (and `rescaled_range_est`/`rescaled_range_plot`): the R/S statistic of `x[1:n]` now uses the subsample mean instead of the full-sample mean.
- `sstdk` (and `log_variance_est`/`log_variance_plot`): the last complete block was dropped.
- `autocorrelation_plot`: the lag axis now starts at 0.

#### ParametricEstimators.jl
- `csa_var_vals`, `csa_var_matrix`: removed a spurious factor `B(p,q)`; the autocovariance is `B(p + k/2, q - 1)/B(p, q)`. This corrects the `σ` returned by `csa_mle_est`; `p` and `q` estimates are unchanged.
- `har_est`: no longer sorts the user's `m` in place.

#### Forecasters.jl
- `csa_forecast`: the one-step forecast reproduced the last observation (MA coefficients were shifted by one lag).
- `har_forecast`: forecasts used the regressors of the previous period (the path was shifted by one step) and failed for `T < 2max(m)`.
- `hitransform`: works on integer input and on series of length ≤ 1.

#### StructuralChanges.jl
- `lm_change_test`: when the maximum came from the backward statistic, the returned location was the length of the backward subsample instead of the break date `T - τ`.
- `lm_change_test`: the series is demeaned before fractional differencing; the truncated filter applied to a nonzero mean left a deterministic term in the regressions.

#### Package structure
- Each submodule is now included once. Previously the source files were re-`include`d inside other submodules, creating up to eight copies of some modules.

#### Documentation
- Fixed docstring examples that did not run or used invalid parameters (e.g. `csa_var_vals` with `q < 1`). A test now runs every docstring example.


<!-- Links generated by Changelog.jl -->

