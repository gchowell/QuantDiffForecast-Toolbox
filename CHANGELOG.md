# Changelog

Versions follow the tags of this repository. Entries list what changed for a user of the toolbox.

## v2.0.0 (2026-11-26)

The first tagged release. The 2024 tutorial paper and video describe the original workflow, which had no version number. This entry summarizes what the repository has gained since the paper; for individual changes, see the commit history.

### Rerun

- **Results obtained by calling a runner without arguments.** The default options changed. `options_fit.m` now fits the 1918 San Francisco influenza series by Poisson likelihood with 300 bootstrap datasets; at the time of the paper it used least squares, 40 bootstrap datasets and a different data file. Analyses that pass their own options function are not affected by this change.
- **Results from earlier revisions in general** can differ, because fitting, bootstrap and scoring code changed between the paper and this release. Record the release you use.

### Added

- **Profile likelihood** (`Run_ProfileLikelihood_ODEModel.m`, `plot_PL_1918_examples.m`).
- **Practical identifiability simulations** (`Run_PracticalIdentifiability_ODEModel.m`, `plotPracticalIdentifiabilityResults.m`, `summarizeIdentifiability.m`), with identifiability indices and bound-hit fractions (`quantdiffPII.m`, `boundHitFractions.m`).
- **Forecast scores for each lead and comparison of models across forecast origins** (`quantdiffForecastScores.m`, `quantdiffCompareForecasts.m`), with flat and trend baselines (`baselineForecast.m`).
- **Forecast ensembles** (`run_ensemble.m`, `getensemble_linearpool_from_curves.m`, `ensembleWeights.m`).
- **Initial conditions for each rolling window** through `vars.initialFcn` (`quantdiffInitialConditions.m`).
- **Estimation method 6,** the sum of absolute deviations with Laplace noise (`estimateLaplaceScale.m`), and a Pearson estimate of the variance-to-mean factor for least squares with negative-binomial noise (`estimatePearsonDispersion.m`).
- **SEIR model with an unreported infectious compartment** and the Spain COVID-19 example (`SEIR_unreported.m`, `R0s_unreported.m`, `input/curve-cases-covid-spain.txt`).
- **Checks** that run without the test framework: `qdfSmokeTest`, `qdfSmokeTest2`, `qdfSmokeTest3`, `qdfTestInitialConditions` and `qdfTestForecastScores`.
- **`make_readme_figures.m`** regenerates four of the README figures.
- **`CITATION.cff` and this changelog.**

### Removed

- `composite1.m`, `SIR1.m`, `R0s_SIR.m` and `getMeanVarianceRatio.m`, which no runner or example used.
- `input/curve-pandemic influenza-SF 1918.txt`. The 1918 series is `input/curve-flu1918SF.txt`.
- `Run_PracticalIndentifiability_ODEModel.m` (misspelled name). Use `Run_PracticalIdentifiability_ODEModel.m`.

### Known limits

- **The standalone Windows and Mac packages** in `stand alone executable` were last updated in April–May 2025 and do not include later changes to the MATLAB source.
- **`plotForecast_ODEModel` rewrites `Forecast-*.csv`** with five columns, dropping `model_median`. Read the file before plotting if you need that column.
- **`Run_Forecasting_ODEModel2.m`** is a variant of the forecast runner used by the practical-identifiability study. A change to one runner does not carry over to the other.

## Before v2.0.0

See the commit history.
