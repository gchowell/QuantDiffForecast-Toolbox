# QuantDiffForecast: A MATLAB Toolbox for Parameter Estimation and Forecasting with ODE Models

Welcome to **QuantDiffForecast**, a MATLAB toolbox designed to estimate parameters and generate short-term forecasts with **quantified uncertainty** from dynamical models based on ordinary differential equations (ODEs). This toolbox is **user-friendly**, **flexible**, and suitable for applications across various scientific fields, including **epidemiology**, **population dynamics**, and **systems biology**.

📄 **QuantDiffForecast Tutorial**: [https://onlinelibrary.wiley.com/doi/full/10.1002/sim.10036](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.10036)

🎥 **Video Tutorial**: [https://www.youtube.com/watch?v=eyyX63H12sY&t=41s](https://www.youtube.com/watch?v=eyyX63H12sY&t=41s)

---

## Features

- **Parameter estimation**: Provides methods for parameter estimation using nonlinear least squares (NLS) and maximum likelihood estimation (MLE), with support for Poisson, negative binomial, and normal error structures. The modeler can fit a model to one or multiple time series.
- **Forecasting with quantified uncertainty**: Generates forecasts using parametric bootstrapping to provide uncertainty quantification and prediction intervals.
- **Flexible model input**: Users can define their own ODE models and parameter ranges, supported by customizable input files.
- **Rolling window analysis**: Evaluate parameter stability and forecast performance over time.
- **Illustrative examples**: Includes built-in examples such as epidemic models applied to the 1918 influenza pandemic.

## Getting Started

To get started, you'll need to create a `.txt` file containing your time-series data. Place this file in the `input` folder, and then specify the ODE model and related parameters in the MATLAB `.m` files. 

### Example: SEIR Model for Epidemics

The simplest example provided in this repository is an SEIR (Susceptible-Exposed-Infectious-Removed) model, applied to data from the 1918 influenza pandemic in San Francisco.

1. Specify the SEIR model parameters in `options_fit.m` and `options_forecast.m`.
2. Use the provided script `Run_Fit_ODEModel.m` to estimate parameters and fit the model to data:

   ```matlab
   Run_Fit_ODEModel(@options_fit_SEIR_flu1918_dist1_1,1,1,17)
   ```

## Example outputs
   <p align="center">
  <img src="docs/images/model_fit.png" width="48%">
  <img src="docs/images/parameters.png" width="48%">
</p>

3. Visualize the fit and other related outputs:

   ```matlab
   plotFit_ODEModel(@options_fit_SEIR_flu1918,1,1,17)
   ```

## Additional outputs
   <p align="center">
  <img src="docs/images/stateVars.png" width="48%">
  <img src="docs/images/R0.png" width="48%">
</p>

4. Generate a 10-day ahead forecast:

   ```matlab
   Run_Forecasting_ODEModel(@options_forecast_SEIR_flu1918_dist1_1,1,1,17,10)
   ```

## Example outputs
   <p align="center">
  <img src="docs/images/forecast.png" width="48%">
  <img src="docs/images/forecastingPerformance.png" width="48%">
</p>

5. Visualize the 10-day ahead forecast and other related outputs:

   ```matlab
    plotForecast_ODEModel(@options_forecast_SEIR_flu1918_dist1_1,1,1,17,10)
   ```

## Output Files & Naming Conventions

| File prefix | Produced by | Purpose | Key columns / contents |
|---|---|---|---|
| `AICc-… .csv` | Fit/Forecast | Rolling-window model selection metric(s). | `time`, `AICc` (optionally `AIC`, `BIC` if enabled). |
| `parameters-rollingwindow-… .csv` | Fit/Forecast | Parameter estimates and 95% CIs per window. | `time`, then for each parameter *p*: `p mean`, `p 95% CI LB`, `p 95% CI UB`. |
| `MCSEs-rollingwindow-… .csv` | Fit/Forecast | Monte Carlo standard errors for each parameter per window. | `time`, then `p MCSE` columns. |
| `SCIs-rollingwindow-… .csv` | Fit/Forecast | Identifiability span (SCI) for each parameter. | `time`, then `p SCI` where `SCI = log10(UB/LB)` (smaller ⇒ tighter/clearer ID). |
| `parameters-composite-… .csv` | Fit/Forecast | Composite parameter(s) (e.g., $R_0$) derived from estimates. | `time`, `<name> mean`, `<name> 95% CI LB`, `<name> 95% CI UB`, `<name> SCI`. Uses `params.composite` and `params.composite_name`. |
| `parameters-ODEModel-curve-<file>.mat` | Fit/Forecast | MATLAB struct snapshot for downstream plotting/forecasting. | Model, parameter draws, bootstrap artifacts, etc. |
| `StateVars-… .csv` | Fit/Forecast | Deterministic state trajectories over each window. | `time`, then each state in `vars.label` (e.g., `S,E,I,R,C`). |
| `<param>-histogram-rollingwindow-… .csv` | Fit/Forecast | Bootstrap distribution summaries for a parameter. | Typical columns: `time/window`, `bin_center`, `count` (format may vary). |
| `quantile-… .csv` | **Forecast** | Forecast distribution summaries for the fitted observable. | `time`, `q0.025`, `q0.25`, `q0.50`, `q0.75`, `q0.975` (set may vary). Includes `vars.fit_index-<k>` in name. |
| `Forecast-… .csv` | **Forecast** | Point forecasts/central trajectories. | `time`, `mean`, `median` (and optionally `sd`). Includes `vars.fit_index-<k>`. |
| `performance-calibration-… .csv` | **Forecast** | Aggregated calibration/forecast metrics. | `horizon`, `MAE`, `RMSE`, `MAPE`, `PI_coverage`, `PI_width` (exact set depends on metrics enabled). |


## Tutorial and Documentation

For a step-by-step guide and a detailed tutorial on how to use this toolbox, please refer to our paper:

**Chowell G., Bleichrodt A., Luo R. (2024)**: "Parameter Estimation and Forecasting with Quantified Uncertainty for ODE Models using QuantDiffForecast: A MATLAB Toolbox and Tutorial". Statistics in Medicine, 43(9), 1826-1848.
[https://onlinelibrary.wiley.com/doi/full/10.1002/sim.10036]

Additionally, a **YouTube tutorial** demonstrating the functionality of the toolbox is available [here](https://www.youtube.com/watch?v=eyyX63H12sY).

## License

This project is licensed under the terms of the Creative Commons Attribution-NonCommercial-NoDerivs License. See [LICENSE](LICENSE) for more details.

## Contact

For questions or support, please contact **Gerardo Chowell** at [gchowell@gsu.edu](mailto:gchowell@gsu.edu).
