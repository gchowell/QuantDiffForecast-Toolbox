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
   Run_Fit_ODEModel(@options_fit_SEIR_flu1918,1,1,17)
   ```

## Example outputs
   <p align="center">
  <img src="docs/images/model_fit.png" width="48%">
  <img src="docs/images/parameters.png" width="48%">
</p>

3. Visualize the fit and generate forecasts using:

   ```matlab
   plotFit_ODEModel(@options_fit_SEIR_flu1918,1,1,17)
   Run_Forecasting_ODEModel(@options_forecast_SEIR_flu1918,1,1,17,10)
   ```

## Tutorial and Documentation

For a step-by-step guide and a detailed tutorial on how to use this toolbox, please refer to our paper:

**Chowell G., Bleichrodt A., Luo R. (2024)**: "Parameter Estimation and Forecasting with Quantified Uncertainty for ODE Models using QuantDiffForecast: A MATLAB Toolbox and Tutorial". Statistics in Medicine, 43(9), 1826-1848.
[https://onlinelibrary.wiley.com/doi/full/10.1002/sim.10036]

Additionally, a **YouTube tutorial** demonstrating the functionality of the toolbox is available [here](https://www.youtube.com/watch?v=eyyX63H12sY).

## License

This project is licensed under the terms of the Creative Commons Attribution-NonCommercial-NoDerivs License. See [LICENSE](LICENSE) for more details.

## Contact

For questions or support, please contact **Gerardo Chowell** at [gchowell@gsu.edu](mailto:gchowell@gsu.edu).
