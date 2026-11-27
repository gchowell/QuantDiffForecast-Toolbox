"""Quick start: forecast ten days from the first 17 observations.

MATLAB                                         Python
---------------------------------------------  ---------------------------------------------
forecastOptions = @options_forecast_SEIR_...;  opts = options_forecast_SEIR_flu1918_dist1_3()
rng(1, 'twister');                             (seed=1 below)
Run_Forecasting_ODEModel(forecastOptions, ...  res = run_forecasting(opts, 1, 1, 17, 10, seed=1)
    1, 1, 17, 10);
"""
import quantdiffforecast as q
from quantdiffforecast.examples import options_forecast_SEIR_flu1918_dist1_3

if __name__ == "__main__":       # needed on Windows/macOS when n_jobs > 1
    opts = options_forecast_SEIR_flu1918_dist1_3()
    res = q.run_forecasting(opts, 1, 1, 17, 10, seed=1, n_jobs=-1)   # all CPU cores

    series = res.windows[0].series[0]
    print(series.forecast_table.iloc[17:27])          # rows 18-27: the held-out days
    print(series.scores_by_horizon.iloc[-1][["horizon", "observed", "median", "abs_error",
                                             "cover95", "WIS"]])    # day 10 alone
    print(series.performance_forecasting.iloc[-1])    # days 1-10 together
    print(res.parameters.T)                           # estimates and 95% bootstrap intervals
    print("CSV files and figures are in", res.output_dir)
