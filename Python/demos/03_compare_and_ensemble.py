"""Compare two models across forecast origins, then pool them (takes a few minutes)."""
import quantdiffforecast as q
from quantdiffforecast.examples import (options_forecast_SEIR_flu1918_dist1_3,
                                        options_forecast_SEIR_flu1918_IC)

if __name__ == "__main__":
    models = [options_forecast_SEIR_flu1918_dist1_3(), options_forecast_SEIR_flu1918_IC()]
    labels = ["Quick-start SEIR", "Window-specific SEIR"]
    windows = [[1, 14], [1, 17], [1, 20]]          # [first calibration row, number of rows]
    R = q.compare_forecasts(models, windows, 10, labels=labels, baseline="flat", run=True, seed=1,
                            name="flu1918", n_jobs=-1)
    print(R["byHorizon"][R["byHorizon"].horizon == 10])     # day 10 across origins
    print(R["summary"])                                      # all shared targets

    E = q.run_ensemble(models, labels, -1, 17, 10, "flu1918-equal")   # -1 = equal weights
    print(E["weights"], E["relWIS_ens"])           # relative WIS against flat and trend baselines
