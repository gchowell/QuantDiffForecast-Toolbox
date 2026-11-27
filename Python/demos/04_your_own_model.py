"""Template: your own data and your own ODE model.

1. Put a header-free text file in the input folder (or pass input_dir=...):
   column 1 = consecutive time index, remaining columns = observed series.
2. Write the model as a Python function dx = f(t, x, p, extra0).
3. Describe the analysis in an Options object and call a runner.
"""
import numpy as np

import quantdiffforecast as q


def logistic_growth(t, x, p, extra0=None):
    """dC/dt = r C (1 - C/K);  p = [r, K];  one state, C (cumulative cases)."""
    r, K = p[0], p[1]
    return np.array([r * x[0] * (1 - x[0] / K)])


def my_options():
    return q.Options(
        cadfilename1="curve-flu1918SF",      # input/curve-flu1918SF.txt
        caddisease="1918 Flu", datatype="cases",
        method1=3, dist1=3,                  # negative binomial likelihood, VAR = mean + alpha*mean
        numstartpoints=10, B=200,
        model=q.Model(fc=logistic_growth, name="Logistic growth"),
        params=q.Params(label=["r", "K"], LB=[0.001, 100], UB=[5, 1e6], initial=[0.3, 20000],
                        fixed=[0, 0],        # 0 = estimate, 1 = hold at `initial`
                        fixI0=1,             # initial value of the fitted state = first observation
                        logscale=True, polish=True),
        vars=q.Vars(label=["C"], initial=[4],
                    fit_index=[0],           # 0-based: the first state is observed ...
                    fit_diff=[1]),           # ... through its increments (daily new cases)
        windowsize1=35, tstart1=1, tend1=1,
        getperformance=True, forecastingperiod=10)


if __name__ == "__main__":
    res = q.run_forecasting(my_options(), seed=1, n_jobs=-1)
    print(res.parameters.T)
    print(res.windows[0].series[0].performance_forecasting)
