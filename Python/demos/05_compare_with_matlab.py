"""Check the Python port against a MATLAB run of the quick start.

Step 1, in MATLAB (from the toolbox code folder). Do NOT call
plotForecast_ODEModel before copying the files: it rewrites the Forecast CSV.

    rng(1, 'twister');
    Run_Forecasting_ODEModel(@options_forecast_SEIR_flu1918_dist1_3, 1, 1, 17, 10);

Step 2, here:

    python demos/05_compare_with_matlab.py "/path/to/forecasting_odemodels code/output"

The script makes three kinds of comparison, from strictest to loosest.

  A  No randomness, no optimizer. Python evaluates the likelihood and AICc AT
     MATLAB'S parameter estimates, and scores MATLAB'S OWN predictive draws.
     These must agree almost exactly; a difference here is a porting error.
  B  Point estimates. Deterministic in each language, but found by different
     optimizers and ODE solvers, so they agree to a few significant digits.
  C  Bootstrap intervals and forecast scores. MATLAB and NumPy draw different
     random numbers, so these agree only to Monte Carlo error. The script
     measures that error by running Python with several seeds and asks whether
     the MATLAB value is a plausible member of that spread.
"""
import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd

import quantdiffforecast as q
from quantdiffforecast.examples import options_forecast_SEIR_flu1918_dist1_3

W, H = 17, 10
SERIES = ("model_name-SEIR model-fit_index-5-fixI0-1-method-3-dist-3-tstart-1-tend-1-"
          "calibrationperiod-17-horizon-10-1918 Flu-cases.csv")
WINDOW = ("model_name-SEIR model-fixI0-1-method-3-dist-3-tstart-1-tend-1-"
          "calibrationperiod-17-horizon-10-1918 Flu-cases.csv")
MATFILE = ("Forecast-ODEModel-curve-flu1918SF.txt-model_name-SEIR model-fit_index-5-fixI0-1-"
           "method-3-dist-3-tstart-1-tend-1-calibrationperiod-17-forecastingperiod-10.mat")
results = []


def report(tier, name, matlab, python, tol, note=""):
    diff = abs(matlab - python)
    ok = diff <= tol
    results.append(ok)
    print(f"  {'PASS' if ok else 'FAIL'}  {tier}  {name:<34s} MATLAB {matlab:>12.6g}   Python {python:>12.6g}"
          f"   |diff| {diff:.2e}  (tolerance {tol:.2e}) {note}")


def read(folder, name):
    path = os.path.join(folder, name)
    if not os.path.isfile(path):
        sys.exit(f"Missing MATLAB file:\n  {path}\nRun the quick start in MATLAB first (see the "
                 "top of this script).")
    return pd.read_csv(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("matlab_output", help="the MATLAB toolbox's output folder")
    ap.add_argument("--seeds", type=int, default=5, help="Python runs used to measure Monte Carlo error")
    ap.add_argument("--jobs", type=int, default=-1, help="CPU cores for the bootstrap (-1 = all)")
    a = ap.parse_args()
    warnings.simplefilter("ignore")
    mdir = a.matlab_output

    par = read(mdir, "parameters-rollingwindow-" + WINDOW).to_numpy(float)[0]
    aic = read(mdir, "AICc-" + WINDOW).to_numpy(float)[0]
    comp = read(mdir, "parameters-composite-" + WINDOW).to_numpy(float)[0]
    perf = read(mdir, "performance-forecasting-" + SERIES).to_numpy(float)
    fc = read(mdir, "Forecast-" + SERIES)
    byh = read(mdir, "scores-byhorizon-" + SERIES)
    # parameters CSV: time, then (estimate, 2.5%, 97.5%) for beta, kappa, gamma, N, X0, alpha
    m_beta, m_alpha = par[1:4], par[-3:]
    m_P = np.array([par[1], par[4], par[7], par[10], par[13], par[-3], 1.0])

    opts = options_forecast_SEIR_flu1918_dist1_3().copy(output_dir=os.path.join(mdir, "_python_check"))
    data = np.loadtxt(opts.input_path())

    # ------------------------------------------------------------------ A
    print("\nA. Exact checks (no random numbers, no optimizer)")
    ctx = q.ObjectiveContext(opts.model, opts.params, opts.vars, 3, data[:W, 0], data[:W, 1])
    fval = ctx.evaluate(m_P)[0]
    adj = q.aicc_count_adjustment(3, data[:W, 1], 1, 1, m_P[5], m_P[6])
    aicc = q.get_aicc(3, 3, 1, 0, fval, W - 1, adj)[0]
    report("A", "AICc at MATLAB's estimates", aic[1], aicc, 1e-3,
           "<- ODE solver + likelihood + AICc")
    matpath = os.path.join(mdir, MATFILE)
    try:
        from scipy.io import loadmat
        M = loadmat(matpath, variable_names=["forecast2", "datalatest", "data1"])
        S = q.forecast_scores(M["forecast2"], M["datalatest"][:, :2], M["data1"].shape[0], H)["perHorizon"]
        report("A", "WIS, day 10, on MATLAB's draws", byh["WIS"].iloc[-1], S["WIS"].iloc[-1], 1e-6,
               "<- quantile rule + WIS")
        report("A", "mean WIS, days 1-10, MATLAB draws", perf[-1, 4], S["WIS"].mean(), 1e-6)
        report("A", "median, day 10, MATLAB's draws", byh["median"].iloc[-1], S["median"].iloc[-1], 1e-6)
        cov = 100 * S["cover95"].mean()
        report("A", "95% coverage on MATLAB's draws", perf[-1, 3], cov, 1e-6)
    except Exception as exc:
        print(f"  SKIP  A  scoring MATLAB's own draws: could not read {MATFILE}\n"
              f"           ({type(exc).__name__}: {exc}). If the file was saved as -v7.3, re-save "
              "forecast2, datalatest and data1 with save(..., '-v7').")

    # ------------------------------------------------------------------ B, C
    print(f"\nRunning the Python quick start with {a.seeds} seeds to measure Monte Carlo error ...")
    rows = []
    for seed in range(1, a.seeds + 1):
        r = q.run_forecasting(opts, 1, 1, W, H, seed=seed, n_jobs=a.jobs, verbose=False, plot=False,
                              save=False)
        w, s = r.windows[0], r.windows[0].series[0]
        f = s.forecast_table
        rows.append({"beta": w.P[0], "alpha": w.P[5], "R0": w.composite[0], "AICc": w.AICc[0],
                     "beta 2.5%": w.param_estims[0, 1], "beta 97.5%": w.param_estims[0, 2],
                     "alpha 97.5%": w.param_estims[5, 2],
                     "R0 2.5%": w.composite[1], "R0 97.5%": w.composite[2],
                     "median, day 10": f["median"].iloc[-1], "lower limit, day 10": f["LB"].iloc[-1],
                     "upper limit, day 10": f["UB"].iloc[-1],
                     "median, day 5": f["median"].iloc[W + 4],
                     "MAE, days 1-10": s.performance_forecasting["MAE"].iloc[-1],
                     "mean WIS, days 1-10": s.performance_forecasting["WIS"].iloc[-1],
                     "95% coverage (%)": s.performance_forecasting["Coverage 95%PI"].iloc[-1]})
        print(f"  seed {seed} done")
    P = pd.DataFrame(rows)
    matlab = {"beta": m_beta[0], "alpha": m_alpha[0], "R0": comp[1], "AICc": aic[1],
              "beta 2.5%": m_beta[1], "beta 97.5%": m_beta[2], "alpha 97.5%": m_alpha[2],
              "R0 2.5%": comp[2], "R0 97.5%": comp[3],
              "median, day 10": fc["median"].iloc[W + H - 1], "lower limit, day 10": fc["LB"].iloc[W + H - 1],
              "upper limit, day 10": fc["UB"].iloc[W + H - 1], "median, day 5": fc["median"].iloc[W + 4],
              "MAE, days 1-10": perf[-1, 1], "mean WIS, days 1-10": perf[-1, 4],
              "95% coverage (%)": perf[-1, 3]}

    print("\nB. Point estimates (deterministic; different optimizer and ODE solver)")
    report("B", "beta", matlab["beta"], P["beta"].iloc[0], 2e-3 * abs(matlab["beta"]))
    report("B", "R0 = beta/gamma", matlab["R0"], P["R0"].iloc[0], 2e-3 * abs(matlab["R0"]))
    report("B", "alpha (dispersion; flat direction)", matlab["alpha"], P["alpha"].iloc[0],
           2e-2 * abs(matlab["alpha"]))
    report("B", "AICc", matlab["AICc"], P["AICc"].iloc[0], 1e-2)

    print("\nC. Bootstrap intervals and scores (agree to Monte Carlo error)")
    print("   tolerance = 3 x Python's seed-to-seed SD x sqrt(1 + 1/seeds), with a small floor")
    floors = {"95% coverage (%)": 10.0}             # one of ten days in or out
    for k in list(matlab)[4:]:
        sd = P[k].std(ddof=1) if a.seeds > 1 else 0.0
        tol = max(3 * sd * np.sqrt(1 + 1 / a.seeds), floors.get(k, 0.02 * abs(P[k].mean())))
        report("C", k, matlab[k], P[k].mean(), tol, f"[Python SD {sd:.3g}]")

    n = len(results)
    print(f"\n{sum(results)} of {n} comparisons passed.")
    if not all(results):
        print("A failure in A is a porting error. A failure in B points to the optimizer (compare\n"
              "the objective values). A single marginal failure in C is expected now and then:\n"
              "rerun with --seeds 10, or rerun MATLAB with another seed, before concluding anything.")
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
