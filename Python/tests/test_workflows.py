"""End-to-end checks of the main workflows (the role of qdfSmokeTest*.m).

Small bootstrap sizes keep the run to about a minute. Point estimates are
deterministic and are compared with the values the full-size runs give; the
README of the MATLAB toolbox documents the same example (17 calibration rows
of the 1918 influenza series, ten-day forecast).
"""
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pytest

import quantdiffforecast as q
from quantdiffforecast.examples import (options_fit_SEIR_flu1918_dist1_3,
                                        options_forecast_EXP_flu1918_dist1_3,
                                        options_forecast_SEIR_flu1918_dist1_3,
                                        options_forecast_SEIR_flu1918_IC)

warnings.simplefilter("ignore")


@pytest.fixture(scope="module")
def out(tmp_path_factory):
    return str(tmp_path_factory.mktemp("output"))


@pytest.fixture(scope="module")
def quick(out):
    o = options_forecast_SEIR_flu1918_dist1_3().copy(output_dir=out, B=30)
    return q.run_forecasting(o, 1, 1, 17, 10, seed=1, verbose=False)


def test_quickstart_point_fit(quick):
    w = quick.windows[0]
    assert abs(w.P[0] - 0.7718) < 2e-3                 # beta
    assert abs(w.P[5] - 1.561) < 2e-2                  # NB dispersion alpha
    assert abs(w.fval - (-846.4781)) < 1e-3
    assert abs(w.composite[0] - 3.1645) < 1e-2         # R0 = beta / gamma
    assert w.fit_diagnostics["converged"]
    assert quick.accepted_bootstrap == quick.requested_bootstrap == 30


def test_quickstart_outputs(quick, out):
    s = quick.windows[0].series[0]
    assert s.forecast2.shape == (27, 30 * 20)          # 20 predictive draws per refit
    assert list(s.forecast_table.columns) == ["time", "data", "median", "LB", "UB", "model_median"]
    assert s.forecast_table["data"].iloc[17] == 124 and s.forecast_table["data"].iloc[26] == 1104
    pf = s.performance_forecasting
    assert len(pf) == 10 and 30 < pf["WIS"].iloc[-1] < 80 and pf["Coverage 95%PI"].iloc[-1] >= 50
    assert np.isclose(s.scores_by_horizon["WIS"].mean(), pf["WIS"].iloc[-1])
    sfx = ("model_name-SEIR model-fit_index-5-fixI0-1-method-3-dist-3-tstart-1-tend-1-"
           "calibrationperiod-17-horizon-10-1918 Flu-cases.csv")
    for prefix in ("Forecast-", "performance-forecasting-", "scores-byhorizon-", "quantile-"):
        assert os.path.isfile(os.path.join(out, prefix + sfx))
    assert any(f.startswith("fig-forecast") for f in os.listdir(out))


def test_reproducible_with_seed(out):
    o = options_forecast_EXP_flu1918_dist1_3().copy(output_dir=out, B=8)
    a = q.run_forecasting(o, 1, 1, 17, 5, seed=7, verbose=False, plot=False, save=False)
    b = q.run_forecasting(o, 1, 1, 17, 5, seed=7, n_jobs=2, verbose=False, plot=False, save=False)
    assert np.array_equal(a.windows[0].Phatss, b.windows[0].Phatss)
    assert np.array_equal(a.forecast_model12, b.forecast_model12)


def test_fit_only(out):
    o = options_fit_SEIR_flu1918_dist1_3().copy(output_dir=out, B=10)
    r = q.run_fit(o, 1, 1, 17, seed=1, verbose=False, plot=False)
    assert abs(r.windows[0].P[0] - 0.7718) < 2e-3 and r.forecastingperiod == 0
    assert r.windows[0].series[0].forecast2.shape[0] == 17


@pytest.mark.parametrize("method1,dist1", [(0, 0), (0, 1), (0, 2), (1, 1), (4, 4), (5, 5), (6, 6)])
def test_every_method(out, method1, dist1):
    o = options_forecast_SEIR_flu1918_dist1_3().copy(output_dir=out, B=6, method1=method1, dist1=dist1)
    r = q.run_forecasting(o, 1, 1, 17, 10, seed=1, verbose=False, plot=False, save=False)
    assert 0.70 < r.windows[0].P[0] < 0.85 and np.all(np.isfinite(r.forecast_model12))


def test_unsupported_method_and_negative_counts(out, tmp_path):
    o = options_forecast_SEIR_flu1918_dist1_3().copy(output_dir=out, method1=2)
    with pytest.raises(q.QuantDiffForecastError, match="UnsupportedMethod"):
        q.run_forecasting(o, verbose=False, plot=False)
    d = np.loadtxt(o.input_path())
    d[5, 1] = -3
    np.savetxt(tmp_path / "neg.txt", d)
    o = options_forecast_SEIR_flu1918_dist1_3().copy(output_dir=out, input_dir=str(tmp_path),
                                                     cadfilename1="neg")
    with pytest.raises(q.QuantDiffForecastError, match="NegativeCounts"):
        q.run_forecasting(o, verbose=False, plot=False)


def test_rolling_windows_with_initial_fcn(out):
    o = options_forecast_SEIR_flu1918_IC().copy(output_dir=out, B=6)
    r = q.run_forecasting(o, 1, 2, 17, 10, seed=1, verbose=False, plot=False)
    assert len(r.windows) == 2 and list(r.parameters["time"]) == [1, 2]
    assert abs(r.windows[0].P[0] - 0.5765) < 5e-3      # beta of the window-specific model


def test_compare_and_ensemble(quick, out):
    a = options_forecast_SEIR_flu1918_dist1_3().copy(output_dir=out, B=30)
    b = options_forecast_EXP_flu1918_dist1_3().copy(output_dir=out, B=30)
    R = q.compare_forecasts([a, b], [[1, 17]], 10, labels=["SEIR", "EXP"], run=True, seed=1,
                            name="test", verbose=False)
    assert list(R["summary"]["model"]) == ["SEIR", "EXP"] and R["summary"]["nTargets"].iloc[0] == 10
    assert np.all(np.isfinite(R["summary"]["relWIS"]))
    E = q.run_ensemble([a, b], ["SEIR", "EXP"], -1, 17, 10, "test-equal", plot=False, verbose=False)
    assert np.allclose(E["weights"], 0.5) and np.all(np.isfinite(E["relWIS_ens"]))
    with pytest.raises(FileNotFoundError, match="missing"):
        q.compare_forecasts([a], [[46, 7]], 3, save=False, verbose=False)


def test_profile_likelihood(out):
    o = options_fit_SEIR_flu1918_dist1_3().copy(output_dir=out)
    PL = q.run_profile_likelihood(o, 1, 17, seed=1, paramIndex=[0], numstartpoints=2, verbose=False)
    p = PL["param"][0]
    assert PL["referenceStatus"] == "validated" and p["status"].startswith("finite interval")
    assert abs(p["CI"][0] - 0.7416) < 5e-3 and abs(p["CI"][1] - 0.8017) < 5e-3


def test_simulate_and_recover(out, tmp_path):
    from quantdiffforecast.examples import options_forecast_PII_EXPO_r_dist1_3
    o = options_forecast_PII_EXPO_r_dist1_3().copy(output_dir=out, input_dir=str(tmp_path), B=12)
    F, curves, truth = q.run_simulate(o, 30, 5, seed=1)
    assert curves.shape == (60, 2) and truth == {"dist1": 3, "factor1": 5.0, "d": 1.0}
    q.run_practical_identifiability(o, 30, 5, replicates=2, output_dir=str(tmp_path / "study"),
                                    verbose=False)
    S = q.summarize_identifiability(str(tmp_path / "study"), [30], 2, [0, 2], [0.18, 5], ["r", "alpha"])
    assert S["converged"].iloc[0] == 2 and abs(S["bias_MLE"].iloc[0]) < 0.01


def test_two_observed_series(out, tmp_path):
    """Prevalence (a state level) and incidence (increments) fitted together."""
    truth = np.array([0.8, 1 / 1.9, 1 / 4.1, 100000.0])
    F = q.solve_ode(q.models.SEIR1, np.arange(40.0), [100000 - 10, 0, 10, 0, 10], truth)
    rng = np.random.default_rng(3)
    prev = rng.poisson(F[:, 2]).astype(float)
    inc = rng.poisson(np.abs(np.r_[F[0, 4], np.diff(F[:, 4])])).astype(float)
    prev[0] = inc[0] = 10
    np.savetxt(tmp_path / "two-series.txt", np.column_stack([np.arange(40.0), prev, inc]))
    o = q.Options(
        cadfilename1="two-series", caddisease="sim", datatype="cases", method1=1, dist1=1,
        numstartpoints=5, B=8, model=q.Model(q.models.SEIR1, "SEIR two series"),
        params=q.Params(label=["beta", "kappa", "gamma", "N"], LB=[0.01, 0.01, 0.01, 20],
                        UB=[10, 2, 2, 1e6], initial=[0.6, 1 / 1.9, 1 / 4.1, 100000],
                        fixed=[0, 1, 0, 1], fixI0=1, logscale=True, polish=True),
        vars=q.Vars(["S", "E", "I", "R", "C"], [100000 - 10, 0, 10, 0, 10], [2, 4], [0, 1]),
        windowsize1=25, forecastingperiod=7, input_dir=str(tmp_path), output_dir=out)
    r = q.run_forecasting(o, 1, 1, 25, 7, seed=2, verbose=False, plot=False)
    w = r.windows[0]
    assert abs(w.P[0] - 0.8) < 0.05 and abs(w.P[2] - 1 / 4.1) < 0.02      # beta and gamma recovered
    assert [s.fit_index for s in w.series] == [2, 4]
    assert w.series[0].forecast2.shape[0] == 32 and w.series[1].data1[0, 1] == 10
