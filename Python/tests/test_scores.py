"""Port of qdfTestForecastScores.m: per-horizon scores against the reference
values the MATLAB test uses (computed independently with Python/SciPy)."""
import numpy as np
import pytest

from quantdiffforecast import (baseline_forecast, compute_forecast_performance, compute_wis,
                               forecast_scores)

z = np.linspace(-2.5, 2.5, 201)
mu = np.array([10, 12, 15, 20, 26, 36, 45, 3, 70, 90, 115], float)[:, None]
draws = mu + 1.5 * np.sqrt(mu) * z
obs = np.array([10, 11, 16, 19, 27, 40, 44, 0, 60, 130, 100], float)
datalatest = np.column_stack([np.arange(11.0), obs])
nCalib, H = 5, 6
data1 = datalatest[:nCalib]
refWIS = np.array([3.7014498913, 3.81834247523, 1.49707346475, 6.10643287028, 25.3475931643, 8.49112849215])
refLog = np.array([0.103663461034, 0.0859728985552, 0.595357849149, 0.0895407128021, 0.23058620871, 0.0759060403564])
refAgg3 = np.array([8.1258711591, 17.1010805564])
refAgg6 = 24.8591296959
refBase = np.array([3.6442033913, 3.23532411186, 5.25987900115, 3.80994454666, 41.2780745314, 4.95806486139])
refRel = 0.787354413522
tol = 1e-8
S = forecast_scores(draws, datalatest, nCalib, H)
P = S["perHorizon"]


def test_1_running_means_equal_existing_scores():
    perf = compute_forecast_performance(data1, datalatest, draws, draws, H)
    _, WISF = compute_wis(data1, datalatest, draws, H)
    c = S["cumulative"]
    assert np.allclose(c["MAE"], perf["MAEF"][:, 1], atol=1e-9)
    assert np.allclose(c["MSE"], perf["MSEF"][:, 1], atol=1e-9)
    assert np.allclose(c["coverage95"], perf["PIF"][:, 1], atol=1e-9)
    assert np.allclose(c["WIS"], WISF[:, 1], atol=1e-9)


def test_2_per_horizon_wis():
    assert np.max(np.abs(P["WIS"] - refWIS) / refWIS) < tol
    assert np.max(np.abs(P["WIS_log"] - refLog) / refLog) < tol


def test_3_coverage_and_summary():
    assert list(P["cover50"]) == [1, 1, 1, 1, 0, 1]
    assert list(P["cover95"]) == [1, 1, 1, 1, 0, 1]
    assert list(P["observed"]) == [40, 44, 0, 60, 130, 100]
    assert np.allclose(P["median"], [36, 45, 3, 70, 90, 115], atol=1e-9)
    s = S["summary"]
    assert s["nScored"] == 6
    assert abs(s["WIS"] - 8.16033672633) < 1e-8
    assert abs(s["coverage95"] - 500 / 6) < 1e-9
    assert abs(s["MAE"] - 73 / 6) < 1e-9


def test_4_blocks_of_three():
    T = forecast_scores(draws, datalatest, nCalib, H, aggregate=3)["perHorizon"]
    assert len(T) == 2 and list(T["observed"]) == [84, 290]
    assert np.allclose(T["median"], [84, 275], atol=1e-9)
    assert list(T["step_from"]) == [1, 4] and list(T["step_to"]) == [3, 6]
    assert np.max(np.abs(T["WIS"] - refAgg3) / refAgg3) < tol


def test_5_total_over_horizon():
    T6 = forecast_scores(draws, datalatest, nCalib, H, aggregate=6)["perHorizon"]
    Td = forecast_scores(draws[nCalib:].sum(0)[None, :], np.array([[0, 374.0]]), 0, 1)["perHorizon"]
    T4 = forecast_scores(draws, datalatest, nCalib, H, aggregate=4)["perHorizon"]
    assert T6["observed"][0] == 374 and abs(T6["WIS"][0] - refAgg6) < tol * refAgg6
    assert abs(Td["WIS"][0] - T6["WIS"][0]) < 1e-9
    assert abs(T6["abs_error"][0] - 15) < 1e-9 and abs(P["abs_error"].sum() - 73) < 1e-9
    assert len(T4) == 1 and T4["observed"][0] == 144
    with pytest.raises(ValueError, match="aggregate"):
        forecast_scores(draws, datalatest, nCalib, H, aggregate=7)


def test_6_relative_wis():
    same = forecast_scores(draws, datalatest, nCalib, H, baseline_draws=draws)
    assert np.all(same["perHorizon"]["relWIS"] == 1) and same["summary"]["relWIS"] == 1
    Sb = forecast_scores(draws, datalatest, nCalib, H, baseline_draws=0.8 * draws + 5)
    assert np.max(np.abs(Sb["perHorizon"]["WIS_baseline"] - refBase) / refBase) < tol
    assert abs(Sb["summary"]["relWIS"] - refRel) < tol


def test_7_perfect_forecast():
    yP = np.array([5.0, 0, 120])
    Sp = forecast_scores(np.tile(yP[:, None], (1, 50)), np.column_stack([np.arange(3.0), yP]), 0, 3)
    assert np.all(Sp["perHorizon"]["WIS"] == 0) and np.all(Sp["perHorizon"]["WIS_log"] == 0)
    assert Sp["summary"]["coverage50"] == 100 and Sp["summary"]["coverage95"] == 100
    assert Sp["summary"]["MAE"] == 0


def test_8_horizons_beyond_the_data():
    T = forecast_scores(draws, datalatest[:nCalib + 4], nCalib, H)
    t = T["perHorizon"]
    assert np.all(np.isfinite(t["WIS"][:4])) and np.all(np.isnan(t["WIS"][4:]))
    assert np.all(np.isnan(t["WIS_log"][4:])) and np.all(np.isnan(t["cover95"][4:]))
    assert T["summary"]["nScored"] == 4 and abs(T["summary"]["WIS"] - refWIS[:4].mean()) < 1e-8


def test_9_flat_baseline():
    bd, _ = baseline_forecast(data1, H, "flat", 400, 1)
    Sf = forecast_scores(draws, datalatest, nCalib, H, baseline_draws=bd)
    assert np.all(Sf["perHorizon"]["WIS_baseline"] > 0)
    assert np.max(np.abs(np.median(bd[nCalib:], axis=1) - obs[nCalib - 1])) < 1e-9
    assert np.isfinite(Sf["summary"]["relWIS"])
