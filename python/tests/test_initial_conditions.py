"""Port of qdfTestInitialConditions.m: objective values on rows 1-17 of the
1918 influenza series against the reference values used by the MATLAB test."""
import os
from types import SimpleNamespace

import numpy as np
import pytest

from quantdiffforecast import (Model, ObjectiveContext, Params, QuantDiffForecastError, Vars,
                               initial_conditions)
from quantdiffforecast.config import BUNDLED_DATA_DIR, check_initial_states
from quantdiffforecast.examples import options_forecast_SEIR_flu1918_IC
from quantdiffforecast.models import SEIR1

data = np.loadtxt(os.path.join(BUNDLED_DATA_DIR, "curve-flu1918SF.txt"))
W = 17
y1 = data[0, 1]
params = Params(label=["beta", "kappa", "gamma", "N"], LB=[0] * 4, UB=[1e9] * 4, initial=[1] * 4,
                fixed=[0] * 4)


def mkvars(fcn=None):
    return Vars(["S", "E", "I", "R", "C"], [550000 - 4, 0, 4, 0, 4], [4], [1], initialFcn=fcn)


def obj(vars_, z, p=params, model=Model(SEIR1, "SEIR")):
    return ObjectiveContext(model, p, vars_, 1, data[:W, 0], data[:W, 1])(z)


def zOf(beta, N):
    return [beta, 1 / 1.9, 1 / 4.1, N, y1, 1, 1]


coupled = lambda p, x0, v: [p[3] - 4, 0, 4, 0, x0]
fPlain = obj(mkvars(), zOf(0.76, 550000))


def test_1_original_behaviour():
    assert np.array_equal(initial_conditions(zOf(0.76, 550000), params, mkvars()), [549996, 0, 4, 0, y1])


def test_2_identity_function():
    assert obj(mkvars(lambda p, x0, v: v.initial), zOf(0.76, 550000)) == fPlain


def test_3_reference_value():
    assert abs(fPlain - (-836.7267)) < 1e-3 * 836.7267


def test_4_redundancy_without_function():
    assert abs(obj(mkvars(), zOf(1.52, 1100000)) - fPlain) < 1e-4 * abs(fPlain)


def test_5_coupling():
    fBig = obj(mkvars(coupled), zOf(0.76, 550000))
    fSmall = obj(mkvars(coupled), zOf(0.76 / 550, 1000))
    assert abs(fBig - fPlain) < 1e-9 * abs(fPlain)
    assert abs(fSmall - 2670.0961) < 1e-3 * 2670.0961


def test_6_fitted_state_pinned():
    v = mkvars(lambda p, x0, v: [p[3] - 4, 0, 4, 0, 999])
    assert initial_conditions(zOf(0.76, 550000), params, v)[4] == y1
    ctx = ObjectiveContext(Model(SEIR1, "SEIR"), params, v, 1, data[:W, 0], data[:W, 1])
    assert ctx.evaluate(zOf(0.76, 550000))[1][0] == y1


def test_7_invalid_results():
    bad = [lambda p, x0, v: [1, 2, 3], lambda p, x0, v: [p[3] - 4, 0, -1, 0, x0],
           lambda p, x0, v: [np.nan, 0, 4, 0, x0], "notAFunction"]
    for b in bad:
        v = mkvars(b)
        with pytest.raises(QuantDiffForecastError, match="InvalidInitialFcn"):
            initial_conditions(zOf(0.76, 550000), params, v)
        assert np.isinf(obj(v, zOf(0.76, 550000)))
    assert np.isinf(obj(mkvars(coupled), zOf(0.76, 2)))


def test_8_window_context_and_warning():
    v = mkvars(lambda p, x0, v: [p[3] - 4 - v.window.priorCum[0], 0, 4, v.window.priorCum[0], x0])
    assert np.array_equal(initial_conditions(zOf(0.76, 550000), params, v), [549996, 0, 4, 0, y1])
    prior = data[:20, 1].sum()
    v.window = SimpleNamespace(tstart=21, rows=np.arange(21, 38), time0=20, firstObs=[320],
                               priorCum=np.array([prior]))
    IC = initial_conditions([0.76, 1 / 1.9, 1 / 4.1, 550000, 320, 1, 1], params, v)
    assert np.array_equal(IC, [550000 - 4 - prior, 0, 4, prior, 320])
    with pytest.warns(UserWarning):
        assert check_initial_states(mkvars(), 3)
    assert not check_initial_states(mkvars(), 1)
    assert not check_initial_states(v, 3)


def test_9_and_10_example_options():
    o = options_forecast_SEIR_flu1918_IC()
    p, v = o.params, o.vars
    IC = initial_conditions(np.concatenate([p.initial, v.initial[v.fit_index], [1, 1]]), p, v)
    assert np.max(np.abs(IC - v.initial)) < 1e-9 and abs(IC[:4].sum() - p.initial[3]) < 1e-9
    f1 = obj(v, [0.58, 1 / 1.9, 1 / 4.1, 550000, 8, 4, y1, 1, 1], p, o.model)
    f2 = obj(v, [0.58, 1 / 1.9, 1 / 4.1, 550000, 8, 50, y1, 1, 1], p, o.model)
    assert abs(f1 - (-859.7854)) < 1e-3 * 859.7854 and f1 == f2
    v.initialFcn = lambda p, x0, v: [p[3] - p[4] - p[5] - v.window.priorCum[0], p[5], p[4],
                                     v.window.priorCum[0], x0]
    fNested = obj(v, [0.76, 1 / 1.9, 1 / 4.1, 550000, 4, 0, y1, 1, 1], p, o.model)
    fBoth = obj(v, [0.58, 1 / 1.9, 1 / 4.1, 550000, 5.7, 10.4, y1, 1, 1], p, o.model)
    assert abs(fNested - fPlain) < 1e-9 * abs(fPlain)
    assert abs(fBoth - (-859.9134)) < 1e-3 * 859.9134
