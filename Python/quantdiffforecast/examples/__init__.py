"""Ports of the bundled MATLAB options files.

Each MATLAB ``options_*.m`` function becomes a Python function of the same
name that returns an :class:`~quantdiffforecast.config.Options` object.
Copy one of these functions to configure your own analysis.
"""
from __future__ import annotations

import numpy as np

from ..config import Model, Options, Params, Vars
from ..models import EXP, GGM, R0s, R0s_SEIR, R0s_unreported, SEIR1, SEIR_unreported

_SEIR_LABELS = [r"\beta", r"\kappa", r"\gamma", "N"]
_OPTIM = dict(logscale=True, optimOptions=None, useParallel=True, polish=True,
              polish_gradient_tol=1e-2, polish_budget=100)


def _seir_vars(N, I0=4):
    return Vars(label=["S", "E", "I", "R", "C"], initial=[N - I0, 0, I0, 0, I0],
                fit_index=[4], fit_diff=[1])


def _flu(method1, dist1, B, LB, UB, initial, composite=R0s, forecast=True, name="",
         printscreen1=True, numstartpoints=10, optimOptions=None, model_name="SEIR model"):
    params = Params(label=_SEIR_LABELS, LB=LB, UB=UB, initial=initial, fixed=[0, 1, 1, 1],
                    fixI0=1, composite=composite, composite_name="R_0", extra0=None,
                    **{**_OPTIM, "optimOptions": optimOptions})
    return Options(cadfilename1="curve-flu1918SF", caddisease="1918 Flu", datatype="cases",
                   method1=method1, dist1=dist1, numstartpoints=numstartpoints, B=B,
                   model=Model(SEIR1, model_name), params=params, vars=_seir_vars(initial[3]),
                   windowsize1=17, tstart1=1, tend1=1, printscreen1=printscreen1,
                   getperformance=True, forecastingperiod=10 if forecast else 0, name=name)


# --- default options -------------------------------------------------------
def options_fit():
    return _flu(1, 1, 300, [0.01, 0.01, 0.01, 20], [10, 2, 2, 1e6],
                [0.6, 1 / 1.9, 1 / 4.1, 550000], forecast=False, name="options_fit")


def options_forecast():
    return _flu(3, 3, 100, [0.001, 0, 0.1, 550000], [10, 2, 2, 550000],
                [0.76, 1 / 1.9, 1 / 4.1, 550000], composite=R0s_SEIR, name="options_forecast")


# --- 1918 influenza, SEIR, quick start ---------------------------------------
def options_fit_SEIR_flu1918_dist1_3():
    return _flu(3, 3, 300, [0.01, 0.01, 0.01, 20], [10, 2, 2, 1e6],
                [0.6, 1 / 1.9, 1 / 4.1, 550000], forecast=False,
                name="options_fit_SEIR_flu1918_dist1_3")


def options_fit_SEIR_flu1918_dist1_1():
    o = options_fit_SEIR_flu1918_dist1_3()
    o.method1, o.dist1, o.name = 1, 1, "options_fit_SEIR_flu1918_dist1_1"
    return o


def options_fit_SEIR_flu1918_dist1_0():
    o = options_fit_SEIR_flu1918_dist1_3()
    o.method1, o.dist1, o.name = 0, 0, "options_fit_SEIR_flu1918_dist1_0"
    return o


def options_forecast_SEIR_flu1918_dist1_3():
    return _flu(3, 3, 300, [0.001, 0, 0, 20], [10, 2, 2, 1e6],
                [0.76, 1 / 1.9, 1 / 4.1, 550000], name="options_forecast_SEIR_flu1918_dist1_3")


def options_forecast_SEIR_flu1918_dist1_0():
    return _flu(0, 0, 300, [0.001, 0, 0, 20], [10, 2, 2, 1e6],
                [0.76, 1 / 1.9, 1 / 4.1, 550000], composite=R0s_SEIR,
                name="options_forecast_SEIR_flu1918_dist1_0")


def options_forecast_SEIR_flu1918_dist1_1():
    return _flu(1, 1, 200, [0.001, 0, 0.1, 550000], [10, 2, 2, 550000],
                [0.76, 1 / 1.9, 1 / 4.1, 550000], composite=R0s_SEIR,
                optimOptions={"StepTolerance": 1e-8, "FunctionTolerance": 1e-8,
                              "FiniteDifferenceStepSize": 1e-5},
                name="options_forecast_SEIR_flu1918_dist1_1")


def options_forecast_SEIR_flu1918_IC():
    """Window-specific initial conditions (vars.initialFcn) for rolling windows."""
    params = Params(label=[r"\beta", r"\kappa", r"\gamma", "N", "I_0", "E_0"],
                    LB=[0.001, 0, 0, 20, 0.01, 0.01], UB=[10, 2, 2, 1e6, 1e5, 1e5],
                    initial=[0.76, 1 / 1.9, 1 / 4.1, 550000, 4, 4], fixed=[0, 1, 1, 1, 0, 1],
                    fixI0=1, composite=R0s, composite_name="R_0", **_OPTIM)
    p = params.initial
    vars = Vars(label=["S", "E", "I", "R", "C"],
                initial=[p[3] - p[4] - 4 / p[1], 4 / p[1], p[4], 0, 4],
                fit_index=[4], fit_diff=[1],
                # MATLAB: @(p,x0,vars) [p(4)-p(5)-x0/p(2)-vars.window.priorCum, x0/p(2), p(5), vars.window.priorCum, x0]
                initialFcn=lambda p, x0, v: [p[3] - p[4] - x0 / p[1] - v.window.priorCum[0],
                                             x0 / p[1], p[4], v.window.priorCum[0], x0])
    return Options(cadfilename1="curve-flu1918SF", caddisease="1918 Flu", datatype="cases",
                   method1=3, dist1=3, numstartpoints=10, B=300,
                   model=Model(SEIR1, "SEIR model IC"), params=params, vars=vars,
                   windowsize1=17, tstart1=1, tend1=1, printscreen1=True, getperformance=True,
                   forecastingperiod=10, name="options_forecast_SEIR_flu1918_IC")


def options_fit_SEIR_flu1918_PL_C():
    """Profile-likelihood case C: beta and N estimated, S(0) follows N."""
    params = Params(label=_SEIR_LABELS, LB=[0.01, 0.01, 0.01, 400], UB=[10, 2, 2, 1e7],
                    initial=[0.6, 1 / 1.9, 1 / 4.1, 550000], fixed=[0, 1, 1, 0], fixI0=1,
                    composite=R0s, composite_name="R_0")
    vars = _seir_vars(550000)
    vars.initialFcn = lambda p, x0, v: [p[3] - 4, 0, 4, 0, x0]
    return Options(cadfilename1="curve-flu1918SF", caddisease="1918 Flu", datatype="cases",
                   method1=3, dist1=3, numstartpoints=10, B=300, model=Model(SEIR1, "SEIR model"),
                   params=params, vars=vars, windowsize1=17, name="options_fit_SEIR_flu1918_PL_C")


# --- phenomenological models on the 1918 series ------------------------------
def options_forecast_EXP_flu1918_dist1_3():
    params = Params(label=["r"], LB=[0], UB=[10], initial=[0.18], fixed=[0], fixI0=1, **_OPTIM)
    return Options(cadfilename1="curve-flu1918SF", caddisease="1918 Flu", datatype="cases",
                   method1=3, dist1=3, numstartpoints=10, B=300, model=Model(EXP, "EXP model"),
                   params=params, vars=Vars(["C"], [5], [0], [1]), windowsize1=17,
                   forecastingperiod=10, name="options_forecast_EXP_flu1918_dist1_3")


def options_forecast_GGM_flu1918():
    params = Params(label=["r", "p"], LB=[0, 0], UB=[10, 1], initial=[0.18, 0.9], fixed=[0, 0],
                    fixI0=0, **{**_OPTIM, "logscale": False})
    return Options(cadfilename1="curve-flu1918SF", caddisease="1918 Flu", datatype="cases",
                   method1=0, dist1=0, numstartpoints=10, B=300, model=Model(GGM, "GGM"),
                   params=params, vars=Vars(["C"], [5], [0], [1]), windowsize1=17,
                   forecastingperiod=10, name="options_forecast_GGM_flu1918")


# --- COVID-19 Spain, SEIR with under-reporting --------------------------------
def options_fit_SEIR_unreported_covid_spain_dist1_3():
    N = 47332614
    params = Params(label=[r"\beta_0", r"\beta_1", r"\alpha", "q", r"\rho", r"\kappa", r"\gamma", "N"],
                    LB=[0.001, 0, 0.5, 0, 0, 0, 0, N], UB=[4, 2, 1, 5, 1, 2, 2, N],
                    initial=[1.9, 0.15, 0.94, 0.06, 0.82, 1 / 6, 1 / 4, N],
                    fixed=[0, 0, 0, 0, 0, 1, 1, 1], fixI0=1, composite=R0s_unreported,
                    composite_name="R_0", extra0=None, **_OPTIM)
    vars = Vars(["S", "E", "I", "U", "R", "C"], [N - 3, 0, 3, 0, 0, 3], [5], [1])
    return Options(cadfilename1="curve-cases-covid-spain", caddisease="COVID-19", datatype="cases",
                   method1=3, dist1=3, numstartpoints=20, B=300,
                   model=Model(SEIR_unreported, "SEIR covid model with underreporting"),
                   params=params, vars=vars, windowsize1=17,
                   name="options_fit_SEIR_unreported_covid_spain_dist1_3")


def options_fit_SEIR_unreported_covid_spain_dist1_1():
    o = options_fit_SEIR_unreported_covid_spain_dist1_3()
    o.method1, o.dist1, o.name = 1, 1, "options_fit_SEIR_unreported_covid_spain_dist1_1"
    return o


# --- practical identifiability study (simulated exponential growth) -----------
def options_forecast_PII_EXPO_r_dist1_3():
    params = Params(label=["r"], LB=[0], UB=[10], initial=[0.18], fixed=[0], fixI0=1,
                    composite=None, **_OPTIM)
    return Options(cadfilename1="sim-EXP-PII", caddisease="simulated", datatype="cases",
                   method1=3, dist1=3, numstartpoints=10, B=300,
                   model=Model(EXP, "EXP-PII-r-dist1-3"), params=params,
                   vars=Vars(["C"], [5], [0], [1]), windowsize1=30, forecastingperiod=10,
                   printscreen1=False, name="options_forecast_PII_EXPO_r_dist1_3")


def options_forecast_PII_EXPO_r_dist1_5():
    o = options_forecast_PII_EXPO_r_dist1_3()
    o.method1, o.dist1, o.model = 5, 5, Model(EXP, "EXP-PII-r-dist1-5")
    o.name = "options_forecast_PII_EXPO_r_dist1_5"
    o.cadfilename1 = "sim-EXP-PII-NB5"
    return o
