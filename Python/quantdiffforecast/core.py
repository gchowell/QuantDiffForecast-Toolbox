"""Numerical core: ODE solves, observation curves, objectives and AICc.

Python port of quantdiffODEOptions.m, quantdiffInitialConditions.m,
quantdiffObservationCurve.m, quantdiffObjectiveValue.m, parameterSearchODE.m,
get_nparams.m, getAICc.m, quantdiffEffectiveN.m and
quantdiffAICcCountAdjustment.m.

The MATLAB objective read six global variables (model, params, vars,
method1, timevect, ydata). Here they live in an explicit
:class:`ObjectiveContext`, which also makes parallel bootstrap fits safe.
"""
from __future__ import annotations

import contextlib
import warnings
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
from scipy.integrate import odeint
from scipy.special import gammaln

from .config import QuantDiffForecastError

# ---------------------------------------------------------------------------
# ODE tolerances (quantdiffODEOptions): RelTol 1e-6, AbsTol 1e-8, NonNegative.
# The MATLAB global QDF_ODE_TOL override becomes a context manager.
# ---------------------------------------------------------------------------
_ODE_TOL = {"rtol": 1e-6, "atol": 1e-8}
ODE_PROFILE_ID = "quantdiff-consistent-v1"


@contextlib.contextmanager
def ode_tolerances(rtol=None, atol=None):
    """Temporarily replace the ODE tolerances used by every solve."""
    old = dict(_ODE_TOL)
    if rtol is not None:
        _ODE_TOL["rtol"] = float(rtol)
    if atol is not None:
        _ODE_TOL["atol"] = float(atol)
    try:
        yield
    finally:
        _ODE_TOL.update(old)


def ode_profile():
    return {"id": ODE_PROFILE_ID, "solver": "LSODA (scipy.integrate.odeint)",
            "RelTol": _ODE_TOL["rtol"], "AbsTol": _ODE_TOL["atol"], "NonNegative": "all states"}


def solve_ode(fc, tvec, IC, p, extra0=None, rtol=None, atol=None):
    """Integrate dx/dt = fc(t, x, p, extra0) and return states at tvec.

    Replaces ``ode15s(model.fc, timevect, IC, odeset(..., 'NonNegative', ...))``.
    LSODA switches automatically between stiff (BDF) and non-stiff methods.
    MATLAB's NonNegative option is emulated: the right-hand side is evaluated
    at max(x, 0) and a component already at zero is not allowed to decrease.
    Raises QuantDiffForecastError if the solve fails or returns non-finite
    states (the MATLAB code raised 'QuantDiffForecast:InvalidTrajectory').
    """
    tvec = np.asarray(tvec, dtype=float)
    IC = np.asarray(IC, dtype=float)
    rtol = _ODE_TOL["rtol"] if rtol is None else rtol
    atol = _ODE_TOL["atol"] if atol is None else atol

    def rhs(t, x):
        xc = np.maximum(x, 0.0)
        dx = np.asarray(fc(t, xc, p, extra0), dtype=float)
        neg = (x <= 0.0) & (dx < 0.0)
        if neg.any():
            dx = np.where(neg, 0.0, dx)
        return dx

    if tvec.size == 1:
        return IC.reshape(1, -1).copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            F, info = odeint(rhs, IC, tvec, tfirst=True, rtol=rtol, atol=atol,
                             full_output=True, mxstep=20000)
        except Exception as exc:  # model raised
            raise QuantDiffForecastError(
                f"QuantDiffForecast:InvalidTrajectory - ODE solve failed: {exc}") from exc
    if info.get("message", "") != "Integration successful." or not np.all(np.isfinite(F)):
        raise QuantDiffForecastError(
            "QuantDiffForecast:InvalidTrajectory - The ODE solve must return complete "
            "finite real states at the requested times.")
    return np.maximum(F, 0.0)


# ---------------------------------------------------------------------------
# initial conditions and observation curve
# ---------------------------------------------------------------------------
def initial_conditions(z, params, vars):
    """quantdiffInitialConditions: IC for the full optimizer vector
    z = [params, X0 of the fitted states, alpha, d]."""
    z = np.asarray(z, dtype=float)
    nfit = len(vars.fit_index)
    npar = params.num
    if z.size < npar + nfit:
        raise QuantDiffForecastError("QuantDiffForecast:InvalidInitialConditions - "
                                     "parameter vector too short.")
    x0 = z[npar:npar + nfit].copy()
    if vars.initialFcn is not None:
        if not callable(vars.initialFcn):
            raise QuantDiffForecastError("QuantDiffForecast:InvalidInitialFcn - "
                                         "vars.initialFcn must be callable.")
        v = vars
        if getattr(vars, "window", None) is None:
            v = SimpleNamespace(**{k: getattr(vars, k) for k in
                                   ("label", "initial", "fit_index", "fit_diff", "initialFcn")})
            v.window = SimpleNamespace(tstart=1, rows=None, time0=0.0, firstObs=x0.copy(),
                                       priorCum=np.zeros(nfit))
        p = z[:npar].copy()
        x0_arg = x0[0] if nfit == 1 else x0
        try:
            IC = np.asarray(vars.initialFcn(p, x0_arg, v), dtype=float).ravel()
        except QuantDiffForecastError:
            raise
        except Exception as exc:
            raise QuantDiffForecastError(
                f"QuantDiffForecast:InvalidInitialFcn - vars.initialFcn raised: {exc}") from exc
        if IC.size != len(vars.initial):
            raise QuantDiffForecastError(
                f"QuantDiffForecast:InvalidInitialFcn - vars.initialFcn must return "
                f"{len(vars.initial)} values, one per state in vars.label.")
        if not np.all(np.isfinite(IC)) or np.any(IC < 0):
            raise QuantDiffForecastError(
                "QuantDiffForecast:InvalidInitialFcn - vars.initialFcn returned non-finite "
                "or negative initial conditions for this parameter vector.")
    else:
        IC = np.asarray(vars.initial, dtype=float).copy()
    IC[np.asarray(vars.fit_index)] = x0
    return IC


def observation_curve(states, vars):
    """quantdiffObservationCurve: stacked observed series (column-major).
    fit_diff == 1 gives abs([C(1); diff(C)])."""
    states = np.asarray(states, dtype=float)
    if states.ndim != 2 or states.size == 0 or not np.all(np.isfinite(states)):
        raise QuantDiffForecastError("QuantDiffForecast:InvalidTrajectory - "
                                     "Observation extraction requires finite states.")
    parts = []
    for idx, dif in zip(vars.fit_index, vars.fit_diff):
        v = states[:, idx]
        if dif == 1:
            v = np.abs(np.concatenate(([v[0]], np.diff(v))))
        parts.append(v)
    return np.concatenate(parts)


# ---------------------------------------------------------------------------
# objective values
# ---------------------------------------------------------------------------
def objective_value(fitcurve, observations, method1, alpha, d):
    """quantdiffObjectiveValue: the six objectives, including the 0.001
    replacement of exact zeros and the omitted likelihood constants."""
    yfit = np.asarray(fitcurve, dtype=float).ravel().copy()
    ydata = np.asarray(observations, dtype=float).ravel()
    if yfit.size == 0 or yfit.size != ydata.size or not (
            np.all(np.isfinite(yfit)) and np.all(np.isfinite(ydata))):
        raise QuantDiffForecastError("QuantDiffForecast:InvalidScoreInput")
    if yfit.sum() == 0:
        raise QuantDiffForecastError("QuantDiffForecast:ZeroCurve")
    yfit[yfit == 0] = 0.001
    yg = np.maximum(ydata, 0)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        if method1 == 0:
            obj = np.sum((ydata - yfit) ** 2)
        elif method1 == 1:
            obj = -np.sum(ydata * np.log(yfit) - yfit)
        elif method1 == 3:
            m = yfit / alpha
            obj = -np.sum(gammaln(yg + m) - gammaln(m) + ydata * np.log(alpha)
                          - (ydata + m) * np.log1p(alpha))
        elif method1 == 4:
            m = 1.0 / alpha
            obj = -np.sum(gammaln(yg + m) - gammaln(m) + ydata * np.log(alpha * yfit)
                          - (ydata + m) * np.log(1 + alpha * yfit))
        elif method1 == 5:
            m = (1.0 / alpha) * yfit ** (2 - d)
            obj = -np.sum(gammaln(yg + m) - gammaln(m) + ydata * np.log(alpha * yfit ** (d - 1))
                          - (ydata + m) * np.log(1 + alpha * yfit ** (d - 1)))
        elif method1 == 6:
            obj = np.sum(np.abs(ydata - yfit))
        else:
            raise QuantDiffForecastError("QuantDiffForecast:UnsupportedMethod")
    obj = float(obj)
    if not np.isfinite(obj):
        raise QuantDiffForecastError("QuantDiffForecast:InvalidObjective")
    return obj


@dataclass
class ObjectiveContext:
    """Everything parameterSearchODE read from globals."""
    model: object
    params: object
    vars: object
    method1: int
    timevect: np.ndarray
    ydata: np.ndarray

    def evaluate(self, z):
        """Strict evaluation: returns (objective, fitcurve, states, IC); raises
        on an invalid trial (MATLAB: nargout > 1)."""
        z = np.asarray(z, dtype=float)
        nfit = len(self.vars.fit_index)
        alpha = z[self.params.num + nfit]
        d = z[self.params.num + nfit + 1]
        IC = initial_conditions(z, self.params, self.vars)
        states = solve_ode(self.model.fc, self.timevect, IC, z, self.params.extra0)
        fitcurve = observation_curve(states, self.vars)
        obj = objective_value(fitcurve, self.ydata, self.method1, alpha, d)
        return obj, fitcurve, states, IC

    def __call__(self, z):
        """Optimizer interface: Inf for any invalid trial."""
        try:
            return self.evaluate(z)[0]
        except Exception:
            return np.inf


# ---------------------------------------------------------------------------
# sample size, number of parameters, AICc
# ---------------------------------------------------------------------------
def effective_n(nObs, nSeries, fixI0):
    """quantdiffEffectiveN: pinned first points do not count."""
    if nObs % nSeries:
        raise QuantDiffForecastError(
            f"{nObs} observations cannot be split into {nSeries} series of equal length.")
    pinned = nSeries if fixI0 == 1 else 0
    return nObs - pinned, {"nObs": nObs, "nSeries": nSeries, "pinnedPoints": pinned}


def get_nparams(method1, dist1, params_num, fixI0_count):
    n = params_num + fixI0_count
    if method1 in (0, 3, 4, 6):
        n += 1
    elif method1 == 5:
        n += 2
    return n


def get_aicc(method1, dist1, params_num, fixI0_count, fval, n, fvalAdjust=0.0):
    """getAICc: returns (AICc, part1, part2, numparams)."""
    k = get_nparams(method1, dist1, params_num, fixI0_count)
    if n - k - 1 <= 0:
        raise QuantDiffForecastError(
            f"getAICc:tooFewData - AICc needs at least k + 2 data points: n = {n}, k = {k}.")
    corr = 2 * k + (2 * k * (k + 1)) / (n - k - 1)
    if method1 == 0:
        part1 = n * np.log(fval / n)
        part2 = n * np.log(2 * np.pi) + n + corr
    elif method1 == 6:
        part1 = 2 * n * np.log(2 * fval / n)
        part2 = 2 * n + corr
    else:
        part1 = 2 * (fval + fvalAdjust)
        part2 = corr
    return part1 + part2, part1, part2, k


def aicc_count_adjustment(method1, ydata, nSeries, fixI0, alpha=1.0, d=1.0):
    """quantdiffAICcCountAdjustment: add sum(log y!) and remove the pinned
    first points so count-likelihood AICc is a full log-likelihood."""
    if method1 not in (1, 3, 4, 5):
        return 0.0
    y = np.asarray(ydata, dtype=float).ravel()
    nT = y.size // nSeries
    pinned = np.zeros(y.size, dtype=bool)
    if fixI0 == 1:
        pinned[::nT] = True
    logfact = float(np.sum(gammaln(np.maximum(y[~pinned], 0) + 1)))
    first = 0.0
    if pinned.any():
        y1 = y[pinned]
        mu1 = y1.copy()
        mu1[mu1 == 0] = 0.001
        first = objective_value(mu1, y1, method1, alpha, d)
    return logfact - first
