"""Parameter estimation: Python port of fit_model.m.

MATLAB -> Python mapping used here

  fmincon('sqp')               scipy.optimize.minimize(method='L-BFGS-B')
                               The problem has box bounds only, and SciPy's own
                               SQP (SLSQP) has no step-size stopping rule like
                               fmincon's StepTolerance: on a finite-difference
                               objective it can iterate for tens of thousands of
                               evaluations. L-BFGS-B reaches the same optima and
                               stops reliably. Select SLSQP with
                               params.optimOptions = {"Algorithm": "SLSQP"}.
  MultiStart + CustomStartPointSet
                               a loop of local solves over the same start points
  lhsdesign(...,'maximin')     scipy.stats.qmc.LatinHypercube, best of 50 designs
                               by minimum pairwise distance
  fminsearch (polish)          scipy.optimize.minimize(method='Nelder-Mead')

Everything else is kept: the optimizer works only on the free coordinates of
z = [params, X0, alpha, d], each mapped to [0, 1] (optionally after log10),
the objective is divided by |f(seed)|, failed ODE solves are mapped above every
valid value by the same bounded transform (localBound), and the rescale pass,
the stationarity check, the Nelder-Mead polish and the confirmation step are
reproduced, together with the fitDiagnostics fields.
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.optimize import minimize
from scipy.stats import qmc

from .config import QuantDiffForecastError, check_method
from .core import ObjectiveContext, observation_curve, ode_profile, solve_ode

_T = 1e6
_FAILGAP = 1e3


def _bound(v):
    if not np.isfinite(v):
        return _T + _FAILGAP
    if v > _T:
        return _T + math.log1p(v - _T)
    return v


def _unbound(g):
    if not np.isfinite(g) or g >= _T + _FAILGAP:
        return np.inf
    if g > _T:
        return _T + math.expm1(g - _T)
    return g


@dataclass
class FitResult:
    P: np.ndarray                 # full vector [params, X0, alpha, d]
    residual: np.ndarray
    fitcurve: np.ndarray
    forecastcurve: np.ndarray
    timevect2: np.ndarray
    initialguess: np.ndarray
    fval: float
    F1: np.ndarray                # states on the calibration grid
    F2: np.ndarray                # states on the forecast grid
    diagnostics: dict = field(default_factory=dict)


def _lhs_maximin(n, d, rng, iterations=50):
    best, bestscore = None, -np.inf
    for _ in range(iterations):
        X = qmc.LatinHypercube(d=d, seed=rng).random(n)
        if n < 2:
            return X
        diff = X[:, None, :] - X[None, :, :]
        dist = np.sqrt((diff ** 2).sum(-1))
        score = dist[np.triu_indices(n, 1)].min()
        if score > bestscore:
            best, bestscore = X, score
    return best


def _make_starts(nLHS, nJit, d, LB, UB, span, z0, logmask, rng):
    rows = []
    if nLHS > 0:
        X = _lhs_maximin(nLHS, d, rng)
        rows.append(LB + X * span)
    if nJit > 0:
        S2 = np.tile(z0, (nJit, 1)) * 10.0 ** (0.1 * rng.standard_normal((nJit, d)))
        zero_cols = (z0 == 0) & ~logmask
        if zero_cols.any():
            S2[:, zero_cols] = z0[zero_cols] + 0.01 * rng.standard_normal(
                (nJit, zero_cols.sum())) * span[zero_cols]
        if logmask.any():
            S2[:, logmask] = z0[logmask] + 0.1 * rng.standard_normal((nJit, logmask.sum()))
        S2 = np.clip(S2, LB, UB)
        rows.append(S2)
    rows.append(z0[None, :])
    starts = np.vstack(rows)
    # unique rows, stable
    _, idx = np.unique(starts, axis=0, return_index=True)
    starts = starts[np.sort(idx)]
    inb = np.all((starts >= LB - 1e-12) & (starts <= UB + 1e-12), axis=1) & np.all(np.isfinite(starts), axis=1)
    starts = starts[inb]
    if starts.shape[0] == 0:
        starts = z0[None, :].copy()
    return starts


def _scaled_gradient(fun, q, h):
    n = q.size
    try:
        fp = fun(q)
        g = np.zeros(n)
        for j in range(n):
            e = np.zeros(n)
            e[j] = h
            if q[j] + h <= 1:
                g[j] = (fun(q + e) - fp) / h
            else:
                g[j] = (fp - fun(q - e)) / h
        gp = g.copy()
        gp[(q <= 1e-8) & (g > 0)] = 0
        gp[(q >= 1 - 1e-8) & (g < 0)] = 0
        return gp, float(np.max(np.abs(gp))) if n else 0.0
    except Exception:
        return np.full(n, np.nan), np.nan


def fit_model(data1, params0, numstartpoints, DT, model, params, vars, forecastingperiod,
              method1, pin=None, rng=None) -> FitResult:
    """Fit the model to one calibration window.

    data1    window rows [time, series_1, ..., series_s]
    params0  seed for the full vector [params, X0, alpha, d]
    method1  estimation method (was a global in MATLAB)
    pin      optional vector of the same length; non-NaN entries are held fixed
    rng      numpy Generator for the start points (replaces the global RNG)
    """
    rng = np.random.default_rng() if rng is None else rng
    data1 = np.asarray(data1, dtype=float)
    timevect = data1[:, 0] * DT
    I0 = data1[0, 1:].copy()
    z = np.asarray(params0, dtype=float).copy()

    pLB = params.LB.copy()
    pUB = params.UB.copy()
    for i in range(params.num):
        if params.fixed[i]:
            pLB[i] = z[i]
            pUB[i] = z[i]

    check_method(method1, "fit_model")
    if method1 in (0, 1, 6):
        LBe, UBe = [0.0, 0.0], [0.0, 0.0]
    elif method1 in (3, 4):
        LBe, UBe = [1e-8, 1.0], [1e4, 1.0]
    else:  # 5
        LBe, UBe = [1e-8, 0.6], [1e4, 5.0]

    if params.fixI0 == 1:
        LB = np.concatenate([pLB, I0, LBe])
        UB = np.concatenate([pUB, I0, UBe])
    else:
        LB = np.concatenate([pLB, np.zeros(I0.size) + 0.001, LBe])
        UB = np.concatenate([pUB, np.abs(data1[:, 1:]).sum(axis=0), UBe])

    if pin is not None:
        pin = np.asarray(pin, dtype=float).ravel()
        if pin.size != LB.size:
            raise QuantDiffForecastError(
                f"fit_model:pinLength - pin must have {LB.size} entries (params, X0, alpha, d).")
        pinned = ~np.isnan(pin)
        LB[pinned] = pin[pinned]
        UB[pinned] = pin[pinned]
        z[pinned] = pin[pinned]

    ydata = data1[:, 1:].flatten(order="F")
    if not np.all(np.isfinite(data1)):
        raise QuantDiffForecastError("fit_model:invalidData - non-finite data (NaN/Inf) in data1.")

    ctx = ObjectiveContext(model, params, vars, method1, timevect, ydata)
    f = ctx

    d = z.size
    if LB.size != d:
        raise QuantDiffForecastError(f"Length mismatch: numel(z)={d}, numel(LB)={LB.size}")
    flip = LB > UB
    if flip.any():
        LB[flip], UB[flip] = UB[flip].copy(), LB[flip].copy()
    if not (np.all(np.isfinite(LB)) and np.all(np.isfinite(UB))):
        raise QuantDiffForecastError("LB/UB must be finite.")

    z0 = np.minimum(np.maximum(z, LB), UB)
    nanz = np.isnan(z0)
    z0[nanz] = ((LB + UB) / 2)[nanz]

    useLog = bool(params.logscale)
    scaleObj = bool(params.scaleObjective)
    doPolish = bool(params.polish)
    polishGradTol = float(params.polish_gradient_tol)
    polishBudget = float(params.polish_budget)

    free = np.where(UB > LB)[0]
    nFree = free.size
    span = UB - LB
    tf = ((LB > 0) & (UB > LB)) if useLog else np.zeros(d, dtype=bool)

    def toT(Z):
        T = np.array(Z, dtype=float, copy=True)
        if tf.any():
            T[..., tf] = np.log10(np.maximum(T[..., tf], np.finfo(float).tiny))
        return T

    def fromT(T):
        Z = np.array(T, dtype=float, copy=True)
        if tf.any():
            Z[..., tf] = 10.0 ** Z[..., tf]
        return Z

    LBt, UBt = toT(LB), toT(UB)
    spanT = UBt - LBt
    z0T = toT(z0)

    def toQ(Z):
        T = toT(np.atleast_2d(Z))
        U = (T[:, free] - LBt[free]) / spanT[free]
        return np.clip(U, 0, 1)

    def fromQ(U):
        U = np.clip(np.atleast_1d(np.asarray(U, dtype=float)), 0, 1)
        T = z0T.copy()
        T[free] = LBt[free] + U * spanT[free]
        return fromT(T)

    q0 = toQ(z0)[0] if nFree else np.zeros(0)
    f0 = f(z0)
    fscale = abs(f0) if (scaleObj and np.isfinite(f0) and abs(f0) >= 1) else 1.0

    def fq(U):
        return f(fromQ(U))

    # options (fmincon defaults of the toolbox, overridable)
    oo = {"StepTolerance": 1e-10, "FunctionTolerance": 1e-10, "OptimalityTolerance": 1e-8,
          "FiniteDifferenceStepSize": 1e-6, "MaxFunctionEvaluations": 20000,
          "MaxIterations": 20000}
    if params.optimOptions:
        oo.update(params.optimOptions)
    fdStep = float(np.atleast_1d(oo["FiniteDifferenceStepSize"])[0])

    algorithm = str(oo.get("Algorithm", "L-BFGS-B"))
    counter = {"n": 0}

    def local_solve(fun, x0):
        """One bounded local minimization in unit coordinates.

        Returns (x, f, exitflag) with fmincon-like flags: 1 = converged,
        2 = stopped because no further progress was possible (the analogue
        of fmincon's step-tolerance exit), 0 = iteration/evaluation limit.
        """
        best = {"x": np.asarray(x0, float).copy(), "f": np.inf}

        def tracked(u):
            counter["n"] += 1
            v = fun(u)
            if v < best["f"]:
                best["f"], best["x"] = v, np.array(u, float)
            return v

        maxit = int(min(oo["MaxIterations"], 100 * max(nFree, 1) + 200))
        maxfun = int(min(oo["MaxFunctionEvaluations"], 400 * (nFree + 1)))
        if algorithm.upper() == "SLSQP":
            res = minimize(tracked, x0, method="SLSQP", bounds=[(0.0, 1.0)] * nFree,
                           options={"ftol": oo["FunctionTolerance"], "eps": fdStep,
                                    "maxiter": maxit})
            status = int(getattr(res, "status", 1))
            flag = 1 if status == 0 else (2 if status == 8 else 0)
        else:
            res = minimize(tracked, x0, method="L-BFGS-B", bounds=[(0.0, 1.0)] * nFree,
                           options={"ftol": oo["FunctionTolerance"],
                                    "gtol": oo["OptimalityTolerance"], "eps": fdStep,
                                    "maxiter": maxit, "maxfun": maxfun, "maxls": 30})
            status = int(getattr(res, "status", 1))
            flag = 1 if status == 0 else (2 if status == 2 else 0)
        # forward differences evaluate points the line search never accepts;
        # return the best point actually seen (never worse than res.x)
        x, fv = (best["x"], best["f"]) if best["f"] <= res.fun else (res.x, float(res.fun))
        if not np.isfinite(fv):
            flag = 0
        return np.clip(x, 0, 1), float(fv), flag

    # start points
    nStarts = max(0, int(math.floor(numstartpoints)))
    nLHS = nStarts
    nJit = min(3, max(1, nStarts))
    if useLog:
        def makeStarts():
            return fromT(_make_starts(nLHS, nJit, d, LBt, UBt, spanT, z0T, tf, rng))
    else:
        def makeStarts():
            return _make_starts(nLHS, nJit, d, LB, UB, span, z0, np.zeros(d, bool), rng)
    starts = makeStarts()
    initialguess = starts

    maxAttempts = 3
    attempt = 0
    flagg = -np.inf
    bestP, bestFval = None, np.inf
    startAgreement = np.nan
    P = z0.copy()
    fval = np.inf

    if nFree == 0:
        P = z0.copy()
        fval = f(z0)
        flagg = 1
        bestP, bestFval = z0, fval
        startAgreement = 1.0
    else:
        def fqs(U):
            return _bound(fq(U) / fscale)
        while (flagg <= 0 or not np.isfinite(fval)) and attempt < maxAttempts:
            attempt += 1
            if attempt > 1:
                starts = makeStarts()
                initialguess = np.vstack([initialguess, starts])
            runs = []
            for s in toQ(starts):
                x, fv, fl = local_solve(fqs, s)
                runs.append((x, _unbound(fv) * fscale, fl))
            fvs = np.array([r[1] for r in runs])
            flags = np.array([r[2] for r in runs])
            ib = int(np.argmin(fvs))
            P = fromQ(runs[ib][0])
            fval = fvs[ib]
            flagg = 1 if flags.any() else 0
            if np.isfinite(fvs[ib]) and fvs[ib] < bestFval:
                bestFval = fvs[ib]
                bestP = P.copy()
            okRun = (flags > 0) & np.isfinite(fvs)
            if okRun.any():
                startAgreement = float(np.sum(okRun & (fvs <= bestFval + 1e-4 * max(1, abs(bestFval))))
                                       / okRun.sum())

    if flagg > 0 and not np.isfinite(fval):
        flagg = 0
    if flagg <= 0:
        if bestP is not None:
            warnings.warn(f"fit_model:multistart - MultiStart did not converge cleanly after "
                          f"{maxAttempts} attempts; using best local minimum found (fval = {bestFval:g}).")
            P, fval = bestP, bestFval
        else:
            warnings.warn("fit_model:multistart - all attempts failed; falling back to the seed "
                          "parameters. Results from this fit are NOT reliable.")
            P = z0.copy()
            fval = f(z0)

    multistartExitflag = flagg
    converged = flagg > 0
    polishExitflag, polishImproved, confirmExitflag, polishSkipped = np.nan, False, np.nan, False

    rescale = {"run": False, "scaleBefore": fscale, "scaleAfter": fscale,
               "objectiveBefore": fval, "objectiveAfter": fval, "exitflag": np.nan, "improved": False}
    fscaleNew = max(abs(fval), 1) if np.isfinite(fval) else fscale
    if nFree > 0 and np.all(np.isfinite(P)) and np.isfinite(fval) and (
            fscaleNew < 1e-2 * fscale or fscaleNew > 1e2 * fscale):
        rescale.update(run=True, scaleAfter=fscaleNew)

        def fqsNew(U):
            return _bound(fq(U) / fscaleNew)
        try:
            qc, fc2, flc = local_solve(fqsNew, toQ(P)[0])
            fc2 = _unbound(fc2) * fscaleNew
            rescale["exitflag"] = flc
            if np.isfinite(fc2) and fc2 < fval:
                P, fval = fromQ(qc), fc2
                rescale["improved"] = True
                if flc > 0:
                    converged = True
            rescale["objectiveAfter"] = fval
        except Exception:
            pass
        fscale = fscaleNew

    def fqs_final(U):
        return _bound(fq(U) / fscale)

    def gradRef(fv):
        return max(abs(fv), 1) if np.isfinite(fv) else 1.0

    scaledGradInf, scaledGrad = np.nan, np.array([])
    if nFree > 0 and np.all(np.isfinite(P)):
        ref = gradRef(fval)
        scaledGrad, scaledGradInf = _scaled_gradient(lambda U: fq(U) / ref, toQ(P)[0], fdStep)

    if doPolish and nFree > 0 and np.all(np.isfinite(P)) and np.isfinite(scaledGradInf) \
            and scaledGradInf <= polishGradTol:
        polishSkipped = True
    if doPolish and nFree > 0 and np.all(np.isfinite(P)) and not polishSkipped:
        qP = toQ(P)[0]

        def fpol(qf):
            if np.any(qf < 0) or np.any(qf > 1):
                return np.inf
            return fq(qf)
        budget = int(math.ceil(polishBudget * nFree))
        for _pass in range(2):
            try:
                res = minimize(fpol, qP, method="Nelder-Mead",
                               options={"xatol": 1e-8, "fatol": 1e-8, "maxfev": budget,
                                        "maxiter": budget})
                polishExitflag = 1 if res.success else 0
                qfRaw = np.asarray(res.x, dtype=float)
                feasible = bool(np.all(qfRaw >= 0) and np.all(qfRaw <= 1))
                qTry = np.clip(qfRaw, 0, 1)
                fTry = fq(qTry)
                improvedNow = np.isfinite(fTry) and fTry < fval
                if improvedNow:
                    qP, P, fval = qTry, fromQ(qTry), fTry
                    polishImproved = True
                    converged = (polishExitflag == 1) and feasible
                elif polishExitflag == 1 and feasible:
                    converged = converged or flagg > 0
                if converged or not improvedNow:
                    break
            except Exception:
                break
        if polishImproved and not converged:
            try:
                qc, fc2, flc = local_solve(fqs_final, qP)
                fc2 = _unbound(fc2) * fscale
                if flc > 0 and np.isfinite(fc2) and fc2 <= fval + 1e-9 * max(1, abs(fval)):
                    qP, P, fval = qc, fromQ(qc), min(fc2, fval)
                    converged = True
                    confirmExitflag = flc
            except Exception:
                pass

    if polishImproved and nFree > 0 and np.all(np.isfinite(P)):
        ref = gradRef(fval)
        scaledGrad, scaledGradInf = _scaled_gradient(lambda U: fq(U) / ref, toQ(P)[0], fdStep)

    optimizerFval = fval
    fval, fitcurve, F1, IC = ctx.evaluate(P)
    residual = fitcurve - ydata

    diagnostics = {
        "numericalProfile": ode_profile(), "optimizerObjective": optimizerFval,
        "finalObjective": fval, "objectiveDifference": fval - optimizerFval,
        "parameters": P.copy(), "initialConditions": IC.copy(), "calibrationTimes": timevect,
        "numberOfObservations": ydata.size, "method": method1, "zeroReplacement": 0.001,
        "optimizerExitflag": flagg, "optimizationAttempts": attempt,
        "LB": LB.copy(), "UB": UB.copy(), "multistartExitflag": multistartExitflag,
        "polishExitflag": polishExitflag, "polishImproved": polishImproved,
        "polishSkipped": polishSkipped, "polishGradientTol": polishGradTol,
        "polishBudget": polishBudget, "confirmExitflag": confirmExitflag,
        "converged": bool(converged), "logscale": useLog, "freeIndex": free,
        "nFree": nFree, "objectiveScale": fscale, "rescalePass": rescale,
        "startAgreement": startAgreement, "scaledGradient": scaledGrad,
        "scaledGradientInfNorm": scaledGradInf,
    }

    if forecastingperiod < 1:
        forecastcurve = fitcurve
        timevect2 = timevect
        F2 = F1
    else:
        timevect2 = np.arange(data1[0, 0], data1[-1, 0] + forecastingperiod + 1) * DT
        F2 = solve_ode(model.fc, timevect2, IC, P, params.extra0)
        forecastcurve = observation_curve(F2, vars)

    return FitResult(P=P, residual=residual, fitcurve=fitcurve, forecastcurve=forecastcurve,
                     timevect2=timevect2, initialguess=initialguess, fval=fval, F1=F1, F2=F2,
                     diagnostics=diagnostics)
