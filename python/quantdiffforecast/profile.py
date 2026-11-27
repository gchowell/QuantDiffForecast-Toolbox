"""Profile likelihood (Run_ProfileLikelihood_ODEModel.m).

For every estimated quantity the profiler walks outward from the maximum-
likelihood estimate with a curvature-adapted step (Raue et al., 2009),
re-optimizes everything else with that coordinate pinned, and records the
deviance D(v) = -2 [log L(v) - log L(theta_hat)]. The (1 - alpha) interval is
{v : D(v) <= chi2inv(1 - alpha, 1)}; crossings are refined by
interpolation/bisection. The status rules, bound-activity checks, flatness
index and ridge slopes follow the MATLAB function.

MATLAB -> Python: ``opts`` struct -> keyword arguments; the returned struct
``PL`` -> a dict with the same field names (``PL['param'][k]['CI']`` ...);
``chi2inv`` -> ``scipy.stats.chi2.ppf``; paramIndex is 0-based.
"""
from __future__ import annotations

import os
import re
import warnings

import numpy as np
import pandas as pd
from scipy.stats import chi2

from .config import (Options, QuantDiffForecastError, check_initial_states, check_negative_data,
                     load_data, make_window, validate_options)
from .core import effective_n, get_nparams, ode_tolerances
from .optimizer import fit_model

DT = 1.0
_DEFAULTS = dict(paramIndex=None, numstartpoints=None, alpha=0.05, targetFirstD=None, h0=0.01,
                 maxPoints=40, maxRetract=6, maxFailures=3, stopFactor=2, refineIters=10,
                 refineTolD=0.05, refTol=0.1, flatTol=0.1, refitFactor=5, noiseRefits=3,
                 minCleanPts=5, boundTol=1e-3, maxRestarts=2, logSpacing=True, printscreen=True,
                 innerLogScale=True, innerScaleObjective=False, innerPolish=True, odeTol=None)


def _converged(dg):
    return bool(dg.get("converged", dg["optimizerExitflag"] > 0))


def _transforms(use_log):
    if use_log:
        return np.log10, lambda x: 10.0 ** x
    return (lambda v: v), (lambda x: x)


def _at_bound(P, jprof, fa, tol):
    hits = []
    nAlpha = fa["nModelAndX0"]
    for i in range(min(len(P), fa["nz"])):
        if i == jprof or fa["fixedz"][i] or fa["LBz"][i] == fa["UBz"][i]:
            continue
        atLB = P[i] <= fa["LBz"][i] + tol * abs(fa["LBz"][i])
        atUB = P[i] >= fa["UBz"][i] - tol * abs(fa["UBz"][i])
        if i == nAlpha:
            atLB = False
        if atLB or atUB:
            hits.append(fa["znames"][i])
    return hits


def _pinned_fit(j, v, seed, o, fa):
    pin = np.full(fa["nz"], np.nan)
    pin[j] = v
    seed = np.array(seed, float)
    seed[j] = v
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = fit_model(fa["data1"], seed, o["numstartpoints"], DT, fa["model"], fa["params"],
                          fa["vars"], 0, fa["method1"], pin=pin, rng=fa["rng"])
        ok = _converged(r.diagnostics) and np.isfinite(r.fval)
        return r.P, (r.fval if ok else np.nan), ok
    except Exception as exc:
        warnings.warn(f"ProfileLikelihood:fitFailed - {fa['znames'][j]} = {v:.4g} failed: {exc}")
        return seed, np.nan, False


def _walk(j, direction, Pseed, fmin, devfun, thr, o, fa, hStart):
    LBj, UBj, mle = fa["LBz"][j], fa["UBz"][j], Pseed[j]
    T, Tinv = _transforms(o["logSpacing"] and LBj > 0)
    width = abs(T(UBj) - T(LBj))
    x = T(mle)
    xb = T(UBj) if direction > 0 else T(LBj)
    rec = {"theta": [], "fval": [], "P": []}
    info = {"walkStop": "budget", "endpoint": "", "explored": mle, "nFailed": 0, "h0used": np.nan,
            "ci": np.nan, "endpointQuality": "", "boundActivity": [], "beyondHits": [],
            "flatToBound": False, "flatObserved": False, "bracket": None}
    if abs(xb - x) <= 1e-6 * width:
        info["walkStop"] = info["endpoint"] = "mleAtBound"
        return rec, info
    hmin, hmax = 1e-6 * width, width / 4
    h = hStart if np.isfinite(hStart) else o["h0"] * width
    h = min(max(h, hmin), hmax)
    info["h0used"] = h
    seed, Dprev, npts, calibrated, nfail, nretract = Pseed, 0.0, 0, False, 0, 0
    thetaPrev, reachedBox = mle, False
    while npts < o["maxPoints"]:
        xn = x + direction * h
        atEdge = direction * (xn - xb) >= 0
        if atEdge:
            xn = xb
        v = float(Tinv(xn))
        P, f, ok = _pinned_fit(j, v, seed, o, fa)
        npts += 1
        if not ok:
            nfail += 1
            info["nFailed"] += 1
            if nfail >= o["maxFailures"] or h <= hmin:
                info["walkStop"] = "solverFail"
                break
            h = max(h / 2, hmin)
            continue
        nfail = 0
        D = devfun(f, min(f, fmin))
        dD = D - Dprev
        if not calibrated and dD > 0 and not atEdge:
            calibrated = True
            hcal = h * np.sqrt(o["targetFirstD"] / dD)
            if dD > thr / 2 and hcal < h and nretract < o["maxRetract"]:
                h = min(max(hcal, hmin), hmax)
                nretract += 1
                calibrated = False
                continue
            h = min(max(hcal, hmin), hmax)
            info["h0used"] = h
        elif dD > thr / 2 and h > hmin and not atEdge and nretract < o["maxRetract"]:
            h = min(max(h * np.sqrt(o["targetFirstD"] / dD), hmin), h / 2)
            nretract += 1
            continue
        nretract = 0
        rec["theta"].append(v)
        rec["fval"].append(f)
        rec["P"].append(P[:fa["nz"]].copy())
        info["explored"] = v
        if info["bracket"] is None and Dprev < thr <= D:
            info["bracket"] = (thetaPrev, v)
        thetaPrev, x, Dprev, seed = v, xn, D, P
        if dD < thr / 10:
            h = min(2 * h, hmax)
        if D > o["stopFactor"] * thr:
            info["walkStop"] = "threshold"
            break
        if atEdge:
            reachedBox = True
            info["walkStop"] = "bound"
            break
    info["endpoint"] = "crossing" if info["bracket"] is not None else ("bound" if reachedBox else "")
    return rec, info


def _refine(rec, info, j, Phat, fmin, devfun, thr, o, fa):
    if info["endpoint"] != "crossing" or info["bracket"] is None:
        return rec, info
    T, Tinv = _transforms(o["logSpacing"] and fa["LBz"][j] > 0)
    allTheta = [Phat[j]] + rec["theta"]
    allF = [fmin] + rec["fval"]
    allP = [Phat[:fa["nz"]]] + rec["P"]
    a, b = info["bracket"]
    ia, ib = allTheta.index(a), allTheta.index(b)
    Da, Db = devfun(allF[ia], fmin), devfun(allF[ib], fmin)
    Pa, Pb = allP[ia], allP[ib]
    quality, Pjudge, ci = "approximate", Pb, np.nan
    width0 = abs(T(b) - T(a))
    lastSide = 0
    for _ in range(o["refineIters"]):
        lo, hi = min(T(a), T(b)), max(T(a), T(b))
        mg = 0.05 * (hi - lo)
        xc = T(a) + (thr - Da) * (T(b) - T(a)) / (Db - Da) if Db > Da else (lo + hi) / 2
        if xc < lo + mg or xc > hi - mg or lastSide == 2:
            xc = (lo + hi) / 2
            lastSide = 0
        c = float(Tinv(xc))
        Pc, fc, ok = _pinned_fit(j, c, Pa, o, fa)
        if not ok:
            info["nFailed"] += 1
            break
        rec["theta"].append(c)
        rec["fval"].append(fc)
        rec["P"].append(Pc[:fa["nz"]].copy())
        Dc = devfun(fc, min(fc, fmin))
        if abs(Dc - thr) < o["refineTolD"]:
            ci, Pjudge, quality = c, Pc, "refined"
            break
        if Dc >= thr:
            same = lastSide == 1
            b, Db, Pb, lastSide = c, Dc, Pc, 1
        else:
            same = lastSide == -1
            a, Da, Pa, lastSide = c, Dc, Pc, -1
        if same:
            lastSide = 2
        if abs(T(b) - T(a)) < 1e-3 * width0:
            ci = a + (thr - Da) * (b - a) / (Db - Da) if Db > Da else (a + b) / 2
            Pjudge = Pb if abs(ci - b) < abs(ci - a) else Pa
            quality = "refined"
            break
    if np.isnan(ci):
        ci = a + (thr - Da) * (b - a) / (Db - Da) if Db > Da else b
        Pjudge = Pb
    info["ci"] = ci
    info["endpointQuality"] = quality
    info["boundActivity"] = _at_bound(Pjudge, j, fa, o["boundTol"])
    info["beyondHits"] = _at_bound(Pb, j, fa, o["boundTol"])
    return rec, info


def _ridge(theta, D, Pg, clean, j, fa, flatTol):
    slopes = np.full(fa["nz"], np.nan)
    sel = clean & np.isfinite(D) & (D < flatTol) & (theta > 0)
    if sel.sum() < 3:
        return slopes, ""
    lt = np.log10(theta[sel])
    if lt.max() - lt.min() < 0.05:
        return slopes, ""
    parts = []
    for i in range(fa["nModelAndX0"]):
        if i == j or fa["fixedz"][i] or fa["LBz"][i] == fa["UBz"][i]:
            continue
        col = Pg[sel, i]
        if np.any(col <= 0):
            continue
        b = np.polyfit(lt, np.log10(col), 1)
        slopes[i] = b[0]
        if abs(b[0]) > 0.05:
            parts.append(f"{fa['znames'][i]} ~ {fa['znames'][j]}^{b[0]:.2f}")
    return slopes, "; ".join(parts)


def _reference_fit(data1, seed, nstart, o, fa, devfun, datac):
    def fit(s, n):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return fit_model(data1, s, n, DT, fa["model"], fa["params"], fa["vars"], 0,
                             fa["method1"], rng=fa["rng"])
    r = fit(seed, nstart)
    Phat, fmin, dg = r.P, r.fval, r.diagnostics
    if not np.isfinite(fmin):
        raise QuantDiffForecastError("ProfileLikelihood:referenceFit - non-finite objective.")
    if not _converged(dg):
        r2 = fit(Phat, 2 * max(nstart, 2))
        if np.isfinite(r2.fval) and r2.fval <= fmin:
            Phat, fmin, dg = r2.P, r2.fval, r2.diagnostics
    refConverged = _converged(dg)
    if not refConverged:
        warnings.warn("ProfileLikelihood:referenceFit - the reference fit is not reported as "
                      "converged; all statuses are reported as unresolved.")
    fv = [fmin]
    for _ in range(3):
        fv, improved = [fmin], False
        for _r in range(o["noiseRefits"]):
            r2 = fit(Phat, o["numstartpoints"])
            if _converged(r2.diagnostics) and np.isfinite(r2.fval):
                fv.append(r2.fval)
                if r2.fval < fmin and devfun(fmin, r2.fval) > 1e-6:
                    Phat, fmin, dg, improved = r2.P, r2.fval, r2.diagnostics, True
        if not improved:
            break
    refitVar = devfun(max(fv), min(fv))
    m = fa["method1"]
    if m in (0, 6) and fmin <= 1e-10 * np.sum(np.abs(datac) ** (2 - (m == 6))):
        raise QuantDiffForecastError("ProfileLikelihood:exactFit - the objective at the optimum is "
                                     "~0, so the profiled-scale deviance is undefined.")
    return Phat, fmin, dg, refitVar, refConverged


def _run_all(Phat, fmin, flatTolEff, devfun, thr, o, fa):
    betterP, betterF, out = None, fmin, []
    for j in o["paramIndex"]:
        mle = Phat[j]
        recR, infR = _walk(j, +1, Phat, fmin, devfun, thr, o, fa, np.nan)
        recL, infL = _walk(j, -1, Phat, fmin, devfun, thr, o, fa, infR["h0used"])
        recR, infR = _refine(recR, infR, j, Phat, fmin, devfun, thr, o, fa)
        recL, infL = _refine(recL, infL, j, Phat, fmin, devfun, thr, o, fa)
        theta = np.array(recL["theta"] + [mle] + recR["theta"], float)
        fval = np.array(recL["fval"] + [fmin] + recR["fval"], float)
        Pg = np.array(recL["P"] + [Phat[:fa["nz"]]] + recR["P"], float)
        order = np.argsort(theta, kind="stable")
        theta, fval, Pg = theta[order], fval[order], Pg[order]
        fin = np.isfinite(fval)
        if fin.any():
            ib = np.where(fin)[0][np.argmin(fval[fin])]
            if devfun(fmin, fval[ib]) > o["refTol"] and fval[ib] < betterF:
                betterF, betterP = fval[ib], Pg[ib].copy()
        D = np.array([devfun(f, fmin) if np.isfinite(f) else np.nan for f in fval])
        clean = np.array([np.isfinite(D[g]) and not _at_bound(Pg[g], j, fa, o["boundTol"])
                          for g in range(theta.size)])
        maxD = float(np.nanmax(D))
        maxDclean = float(np.max(D[clean])) if clean.any() else np.nan
        sides = [infL, infR]
        for s, inf1 in enumerate(sides):
            if inf1["endpoint"] == "crossing" and (inf1["boundActivity"] or inf1["beyondHits"]):
                onSide = theta < mle if s == 0 else theta > mle
                sel = (onSide & clean & np.isfinite(D) & (theta >= min(inf1["ci"], mle))
                       & (theta <= max(inf1["ci"], mle)))
                flat = not np.any(D[sel] >= flatTolEff)
                if flat and not inf1["boundActivity"]:
                    inf1["boundActivity"] = inf1["beyondHits"]
                inf1["flatToBound"] = flat and bool(inf1["boundActivity"])
                inf1["flatObserved"] = flat and sel.sum() >= 2
        CI, explored = [np.nan, np.nan], [np.nan, np.nan]
        stop, endQ, act, note = ["", ""], ["", ""], [[], []], ["", ""]
        for s, inf1 in enumerate(sides):
            explored[s] = inf1["explored"]
            act[s] = inf1["boundActivity"]
            ep = inf1["endpoint"]
            if ep == "crossing":
                stop[s], CI[s], endQ[s] = "crossing", inf1["ci"], inf1["endpointQuality"]
                if inf1["walkStop"] in ("solverFail", "budget"):
                    note[s] = f"exploration beyond the crossing ended in {inf1['walkStop']}"
            elif ep == "bound":
                stop[s], CI[s], endQ[s] = "bound", inf1["explored"], "search bound"
            elif ep == "mleAtBound":
                stop[s], CI[s], endQ[s] = "mleAtBound", mle, "search bound"
            else:
                stop[s], endQ[s] = inf1["walkStop"], "unresolved"
        unresolved = [x in ("solverFail", "budget") for x in stop]
        crossing = [x == "crossing" for x in stop]
        atBox = [x in ("bound", "mleAtBound") for x in stop]
        activity = [bool(a) for a in act]
        names, reasons = ("lower", "upper"), []
        for s in range(2):
            if stop[s] == "solverFail":
                reasons.append(f"{names[s]}: solver failure at {explored[s]:.4g}")
            elif stop[s] == "budget":
                reasons.append(f"{names[s]}: fit budget exhausted at {explored[s]:.4g}")
            elif stop[s] == "mleAtBound":
                reasons.append(f"{names[s]}: boundary MLE")
            elif stop[s] == "bound":
                reasons.append(f"{names[s]}: open at search bound")
            elif stop[s] == "crossing":
                if activity[s]:
                    reasons.append(f"{names[s]}: {','.join(act[s])} on bound"
                                   + (" with flat profile up to the crossing"
                                      if sides[s]["flatObserved"] else ""))
                if endQ[s] == "approximate":
                    reasons.append(f"{names[s]}: approximate endpoint")
        why = "; ".join(reasons)
        nclean = int(clean.sum())
        if any(unresolved):
            status = f"unresolved ({why})"
        elif all(crossing) and not any(activity):
            status = "finite interval" + (f" ({why})" if "approximate" in endQ else "")
        elif (all(atBox[s] or (crossing[s] and activity[s] and sides[s]["flatToBound"])
                  for s in range(2)) and nclean >= o["minCleanPts"] and maxDclean < flatTolEff):
            status = "numerically flat over the explored range"
        elif all(crossing):
            status = f"finite interval, nuisance-bound activity ({why}) - sensitivity check needed"
        else:
            if nclean < o["minCleanPts"]:
                reasons.append("few bound-free points")
            elif maxDclean < flatTolEff:
                reasons.append("bound-free profile flat")
            status = f"bound-limited ({'; '.join(reasons)})"
        slopes, relation = _ridge(theta, D, Pg, clean, j, fa, flatTolEff)
        out.append({"label": fa["znames"][j], "index": j, "mle": mle, "CI": np.array(CI),
                    "explored": np.array(explored), "stopReason": stop, "endpointQuality": endQ,
                    "explorationNote": note, "boundActivity": act, "status": status,
                    "FI": maxDclean / thr, "ridge": slopes, "candidate_ridge_relation": relation,
                    "theta": theta, "D": D, "fval": fval, "P": Pg, "maxD": maxD,
                    "maxDclean": maxDclean, "nFailed": infL["nFailed"] + infR["nFailed"],
                    "h0used": infR["h0used"]})
    return out, betterP, betterF


def run_profile_likelihood(opts: Options, tstart1=None, windowsize1=None, *, seed=None, rng=None,
                           verbose=True, save=True, **kwargs):
    """Run_ProfileLikelihood_ODEModel(options, tstart1, windowsize1, opts).

    Keyword arguments are the fields of the MATLAB ``opts`` struct, e.g.
    ``paramIndex=[0]`` (0-based slots of [params, X0, alpha, d]), ``alpha=0.05``,
    ``numstartpoints=3``. Returns the dict ``PL``.
    """
    o = dict(_DEFAULTS)
    unknown = set(kwargs) - set(o)
    if unknown:
        warnings.warn(f"ProfileLikelihood:unknownOption - ignoring {sorted(unknown)}")
    o.update({k: v for k, v in kwargs.items() if k in o and v is not None})
    validate_options(opts, "Run_ProfileLikelihood_ODEModel")
    if o["numstartpoints"] is None:
        o["numstartpoints"] = min(max(2, opts.numstartpoints), 3)
    opts = opts.copy()
    params, vars, model = opts.params, opts.vars, opts.model
    method1, dist1 = opts.method1, opts.effective_dist1()
    tstart1 = opts.tstart1 if tstart1 is None else tstart1
    windowsize1 = opts.windowsize1 if windowsize1 is None else windowsize1
    params.logscale = bool(o["innerLogScale"])
    params.scaleObjective = bool(o["innerScaleObjective"])
    params.polish = bool(o["innerPolish"])
    rng = np.random.default_rng(seed) if rng is None else rng
    tol = {} if o["odeTol"] is None else {"rtol": o["odeTol"][0], "atol": o["odeTol"][1]}

    data = load_data(opts)
    if tstart1 + windowsize1 - 1 > data.shape[0]:
        raise QuantDiffForecastError("Calibration window exceeds the available rows.")
    check_negative_data(method1, data, range(tstart1, tstart1 + windowsize1), "ProfileLikelihood")
    i0 = tstart1 - 1
    datac = data[i0:i0 + windowsize1, 1:]
    data1 = data[i0:i0 + windowsize1, :]
    I0 = datac[0]
    nfit = len(vars.fit_index)
    vars.window = make_window(data, tstart1, windowsize1, DT)
    check_initial_states(vars, tstart1, "ProfileLikelihood")
    nEff, _ = effective_n(datac.size, nfit, params.fixI0)
    params0 = np.concatenate([params.initial, I0 if params.fixI0 else vars.initial[vars.fit_index],
                              [1, 1]])
    if get_nparams(method1, dist1, int(np.sum(params.fixed == 0)),
                   nfit * (params.fixI0 == 0)) >= nEff:
        raise QuantDiffForecastError("Number of estimated parameters must be smaller than the "
                                     "number of informative data points.")
    znames = list(params.label) + [f"X0({i + 1})" for i in vars.fit_index] + ["alpha", "d"]
    nz = len(znames)
    nModelAndX0 = params.num + nfit
    if method1 in (1, 3, 4, 5):
        devfun = lambda f, f0: 2 * (f - f0)
        objname = "negative log-likelihood"
    elif method1 == 0:
        devfun = lambda f, f0: nEff * np.log(f / f0)
        objname = "sum of squared residuals (Gaussian, variance profiled out)"
    else:
        devfun = lambda f, f0: 2 * nEff * np.log(f / f0)
        objname = "sum of absolute deviations (Laplace, scale profiled out)"
    thr = float(chi2.ppf(1 - o["alpha"], 1))
    if o["targetFirstD"] is None:
        o["targetFirstD"] = thr / 4
    if verbose:
        print(f"Profile likelihood: {model.name}, rows {tstart1}-{tstart1 + windowsize1 - 1}, "
              f"method {method1}; {100 * (1 - o['alpha']):.0f}% threshold = {thr:.3f}")

    fa = {"data1": data1, "model": model, "params": params, "vars": vars, "method1": method1,
          "znames": znames, "nModelAndX0": nModelAndX0, "nz": nz, "rng": rng}
    with ode_tolerances(**tol):
        Phat, fmin, dg, refitVar, refConv = _reference_fit(data1, params0, opts.numstartpoints, o,
                                                           fa, devfun, datac)
        LBz, UBz = dg["LB"][:nz], dg["UB"][:nz]
        fixedz = np.concatenate([params.fixed == 1, np.zeros(nfit + 2, bool)])
        if params.fixI0:
            fixedz[params.num:params.num + nfit] = True
        fixedz |= LBz == UBz
        if o["paramIndex"] is None:
            o["paramIndex"] = [int(i) for i in np.where(~fixedz[:nModelAndX0])[0]]
        o["paramIndex"] = [int(i) for i in np.atleast_1d(o["paramIndex"])]
        if any(fixedz[i] for i in o["paramIndex"]):
            raise QuantDiffForecastError("paramIndex must point to estimated slots of "
                                         "[params, X0, alpha, d] (0-based).")
        fa.update(LBz=LBz, UBz=UBz, fixedz=fixedz)
        quality = "ok"
        flatTolEff = max(o["flatTol"], o["refitFactor"] * refitVar)
        if flatTolEff > thr / 2:
            warnings.warn("ProfileLikelihood:refitVariability - repeated fits at the MLE disagree; "
                          "raise numstartpoints. Numerical quality flagged as poor.")
            flatTolEff, quality = thr / 2, "poor"
        restarts, refStatus = 0, "validated"
        bestKnown = {"P": Phat, "fval": fmin}
        while True:
            paramOut, betterP, betterF = _run_all(Phat, fmin, flatTolEff, devfun, thr, o, fa)
            if betterP is None:
                if not refConv:
                    for p in paramOut:
                        p["status"] = "unresolved (reference fit not converged)"
                break
            bestKnown = {"P": betterP, "fval": betterF}
            if restarts >= o["maxRestarts"]:
                refStatus = "unresolved"
                warnings.warn("ProfileLikelihood:restartLimit - a constrained fit still found a "
                              "lower objective than the reference; statuses are unresolved.")
                for p in paramOut:
                    p["status"] = "unresolved (reference optimum not validated)"
                break
            restarts += 1
            seedP = np.concatenate([betterP, Phat[nz:]]) if betterP.size < Phat.size else betterP
            Phat, fmin, dg, refitVar, refConv = _reference_fit(data1, seedP, opts.numstartpoints, o,
                                                               fa, devfun, datac)
            flatTolEff = min(max(o["flatTol"], o["refitFactor"] * refitVar), thr / 2)
    if not refConv and refStatus == "validated":
        refStatus = "unresolved (reference fit not converged)"

    suffix = re.sub(r'[\s/\\:*?"<>|@]+', "_",
                    f"model_name-{model.name}-opts-{opts.name}-data-"
                    f"{os.path.splitext(opts.cadfilename_txt())[0]}-fixI0-{params.fixI0}-method-"
                    f"{method1}-dist-{dist1}-tstart-{tstart1}-calibrationperiod-{windowsize1}-level-"
                    f"{100 * (1 - o['alpha']):g}-{opts.caddisease}-{opts.datatype}.csv")
    PL = {"model": model.name, "method1": method1, "dist1": dist1, "tstart1": tstart1,
          "windowsize1": windowsize1, "threshold": thr, "alpha": o["alpha"], "objective": objname,
          "fmin": fmin, "Phat": Phat, "znames": znames, "LBz": LBz, "UBz": UBz, "referenceFit": dg,
          "referenceStatus": refStatus, "bestKnown": bestKnown, "refitVariability": refitVar,
          "numericalQuality": quality, "flatTolEffective": flatTolEff, "restarts": restarts,
          "mleAtBound": _at_bound(Phat, -1, fa, o["boundTol"]), "optionsName": opts.name,
          "opts": o, "param": paramOut, "suffix": suffix}
    rows = []
    for p in paramOut:
        rows.append({"parameter": p["label"], "index": p["index"], "MLE": p["mle"],
                     "PL_CI_LB": p["CI"][0], "PL_CI_UB": p["CI"][1],
                     "explored_LB": p["explored"][0], "explored_UB": p["explored"][1],
                     "stop_LB": p["stopReason"][0], "stop_UB": p["stopReason"][1],
                     "endpoint_LB": p["endpointQuality"][0], "endpoint_UB": p["endpointQuality"][1],
                     "LB_bound_activity": ",".join(p["boundActivity"][0]),
                     "UB_bound_activity": ",".join(p["boundActivity"][1]),
                     "maxD": p["maxD"], "maxD_boundfree": p["maxDclean"], "FI": p["FI"],
                     "candidate_ridge_relation": p["candidate_ridge_relation"],
                     "failed_fits": p["nFailed"], "status": p["status"]})
        if verbose:
            print(f"  {p['label']}: MLE {p['mle']:.4g}, {100 * (1 - o['alpha']):.0f}% PL interval "
                  f"[{p['CI'][0]:.4g}, {p['CI'][1]:.4g}], FI = {p['FI']:.2g}, {len(p['theta'])} "
                  f"points -> {p['status']}")
    PL["summary"] = pd.DataFrame(rows)
    if save:
        outdir = opts.outdir()
        PL["summary"].to_csv(os.path.join(outdir, "PLsummary-" + suffix), index=False)
        for p in paramOut:
            T = pd.DataFrame(np.column_stack([p["theta"], p["D"], p["fval"], p["P"]]),
                             columns=["theta", "D", "objective"] + znames)
            T.to_csv(os.path.join(outdir, f"PLprofile-param-{p['index'] + 1}-" + suffix), index=False)
        if o["printscreen"]:
            from .plotting import plot_profiles
            plot_profiles(PL, outdir)
    return PL
