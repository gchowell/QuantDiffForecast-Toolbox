"""Simulation and parameter-recovery studies.

Python port of Run_simulate_ODEModel.m, Run_PracticalIdentifiability_ODEModel.m
(with Run_Forecasting_ODEModel2.m) and summarizeIdentifiability.m.
"""
from __future__ import annotations

import json
import os
import warnings

import numpy as np
import pandas as pd

from .config import DEFAULT_INPUT_DIR, Options, validate_options
from .core import initial_conditions, solve_ode
from .noise import add_error_structure
from .runners import _run_windows

_DEFAULT_FACTOR = {0: 4, 1: 1, 2: 3, 3: 3, 4: 0.1, 5: (0.5, 1.5), 6: 20}


def _parse_noise(dist1, spec):
    d = 1.0
    spec = np.atleast_1d(np.asarray(spec, float))
    if dist1 == 5:
        if spec.size == 2:
            factor1, d = float(spec[0]), float(spec[1])
        elif spec.size == 1:
            factor1 = float(spec[0])
            warnings.warn("Run_simulate_ODEModel:dEqualsOne - dist1=5 with a scalar factor1 "
                          "simulates with d = 1 (identical to dist1=3). Pass (alpha, d).")
        else:
            raise ValueError("For dist1=5 factor1 must be a scalar alpha or the pair (alpha, d).")
        if not (0.6 <= d <= 5):
            raise ValueError(f"d = {d} lies outside the estimation box [0.6, 5] of method 5.")
    else:
        if spec.size != 1:
            raise ValueError(f"factor1 must be a scalar for dist1 = {dist1}.")
        factor1 = float(spec[0])
    if not np.isfinite(factor1) or factor1 < 0:
        raise ValueError("factor1 must be a finite nonnegative scalar.")
    return factor1, d


def run_simulate(opts: Options, windowsize1=None, factor1=None, *, seed=None, rng=None):
    """Run_simulate_ODEModel: simulate one noisy dataset at params.initial.

    Returns (states, curves, noiseTruth): the noise-free states (rows = time),
    the matrix [time, noisy observations...] and the generating noise settings.
    The series is windowsize1 + max(30, forecastingperiod) rows long.
    """
    validate_options(opts, "Run_simulate_ODEModel")
    rng = np.random.default_rng(seed) if rng is None else rng
    dist1 = opts.effective_dist1()
    if factor1 is None:
        factor1 = _DEFAULT_FACTOR[dist1]
    factor1, d = _parse_noise(dist1, factor1)
    noiseTruth = {"dist1": dist1, "factor1": factor1, "d": d}
    nExtra = max(30, opts.forecastingperiod)
    n = (windowsize1 or opts.windowsize1) + nExtra
    p, v = opts.params, opts.vars
    timevect = np.arange(n, dtype=float)
    z = np.concatenate([p.initial, v.initial[v.fit_index], [1, 1]])
    IC = initial_conditions(z, p, v)
    F = solve_ode(opts.model.fc, timevect, IC, p.initial, p.extra0, rtol=1e-8, atol=1e-10)
    cols = [timevect]
    for idx, dif in zip(v.fit_index, v.fit_diff):
        c = F[:, idx]
        if dif == 1:
            c = np.abs(np.concatenate(([c[0]], np.diff(c))))
        noisy = add_error_structure(np.cumsum(c), 1, dist1, factor1, d, rng)[:, 0]
        if not np.all(np.isfinite(noisy)):
            raise RuntimeError("Run_simulate_ODEModel:invalidNoise - non-finite noisy values.")
        cols.append(np.maximum(noisy, 0))
    return F, np.column_stack(cols), noiseTruth


def run_practical_identifiability(opts: Options, windowsize1=None, factor1=None, replicates=1,
                                  replicate_ids=None, output_dir=None, n_jobs=1, verbose=True):
    """Run_PracticalIdentifiability_ODEModel: simulate replicate datasets with
    known parameters (params.initial), refit each one with the full bootstrap,
    and save one results file per replicate.

    Replicate r uses seed r for the simulation and the fit, so replicates are
    paired across window lengths, as in MATLAB. Returns the list of per-
    replicate result dicts (also saved as results_W<W>_rep<r>.npz/.json).
    """
    W = windowsize1 or opts.windowsize1
    ids = list(replicate_ids) if replicate_ids is not None else list(range(1, replicates + 1))
    input_dir = opts.input_dir or DEFAULT_INPUT_DIR
    os.makedirs(input_dir, exist_ok=True)
    out_dir = output_dir or opts.outdir()
    os.makedirs(out_dir, exist_ok=True)
    cadbase = os.path.splitext(os.path.basename(opts.cadfilename1))[0]
    results = []
    for r in ids:
        rng = np.random.default_rng(r)
        _, curves, noiseTruth = run_simulate(opts, W, factor1, rng=rng)
        fname = f"{cadbase}-windowsize-{W}-replicate-{r}.txt"
        np.savetxt(os.path.join(input_dir, fname), curves, fmt="%.7e")
        with open(os.path.join(input_dir, fname.replace(".txt", "-noise.json")), "w") as fh:
            json.dump({"noiseTruth": noiseTruth, "windowsize1": W, "replicate": r, "seed": r}, fh)
        o = opts.copy(cadfilename1=fname, output_dir=out_dir, input_dir=input_dir)
        res = _run_windows(o, 1, 1, W, opts.forecastingperiod, bool(opts.getperformance),
                           "Run_PracticalIdentifiability_ODEModel", rng=rng, n_jobs=n_jobs,
                           verbose=False, plot=False, save=False)
        w = res.windows[0]
        rec = {"W": W, "replicate": r, "P": w.P, "param_estims": w.param_estims,
               "param_bootmedian": w.param_bootmedian, "converged": bool(w.fit_diagnostics["converged"]),
               "noiseTruth": noiseTruth, **{k: w.pii[k] for k in
                                            ("SCI", "PII_rangeFrac", "PII_scale", "fracLB", "fracUB",
                                             "MCSE_PII")}}
        np.savez_compressed(os.path.join(out_dir, f"results_W{W}_rep{r}.npz"),
                            **{k: np.asarray(val) for k, val in rec.items() if k != "noiseTruth"},
                            noiseTruth=json.dumps(noiseTruth))
        results.append(rec)
        if verbose:
            print(f"[practical identifiability] W = {W}, replicate {r}: done")
    return results


def summarize_identifiability(output_dir, windowsizes, num_replicates, estimated_indices,
                              true_values, labels, save=True):
    """summarizeIdentifiability: recovery metrics per window length and slot.

    estimated_indices are 0-based slots into [params, X0..., alpha, d].
    Returns a DataFrame (also written to recovery_summary.csv).
    """
    rows = []
    for W in np.atleast_1d(windowsizes):
        recs = []
        for r in range(1, num_replicates + 1):
            path = os.path.join(output_dir, f"results_W{int(W)}_rep{r}.npz")
            if os.path.isfile(path):
                with np.load(path) as z:
                    recs.append({k: z[k] for k in z.files})
        ok = [rec for rec in recs if bool(rec["converged"])]
        for slot, truth, lab in zip(estimated_indices, true_values, labels):
            row = {"windowsize": int(W), "parameter": lab, "truth": truth, "files": len(recs),
                   "converged": len(ok)}
            if ok:
                mle = np.array([rec["param_estims"][slot, 0] for rec in ok])
                lo = np.array([rec["param_estims"][slot, 1] for rec in ok])
                hi = np.array([rec["param_estims"][slot, 2] for rec in ok])
                bm = np.array([rec["param_bootmedian"][slot] for rec in ok])
                cov = float(np.mean((lo <= truth) & (truth <= hi)))

                def col(key):
                    return np.array([rec[key][slot] if slot < rec[key].size else np.nan for rec in ok])
                sci, rf = col("SCI"), col("PII_rangeFrac")
                hits = col("fracLB") + col("fracUB")
                row.update(bias_MLE=mle.mean() - truth,
                           relbias_MLE=(mle.mean() - truth) / truth if truth != 0 else np.nan,
                           bias_bootmedian=bm.mean() - truth,
                           RMSE_MLE=float(np.sqrt(np.mean((mle - truth) ** 2))),
                           coverage95=cov, coverage_MCSE=float(np.sqrt(cov * (1 - cov) / len(ok))),
                           PII_scale=float(np.nanmedian(col("PII_scale"))),
                           median_PII=float(np.nanmedian(sci)) if np.any(np.isfinite(sci)) else np.nan,
                           median_MCSE_PII=float(np.nanmedian(col("MCSE_PII")))
                           if np.any(np.isfinite(col("MCSE_PII"))) else np.nan,
                           n_PII_undefined=int(np.sum(~np.isfinite(sci))),
                           median_rangeFrac=float(np.nanmedian(rf)) if np.any(np.isfinite(rf)) else np.nan,
                           frac_bound_hits=float(np.mean(np.nan_to_num(hits) > 0.05)))
            rows.append(row)
    S = pd.DataFrame(rows)
    if save:
        S.to_csv(os.path.join(output_dir, "recovery_summary.csv"), index=False)
    return S
