"""Configuration objects: the Python replacement for the MATLAB options files.

In MATLAB an analysis is configured by an options *function* that returns
13 (fit-type) or 15 (forecast-type) positional outputs and, as a side effect,
sets the global variable ``method1``. In Python a single :class:`Options`
dataclass holds everything, including ``method1``, so there are no globals and
no output-order rules to remember. Fitting runners simply ignore the
forecasting fields.

Indexing conventions (the one place where the port deliberately differs from
MATLAB):

* positions inside arrays are 0-based: ``vars.fit_index = [4]`` selects the
  fifth state (``C`` in the SEIR example), ``p[0]`` is the first parameter in
  a model or ``initialFcn``;
* window starts ``tstart1``/``tend1`` stay 1-based *row numbers* of the data
  file, exactly as in the MATLAB toolbox and its README ("rows 1-17");
* output file names reproduce the MATLAB names, so ``fit_index`` is written
  1-based there (``...-fit_index-5-...``) and files can be compared directly.
"""
from __future__ import annotations

import copy
import os
from dataclasses import dataclass, field
from types import SimpleNamespace
from typing import Any, Callable, Optional, Sequence

import numpy as np

PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(PACKAGE_DIR)
DEFAULT_INPUT_DIR = os.path.join(PROJECT_DIR, "input")
DEFAULT_OUTPUT_DIR = os.path.join(PROJECT_DIR, "output")

SUPPORTED_METHODS = (0, 1, 3, 4, 5, 6)


class QuantDiffForecastError(ValueError):
    """Raised for invalid configurations or data (MATLAB error identifiers
    are kept in the message, e.g. 'QuantDiffForecast:UnsupportedMethod')."""


@dataclass
class Model:
    fc: Callable                      # dx = fc(t, x, p, extra0)
    name: str = "model"


@dataclass
class Params:
    label: Sequence[str]
    LB: Sequence[float]
    UB: Sequence[float]
    initial: Sequence[float]
    fixed: Sequence[int]
    fixI0: int = 1
    composite: Optional[Callable] = None
    composite_name: str = ""
    extra0: Any = None
    # optimizer settings (all optional, same meaning as in fit_model.m)
    logscale: bool = False
    scaleObjective: bool = True
    optimOptions: Optional[dict] = None
    useParallel: bool = False         # kept for compatibility; use n_jobs in the runners
    polish: bool = False
    polish_gradient_tol: float = 1e-2
    polish_budget: float = 100
    PIIscale: Optional[Sequence[float]] = None

    def __post_init__(self):
        self.label = list(self.label)
        self.LB = np.atleast_1d(np.asarray(self.LB, dtype=float)).copy()
        self.UB = np.atleast_1d(np.asarray(self.UB, dtype=float)).copy()
        self.initial = np.atleast_1d(np.asarray(self.initial, dtype=float)).copy()
        self.fixed = np.atleast_1d(np.asarray(self.fixed, dtype=int)).copy()

    @property
    def num(self) -> int:
        return len(self.label)


@dataclass
class Vars:
    label: Sequence[str]
    initial: Sequence[float]
    fit_index: Sequence[int]          # 0-based state positions
    fit_diff: Sequence[int]
    initialFcn: Optional[Callable] = None   # IC = initialFcn(p, x0, vars)
    window: Any = None                # set by the runners (see make_window)

    def __post_init__(self):
        self.label = list(self.label)
        self.initial = np.atleast_1d(np.asarray(self.initial, dtype=float)).copy()
        self.fit_index = np.atleast_1d(np.asarray(self.fit_index, dtype=int)).copy()
        self.fit_diff = np.atleast_1d(np.asarray(self.fit_diff, dtype=int)).copy()

    @property
    def num(self) -> int:
        return len(self.label)


@dataclass
class Options:
    """All settings of one analysis (fit-type and forecast-type in one)."""
    cadfilename1: str
    caddisease: str
    datatype: str
    method1: int
    dist1: int
    numstartpoints: int
    B: int
    model: Model
    params: Params
    vars: Vars
    windowsize1: int
    tstart1: int = 1
    tend1: int = 1
    printscreen1: bool = True
    getperformance: bool = True
    forecastingperiod: int = 0
    input_dir: Optional[str] = None
    output_dir: Optional[str] = None
    name: str = ""                    # label of the options set (like func2str)

    def copy(self, **changes) -> "Options":
        o = copy.deepcopy(self)
        for k, v in changes.items():
            setattr(o, k, v)
        return o

    # ---------------------------------------------------------------
    def effective_dist1(self) -> int:
        """For positive method codes the runners set dist1 = method1."""
        return self.method1 if self.method1 > 0 else self.dist1

    def input_path(self) -> str:
        fname = self.cadfilename1
        if not fname.lower().endswith(".txt"):
            fname = fname + ".txt"
        return os.path.join(self.input_dir or DEFAULT_INPUT_DIR, fname)

    def cadfilename_txt(self) -> str:
        f = self.cadfilename1
        return f if f.lower().endswith(".txt") else f + ".txt"

    def outdir(self) -> str:
        d = self.output_dir or DEFAULT_OUTPUT_DIR
        os.makedirs(d, exist_ok=True)
        return d


# ---------------------------------------------------------------------------
# validation helpers (quantdiffCheckMethod, quantdiffCheckNegativeData, ...)
# ---------------------------------------------------------------------------
def check_method(method1, caller="QuantDiffForecast"):
    if method1 is None:
        raise QuantDiffForecastError(
            f"QuantDiffForecast:UnsupportedMethod - {caller}: method1 is empty.")
    if method1 not in SUPPORTED_METHODS:
        hint = ""
        if method1 == 2:
            hint = (" For least squares with negative binomial bootstrap noise "
                    "(VAR = factor*mean) set method1 = 0 and dist1 = 2.")
        raise QuantDiffForecastError(
            f"QuantDiffForecast:UnsupportedMethod - {caller}: method1 = {method1} is "
            "not supported; use 0 (least squares), 1 (Poisson), 3, 4 or 5 "
            f"(negative binomial) or 6 (SAD/Laplace).{hint}")


def check_negative_data(method1, data, rows_1based, caller="QuantDiffForecast"):
    """Count likelihoods cannot fit negative observations."""
    if method1 not in (1, 3, 4, 5):
        return
    rows = [r for r in rows_1based if 1 <= r <= data.shape[0]]
    bad = []
    for r in rows:
        for c in range(1, data.shape[1]):
            if data[r - 1, c] < 0:
                bad.append(f"  row {r} (time {data[r-1,0]:g}), series {c}: {data[r-1,c]:g}")
    if bad:
        mname = "Poisson" if method1 == 1 else "negative binomial"
        more = f"\n  ... and {len(bad) - 10} more" if len(bad) > 10 else ""
        raise QuantDiffForecastError(
            f"QuantDiffForecast:NegativeCounts - {caller}: the calibration data contain "
            f"{len(bad)} negative value(s), which the {mname} likelihood (method1 = {method1}) "
            "cannot fit.\n" + "\n".join(bad[:10]) + more +
            "\nCorrect the data before fitting. Least squares (method1 = 0) and SAD "
            "(method1 = 6) accept negative values.")


def check_initial_states(vars: Vars, last_window_start, caller="QuantDiffForecast"):
    """Warn when later windows restart unfitted states from vars.initial."""
    import warnings
    flagged = (last_window_start > 1 and len(vars.initial) > len(vars.fit_index)
               and vars.initialFcn is None)
    if flagged:
        warnings.warn(
            f"QuantDiffForecast:rollingInitialStates - {caller}: calibration windows start "
            f"as late as row {last_window_start}, but the "
            f"{len(vars.initial) - len(vars.fit_index)} state(s) that are not fitted start "
            "every window from vars.initial. Set vars.initialFcn so the initial "
            "conditions follow the window (vars.window.priorCum holds the counts observed "
            "before it), or use windows that all start at row 1.", stacklevel=2)
    return flagged


def validate_options(opts: Options, caller="QuantDiffForecast"):
    p = opts.params
    check_method(opts.method1, caller)
    n = p.num
    if not (len(p.fixed) == n and len(p.LB) == n and len(p.UB) == n and len(p.initial) == n):
        raise QuantDiffForecastError(
            "one or more parameter vectors do not match the number of parameters (len(params.label))")
    for j in range(n):
        if p.initial[j] < p.LB[j] or p.initial[j] > p.UB[j]:
            raise QuantDiffForecastError(
                "values in params.initial should lie within their parameter bounds "
                "defined by params.LB and params.UB")
    v = opts.vars
    if len(v.fit_index) != len(v.fit_diff):
        raise QuantDiffForecastError("vars.fit_index and vars.fit_diff must have the same length.")
    if np.any(v.fit_index < 0) or np.any(v.fit_index >= v.num):
        raise QuantDiffForecastError("vars.fit_index must hold 0-based state positions.")


def load_data(opts: Options) -> np.ndarray:
    path = opts.input_path()
    if not os.path.isfile(path):
        raise FileNotFoundError(f'File "{path}" not found in the specified directory.')
    data = np.loadtxt(path, ndmin=2)
    if data.size == 0:
        raise QuantDiffForecastError("The dataset is empty")
    nser = data.shape[1] - 1
    if not (len(opts.vars.fit_index) == len(opts.vars.fit_diff) == nser):
        raise QuantDiffForecastError(
            "The number of state variables in the data file should be consistent with "
            "the dimensions of vars.fit_index and vars.fit_diff")
    return data


def make_window(data: np.ndarray, tstart_row: int, windowsize: int, DT: float = 1.0):
    """vars.window for the calibration window starting at 1-based row tstart_row."""
    i0 = tstart_row - 1
    rows = np.arange(tstart_row, tstart_row + windowsize)
    datac = data[i0:i0 + windowsize, 1:]
    return SimpleNamespace(
        tstart=tstart_row, rows=rows, time0=DT * data[i0, 0],
        firstObs=datac[0, :].copy(),
        priorCum=data[:i0, 1:].sum(axis=0) if i0 > 0 else np.zeros(data.shape[1] - 1))
