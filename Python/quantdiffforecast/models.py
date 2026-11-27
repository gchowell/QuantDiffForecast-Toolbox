"""Bundled ODE models and composite (derived-quantity) functions.

Python port of SEIR1.m, EXP.m, GGM.m, GLM.m, SEIR_unreported.m, R0s.m,
R0s_SEIR.m and R0s_unreported.m.

Every model has the MATLAB signature, with 0-based indexing:

    dx = model(t, x, p, extra0)

    t       time (scalar)
    x       state vector, 1-D numpy array, in the order of vars.label
    p       full parameter vector, in the order of params.label
    extra0  params.extra0 (anything; None when unused)

and returns a 1-D array of derivatives in state order. A composite function
receives a 2-D array with one row per parameter draw and returns one value per
row, exactly like the MATLAB versions.
"""
from __future__ import annotations

import math

import numpy as np


def SEIR1(t, x, p, extra0=None):
    """SEIR model with a cumulative-infections state.

    p = [beta, kappa, gamma, N]; states S, E, I, R, C.
    """
    beta, k, gamma, N = p[0], p[1], p[2], p[3]
    S, E, I = x[0], x[1], x[2]
    inf = beta * S * I / N
    return np.array([-inf, inf - k * E, k * E - gamma * I, gamma * I, k * E])


def EXP(t, x, p, extra0=None):
    """Exponential growth dC/dt = r C. p = [r]."""
    return np.array([p[0] * x[0]])


def GGM(t, x, p, extra0=None):
    """Generalized growth model dC/dt = r C^p. p = [r, p]."""
    return np.array([p[0] * x[0] ** p[1]])


def GLM(t, x, p, extra0=None):
    """Generalized logistic model dC/dt = r C^p (1 - C/K). p = [r, p, K]."""
    return np.array([p[0] * (x[0] ** p[1]) * (1.0 - x[0] / p[2])])


def SEIR_unreported(t, x, p, extra0=None):
    """SEIR model with reported (I) and unreported (U) infectious states.

    p = [beta0, beta1, alpha, q1, rho, k, gamma1, N]; states S, E, I, U, R, C.
    Transmission declines exponentially from beta0 to beta1 after t0 = 35.
    """
    beta0, beta1, alpha, q1, rho, k, gamma1, N = (p[0], p[1], p[2], p[3],
                                                   p[4], p[5], p[6], p[7])
    t0 = 35.0
    if t < t0:
        betaf = beta0
    else:
        betaf = beta1 + (beta0 - beta1) * math.exp(-q1 * (t - t0))
    S, E, I, U = x[0], x[1], x[2], x[3]
    force = betaf * S * (I + U) ** alpha / N
    return np.array([
        -force,
        force - k * E,
        k * rho * E - gamma1 * I,
        k * (1 - rho) * E - gamma1 * U,
        gamma1 * (I + U),
        k * rho * E,
    ])


# --------------------------------------------------------------------------
# Composite functions: rows of parameter draws -> one value per row
# --------------------------------------------------------------------------
def R0s(params0):
    """R0 = beta / gamma for SEIR1 (columns beta, kappa, gamma, N)."""
    P = np.atleast_2d(np.asarray(params0, dtype=float))
    return P[:, 0] / P[:, 2]


R0s_SEIR = R0s


def R0s_unreported(params0):
    """R0 = beta0 / gamma1 for SEIR_unreported."""
    P = np.atleast_2d(np.asarray(params0, dtype=float))
    return P[:, 0] / P[:, 6]
