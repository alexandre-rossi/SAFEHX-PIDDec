"""Benchmark TITO plants from Hagglund et al. (2022), with their PI tunings."""
from math import factorial
import numpy as np
from .model import PIDCompact, LoopConfig, simulate_step


def pade(T, n=6):
    if T == 0:
        return np.array([1.0]), np.array([1.0])
    c = [factorial(2 * n - k) * factorial(n) / (factorial(2 * n) * factorial(k) * factorial(n - k)) for k in range(n + 1)]
    num = np.array([c[k] * (-T) ** k for k in range(n + 1)])[::-1]
    den = np.array([c[k] * T ** k for k in range(n + 1)])[::-1]
    return num, den


def foptd(K, T, L):
    pn, pd = pade(L)
    return np.polymul([K], pn), np.polymul([T, 1.0], pd)


# (K, T, L) of each FOPTD element; time unit noted in "unit"
FOPTD = {
    "Hylte":     {"G11": (1.8, 0.42, 2.0), "G12": (-2.0, 0.90, 1.4), "G21": (0.41, 3.0, 2.64), "G22": (1.2, 3.0, 4.0),
                  "K1": (0.185, 1.087), "K2": (0.278, 4.333), "H": 60.0, "dt": 0.01, "unit": "s"},
    "WoodBerry": {"G11": (12.8, 16.7, 1.0), "G12": (-18.9, 21.0, 3.0), "G21": (6.6, 10.9, 7.0), "G22": (-19.4, 14.4, 3.0),
                  "K1": (0.0732, 16.7), "K2": (-0.0427, 14.4), "H": 300.0, "dt": 0.1, "unit": "min"},
}


def setup(name):
    P = FOPTD[name]
    tfs = {k: foptd(*P[k]) for k in ["G11", "G12", "G21", "G22"]}
    K11 = PIDCompact(P["K1"][0], P["K1"][1], 0.0, 1.0, 1.0, 1.0)
    K22 = PIDCompact(P["K2"][0], P["K2"][1], 0.0, 1.0, 1.0, 1.0)
    return LoopConfig(K11=K11, K22=K22, plant=tfs), P


def static_gains(name):
    P = FOPTD[name]
    KF1 = -P["G12"][0] / P["G11"][0]      # D12: v2 -> u1
    KF2 = -P["G21"][0] / P["G22"][0]      # D21: v1 -> u2
    return KF1, KF2


def ideal_dynamic(name, which):
    """Ideal simplified decoupler -Gij/Gii as lead-lag + delay. Returns (num, den, delay)."""
    P = FOPTD[name]
    if which == "D21":
        (Kn, Tn, Ln), (Kd, Td, Ld) = P["G21"], P["G22"]
    else:
        (Kn, Tn, Ln), (Kd, Td, Ld) = P["G12"], P["G11"]
    k = -Kn / Kd; L = Ln - Ld
    num, den = np.array([k * Td, k]), np.array([Tn, 1.0])
    return num, den, L


def iae_cross(cfg, P):
    """Hagglund's IAE: off-diagonal responses to setpoint steps in the opposite loop."""
    H, dt = P["H"], P["dt"]
    _, Y1 = simulate_step(cfg, 0, H, dt)   # step r1 -> error in loop 2
    _, Y2 = simulate_step(cfg, 1, H, dt)   # step r2 -> error in loop 1
    return float(np.sum(np.abs(Y1[:, 1])) * dt + np.sum(np.abs(Y2[:, 0])) * dt)
