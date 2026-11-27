"""Reference decouplers and plant-perturbation helpers."""
from __future__ import annotations
import numpy as np
from .model import PLANT_TF


def ideal_D21(tfs=PLANT_TF):
    """Ideal one-way decoupler D21 = -G21/G22 (improper by one degree)."""
    n21, d21 = map(np.asarray, tfs["G21"]); n22, d22 = map(np.asarray, tfs["G22"])
    return -np.polymul(n21, d22), np.polymul(d21, n22)


def filtered_D21(Tf: float, tfs=PLANT_TF):
    """Realizable approximation D21(s)/(Tf s + 1)."""
    num, den = ideal_D21(tfs)
    return num, np.polymul(den, [Tf, 1.0])


def static_D21(tfs=PLANT_TF):
    """Static decoupler -G21(0)/G22(0)."""
    g21 = tfs["G21"][0][-1] / tfs["G21"][1][-1]
    g22 = tfs["G22"][0][-1] / tfs["G22"][1][-1]
    return np.array([-g21 / g22]), np.array([1.0])


def scale_tf(num, den, gain=1.0, tscale=1.0):
    """Return gain * G(tscale * s): tscale > 1 makes the dynamics slower."""
    num = np.asarray(num, float); den = np.asarray(den, float)
    pn = np.arange(len(num) - 1, -1, -1); pd = np.arange(len(den) - 1, -1, -1)
    num_s = num * tscale ** pn; den_s = den * tscale ** pd
    lead = den_s[0]
    return gain * num_s / lead, den_s / lead


def perturbed_plant(gains: dict, tscales: dict, tfs=PLANT_TF):
    return {k: scale_tf(*tfs[k], gain=gains.get(k, 1.0), tscale=tscales.get(k, 1.0)) for k in tfs}
