"""Performance indices used in the paper.

All responses are deviations from the operating point (OPERATING_POINT).

* Flow disturbance rejection cost J_FD: integral of squared flow deviation
  after a unit step (+1 bar) on the pressure setpoint, over a fixed horizon:
      J_FD = sum_k  y2(k)^2 * dt     [(m3/h)^2 s]
* Interaction overshoot: maximum absolute flow deviation caused by the pressure
  setpoint step, relative to the nominal flow x2_bar, in percent.
* T95: time for the flow to first reach 95% of a unit (+1 m3/h) flow setpoint step.
"""
from __future__ import annotations
import numpy as np
from .model import LoopConfig, simulate_step, OPERATING_POINT

X2BAR = OPERATING_POINT["x2"]
FD_HORIZON = 900.0      # s, same simulated horizon as the original study
FD_DT = 0.5             # s, ZOH discretization is exact at the samples for step inputs


def disturbance_response(cfg: LoopConfig, t_end=FD_HORIZON, dt=FD_DT):
    return simulate_step(cfg, channel=0, t_end=t_end, dt=dt)


def J_FD(cfg: LoopConfig, t_end=FD_HORIZON, dt=FD_DT) -> float:
    _, Y = disturbance_response(cfg, t_end, dt)
    y2 = Y[:, 1]
    if not np.all(np.isfinite(y2)) or np.max(np.abs(y2)) > 1e3:
        return np.inf
    return float(np.sum(y2 ** 2) * dt)


def interaction_overshoot(cfg: LoopConfig, t_end=FD_HORIZON, dt=FD_DT) -> float:
    _, Y = disturbance_response(cfg, t_end, dt)
    return float(np.max(np.abs(Y[:, 1])) / X2BAR * 100.0)


def t95(cfg: LoopConfig, t_end=20000.0, dt=1.0) -> float:
    t, Y = simulate_step(cfg, channel=1, t_end=t_end, dt=dt)
    idx = np.where(Y[:, 1] >= 0.95)[0]
    return float(t[idx[0]]) if idx.size else np.nan


def summary(cfg: LoopConfig) -> dict:
    t, Yd = disturbance_response(cfg)
    ts, Ys = simulate_step(cfg, channel=1, t_end=20000.0, dt=1.0)
    idx = np.where(Ys[:, 1] >= 0.95)[0]
    return {
        "J_FD": float(np.sum(Yd[:, 1] ** 2) * FD_DT),
        "overshoot_%": float(np.max(np.abs(Yd[:, 1])) / X2BAR * 100.0),
        "T95_s": float(ts[idx[0]]) if idx.size else np.nan,
        "pressure_overshoot_%": float((np.max(Yd[:, 0]) - 1.0) * 100.0),
        "u1_range_dist": (float(Yd[:, 2].min()), float(Yd[:, 2].max())),
        "u2_range_dist": (float(Yd[:, 3].min()), float(Yd[:, 3].max())),
        "u1_range_step": (float(Ys[:, 2].min()), float(Ys[:, 2].max())),
        "u2_range_step": (float(Ys[:, 3].min()), float(Ys[:, 3].max())),
    }
