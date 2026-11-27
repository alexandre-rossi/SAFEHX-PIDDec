"""Identified SAFEHX hot-water circuit model and PID_Compact-based control loop.

Conventions
-----------
Outputs : y1 = pressure deviation [bar], y2 = flow deviation [m3/h]
Inputs  : u1 = valve opening deviation [%], u2 = VFD speed deviation [%]
G_ij    : transfer function from input j to output i (standard MIMO convention).

Every controller follows the 2-DoF PID_Compact law (Siemens), Eq. (1) of the paper:

    u = Kp [ (b w - x) + 1/(Ti s) (w - x) + Td s / (a Td s + 1) (c w - x) ]

with w the setpoint and x the measured value. The decentralized controllers are
K11 (pressure -> valve) and K22 (flow -> VFD). The cross-coupling compensator
PIDDec (D21) uses the pressure setpoint and measurement (w = r1, x = y1) and its
output is added to u2, so the overall controller is lower triangular.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np
from scipy import signal
from scipy.linalg import expm

# ---------------------------------------------------------------------------
# Identified plant (coefficients as identified in Octave)
# ---------------------------------------------------------------------------
PLANT_TF = {
    # valve -> pressure (5th order, subspace identification)
    "G11": ([-0.0002162, -0.008813, 0.0540003, -1.09254, 4.18853, -14.0799],
            [1, 8.32633, 103.285, 459.092, 1506.89, 152.834]),
    # VFD -> pressure (2nd order, OE)
    "G12": ([-0.0265788, 0.408094],
            [1, 18.9276, 4.01789]),
    # valve -> flow (5th order, subspace identification)
    "G21": ([0.0002468, 0.0197366, -0.106446, 2.73192, -8.62702, 35.1111],
            [1, 8.29138, 100.878, 454.134, 1133.01, 186.79]),
    # VFD -> flow (2nd order, OE)
    "G22": ([0.0002523, 0.0285995],
            [1, 11.9588, 1.75382]),
}

OPERATING_POINT = {"x1": 3.208, "x2": 0.7617, "u1": 25.0, "u2": 57.97}


def _tf2ss(num, den):
    A, B, C, D = signal.tf2ss(num, den)
    return A, B, C, D


def plant_ss(tfs: dict = PLANT_TF):
    """2x2 plant as one state-space model (inputs u1,u2; outputs y1,y2)."""
    blocks = {k: _tf2ss(*v) for k, v in tfs.items()}
    order = [("G11", 0, 0), ("G12", 0, 1), ("G21", 1, 0), ("G22", 1, 1)]
    n = sum(blocks[k][0].shape[0] for k, _, _ in order)
    A = np.zeros((n, n)); B = np.zeros((n, 2)); C = np.zeros((2, n)); D = np.zeros((2, 2))
    i = 0
    for k, out, inp in order:
        a, b, c, d = blocks[k]
        m = a.shape[0]
        A[i:i + m, i:i + m] = a
        B[i:i + m, inp] = b[:, 0]
        C[out, i:i + m] = c[0]
        D[out, inp] += d[0, 0]
        i += m
    return A, B, C, D


# ---------------------------------------------------------------------------
# Controllers
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class PIDCompact:
    Kp: float
    Ti: float
    Td: float = 0.0
    a: float = 0.5
    b: float = 1.0
    c: float = 0.5

    def as_vector(self):
        return np.array([self.Kp, self.Ti, self.Td, self.a, self.b, self.c])

    @staticmethod
    def from_vector(v):
        return PIDCompact(*[float(x) for x in v])

    def ss(self):
        """State space with inputs [w, x] and output u.

        States: integral of (w - x), derivative filter state z.
        """
        Kp, Ti, Td, a, b, c = self.Kp, self.Ti, self.Td, self.a, self.b, self.c
        if Td > 0 and a > 0:
            tau = a * Td
            # z' = ((c w - x) - z)/tau ; D-part = (1/a) * ((c w - x) - z)
            A = np.array([[0.0, 0.0], [0.0, -1.0 / tau]])
            B = np.array([[1.0, -1.0], [c / tau, -1.0 / tau]])
            C = Kp * np.array([[1.0 / Ti, -1.0 / a]])
            D = Kp * np.array([[b + c / a, -1.0 - 1.0 / a]])
        else:
            A = np.array([[0.0]])
            B = np.array([[1.0, -1.0]])
            C = Kp * np.array([[1.0 / Ti]])
            D = Kp * np.array([[b, -1.0]])
        return A, B, C, D

    def tf(self, which="y"):
        """Transfer function from w ('r') or from x ('y', sign included as +Cy)."""
        s = signal.lti
        Kp, Ti, Td, a, b, c = self.Kp, self.Ti, self.Td, self.a, self.b, self.c
        # Build numerator/denominator of Kp*(k0 + 1/(Ti s) + kd*Td s/(a Td s + 1))
        k0, kd = (b, c) if which == "r" else (1.0, 1.0)
        if Td > 0:
            den = np.polymul([Ti, 0], [a * Td, 1])
            num = (np.polyadd(np.polyadd(k0 * np.polymul([Ti, 0], [a * Td, 1]), [a * Td, 1]),
                              kd * np.polymul([Ti, 0], [Td, 0])))
        else:
            den = np.array([Ti, 0.0])
            num = np.array([k0 * Ti, 1.0])
        return Kp * np.atleast_1d(num), np.atleast_1d(den)


K11_NOMINAL = PIDCompact(Kp=-0.17, Ti=0.4, Td=0.8, a=0.5, b=1.0, c=0.5)
K22_NOMINAL = PIDCompact(Kp=0.0034, Ti=0.1, Td=0.035, a=0.5, b=1.0, c=0.5)


@dataclass(frozen=True)
class LoopConfig:
    K11: PIDCompact = K11_NOMINAL
    K22: PIDCompact = K22_NOMINAL
    D21: PIDCompact | None = None          # PIDDec (pressure loop -> VFD input)
    D21_tf: tuple | None = None            # alternative: classical decoupler (num, den) acting on v1 (K11 output)
    D12: PIDCompact | None = None          # PIDDec in the opposite direction (flow loop -> u1)
    D12_tf: tuple | None = None            # classical decoupler acting on v2 (K22 output), added to u1
    plant: dict | None = None              # alternative plant coefficients (sensitivity / Monte Carlo)

    def with_(self, **kw):
        return replace(self, **kw)


def _controller_ss(cfg: LoopConfig):
    """Controller with inputs [r1, r2, y1, y2] and outputs [u1, u2]."""
    blocks = []  # (A,B,C,D, input_map (w_idx, x_idx), output idx)
    blocks.append(cfg.K11.ss() + ((0, 2), 0))
    blocks.append(cfg.K22.ss() + ((1, 3), 1))
    if cfg.D21 is not None:
        blocks.append(cfg.D21.ss() + ((0, 2), 1))
    if cfg.D21_tf is not None:
        # classical decoupler acting on the pressure controller output u1:
        # u2 += D21(s) u1, implemented as the series K11 -> D21 (inputs r1, y1)
        Ad_, Bd_, Cd_, Dd_ = _tf2ss(*cfg.D21_tf)
        Ak, Bk, Ck, Dk = cfg.K11.ss()
        nk, nd = Ak.shape[0], Ad_.shape[0]
        A = np.block([[Ak, np.zeros((nk, nd))], [Bd_ @ Ck, Ad_]])
        B = np.vstack([Bk, Bd_ @ Dk])
        C = np.hstack([Dd_ @ Ck, Cd_])
        D = Dd_ @ Dk
        blocks.append((A, B, C, D, (0, 2), 1))
    if cfg.D12 is not None:
        blocks.append(cfg.D12.ss() + ((1, 3), 0))
    if cfg.D12_tf is not None:
        Ad_, Bd_, Cd_, Dd_ = _tf2ss(*cfg.D12_tf)
        Ak, Bk, Ck, Dk = cfg.K22.ss()
        nk, nd = Ak.shape[0], Ad_.shape[0]
        A = np.block([[Ak, np.zeros((nk, nd))], [Bd_ @ Ck, Ad_]])
        B = np.vstack([Bk, Bd_ @ Dk])
        C = np.hstack([Dd_ @ Ck, Cd_])
        D = Dd_ @ Dk
        blocks.append((A, B, C, D, (1, 3), 0))
    n = sum(blk[0].shape[0] for blk in blocks)
    A = np.zeros((n, n)); B = np.zeros((n, 4)); C = np.zeros((2, n)); D = np.zeros((2, 4))
    i = 0
    for a, b, c, d, (iw, ix), out in blocks:
        m = a.shape[0]
        A[i:i + m, i:i + m] = a
        B[i:i + m, iw] += b[:, 0]; B[i:i + m, ix] += b[:, 1]
        C[out, i:i + m] += c[0]
        D[out, iw] += d[0, 0]; D[out, ix] += d[0, 1]
        i += m
    return A, B, C, D


def closed_loop(cfg: LoopConfig):
    """Closed loop with inputs [r1, r2] and outputs [y1, y2, u1, u2]."""
    Ag, Bg, Cg, Dg = plant_ss(cfg.plant or PLANT_TF)
    Ak, Bk, Ck, Dk = _controller_ss(cfg)
    Bkr, Bky = Bk[:, :2], Bk[:, 2:]
    Dkr, Dky = Dk[:, :2], Dk[:, 2:]
    ng, nk = Ag.shape[0], Ak.shape[0]
    M = np.linalg.inv(np.eye(2) - Dg @ Dky)
    Cy = M @ np.hstack([Cg, Dg @ Ck]); Dyr = M @ Dg @ Dkr
    Cu = np.hstack([np.zeros((2, ng)), Ck]) + Dky @ Cy; Dur = Dkr + Dky @ Dyr
    A = np.block([[Ag, np.zeros((ng, nk))], [np.zeros((nk, ng)), Ak]]) + np.vstack([Bg @ Cu, Bky @ Cy])
    B = np.vstack([Bg @ Dur, Bkr + Bky @ Dyr])
    C = np.vstack([Cy, Cu]); D = np.vstack([Dyr, Dur])
    return A, B, C, D


def max_real_pole(cfg: LoopConfig) -> float:
    A = closed_loop(cfg)[0]
    ev = np.linalg.eigvals(A)
    # integrator modes duplicated in the 2-DoF realization appear at 0 (uncontrollable/unobservable)
    ev = ev[np.abs(ev) > 1e-9]
    return float(np.max(ev.real))


def simulate_step(cfg: LoopConfig, channel: int, t_end: float, dt: float, amplitude: float = 1.0):
    """Exact (ZOH) response to a setpoint step on r1 (channel=0) or r2 (channel=1).

    Returns time vector and array [y1, y2, u1, u2] (deviation variables).
    """
    A, B, C, D = closed_loop(cfg)
    n = A.shape[0]
    M = np.zeros((n + 1, n + 1)); M[:n, :n] = A * dt; M[:n, n] = B[:, channel] * dt
    E = expm(M); Ad, bd = E[:n, :n], E[:n, n]
    N = int(round(t_end / dt)) + 1
    X = np.empty((N, n)); x = np.zeros(n)
    for k in range(N):
        X[k] = x
        x = Ad @ x + bd * amplitude
    Y = X @ C.T + amplitude * D[:, channel]
    t = np.arange(N) * dt
    return t, Y
