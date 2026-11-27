"""SAFEHX RSM (IAE criterion), initialized at the PIDDec equivalent to the static decoupler,
the same initialization used for the benchmark processes."""
import json, numpy as np
from safehx import rsm
from safehx.model import LoopConfig, K11_NOMINAL
from safehx.baselines import static_D21
import run_rsm_iae  # noqa: sets the IAE cost
kst = float(static_D21()[0][0])                       # -G21(0)/G22(0)
c0 = [kst * K11_NOMINAL.Kp, K11_NOMINAL.Ti, 0.0, 0.1, 1.0, 1.0]
print("initial center", c0, flush=True)
rsm.run(c0, LoopConfig(), fx0=0.4, fadd0=0.05 * K11_NOMINAL.Ti, log_path="results/rsm_iae_history_staticinit.json", max_iter=200)
