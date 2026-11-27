"""Two-way decoupling on the Hagglund et al. (2022) benchmark plants."""
import sys, json, numpy as np
from safehx.model import PIDCompact, replace, max_real_pole
from safehx.benchplants import setup, static_gains, ideal_dynamic, iae_cross, pade
from safehx import rsm

name = sys.argv[1]; maxit = int(sys.argv[2]) if len(sys.argv) > 2 else 120
base, P = setup(name); KF1, KF2 = static_gains(name)
def cost(cfg):
    if max_real_pole(cfg) >= 0: return np.inf
    return iae_cross(cfg, P)
I0 = cost(base)
out = {"IAE_none": I0}
# Hagglund static inverted decoupling (as in the paper, SISO gains)
d = 1 - KF1 * KF2
inv = base.with_(K11=replace(base.K11, Kp=base.K11.Kp / d), K22=replace(base.K22, Kp=base.K22.Kp / d),
                 D12_tf=(np.array([KF1]), np.array([1.0])), D21_tf=(np.array([KF2]), np.array([1.0])))
out["rel_static_inverted"] = cost(inv) / I0
# conventional static decoupling, same PI controllers
conv = base.with_(D12_tf=(np.array([KF1]), np.array([1.0])), D21_tf=(np.array([KF2]), np.array([1.0])))
out["rel_static_conventional"] = cost(conv) / I0
# ideal simplified decoupler (lead-lag + delay); non-causal delays are dropped
dyn = {}
for wch in ["D12", "D21"]:
    num, den, L = ideal_dynamic(name, wch)
    out[f"{wch}_ideal_delay"] = L
    if L > 0:
        pn, pd = pade(L); num, den = np.polymul(num, pn), np.polymul(den, pd)
    dyn[wch + "_tf"] = (num, den)
out["rel_ideal_dynamic"] = cost(base.with_(**dyn)) / I0
out["ideal_dynamic_realizable"] = bool(out["D12_ideal_delay"] >= 0 and out["D21_ideal_delay"] >= 0)
# two-way PIDDec by RSM, started at the static decoupler seen from the errors
K11, K22 = base.K11, base.K22
c0 = [KF2 * K11.Kp, K11.Ti, 0.0, 0.1, 1.0, 1.0,     # D21 (e1 -> u2)
      KF1 * K22.Kp, K22.Ti, 0.0, 0.1, 1.0, 1.0]     # D12 (e2 -> u1)
to_cfg = lambda th: base.with_(D21=PIDCompact.from_vector(th[:6]), D12=PIDCompact.from_vector(th[6:]))
out["rel_pid_init"] = cost(to_cfg(np.array(c0))) / I0
hist = rsm.run_general(c0, to_cfg, cost, td_idx=[2, 8], fadd0=[0.05 * K11.Ti, 0.05 * K22.Ti],
                       max_iter=maxit, log_path=f"results/bench2_{name}_hist.json")
b = min(hist, key=lambda h: h["J_center"])
out.update({"rel_pid_rsm": b["J_center"] / I0, "pid": b["center"], "iterations": len(hist)})
json.dump(out, open(f"results/bench2_{name}.json", "w"), indent=1)
print(json.dumps(out, indent=1))
