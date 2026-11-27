"""RSM tuning of the SAFEHX PIDDec with the IAE criterion of Hagglund et al. (2022)."""
import sys, numpy as np
from safehx import rsm
from safehx.model import LoopConfig, max_real_pole
from safehx.benchplants import iae_cross
SAFEHX_IAE = {"H": 2000.0, "dt": 0.5}
def cost(cfg):
    if max_real_pole(cfg) >= 0: return np.inf
    return iae_cross(cfg, SAFEHX_IAE)
rsm.COST["fn"] = cost
if __name__ == "__main__":
    center0 = [6.0, 1000.0, 0.0001, 1.0, 0.01, 1.0]   # initial point (Table I)
    rsm.run(center0, LoopConfig(), log_path=sys.argv[1], max_iter=int(sys.argv[2]))
