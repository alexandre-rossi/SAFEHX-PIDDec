import numpy as np, json
from scipy.optimize import minimize
from safehx.model import replace, max_real_pole
from safehx.benchplants import setup, static_gains, iae_cross
for name, paper in [("Hylte", 0.0652), ("WoodBerry", 0.1512)]:
    base, P = setup(name); I0 = iae_cross(base, P); KF1, KF2 = static_gains(name)
    def inv(k):
        d = 1 - k[0]*k[1]
        c = base.with_(K11=replace(base.K11, Kp=base.K11.Kp/d), K22=replace(base.K22, Kp=base.K22.Kp/d),
                       D12_tf=(np.array([k[0]]), np.array([1.0])), D21_tf=(np.array([k[1]]), np.array([1.0])))
        return iae_cross(c, P)/I0 if max_real_pole(c) < 0 else 10.0
    r = minimize(inv, [KF1, KF2], method="Nelder-Mead", options={"xatol":1e-4,"fatol":1e-6})
    print(name, "optimal static inverted: KF=", np.round(r.x,3), "IAE_rel=%.4f (paper %.4f)" % (r.fun, paper))
    d = json.load(open(f"results/bench2_{name}.json")); d["rel_static_inverted_opt"] = float(r.fun); d["KF_opt"] = r.x.tolist()
    json.dump(d, open(f"results/bench2_{name}.json","w"), indent=1)
