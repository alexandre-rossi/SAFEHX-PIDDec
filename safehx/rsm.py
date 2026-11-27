"""Iterative Response Surface Methodology (RSM) tuning of the PIDDec compensator.

Each iteration:
 1. builds a three-level full factorial design (3^6 = 729 runs) around the current
    center, in coded units (-1, 0, +1). Kp, Ti, a, b, c vary multiplicatively,
    theta*(1 +/- fx); Td varies additively, Td +/- f_add (clipped at Td >= 0);
 2. simulates every run and computes the cost J_FD;
 3. fits a full quadratic surrogate in coded units by Ridge regression, with the
    regularization chosen by efficient leave-one-out cross-validation;
 4. minimizes the surrogate inside the hypercube (L-BFGS-B, started at the best run);
 5. the minimizer becomes the center of the next hypercube, so the next center run
    is also the confirmation simulation of the predicted optimum.
"""
from __future__ import annotations
import itertools, json, time
import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import RidgeCV
from scipy.optimize import minimize
from .model import PIDCompact, LoopConfig
from .metrics import J_FD

PARAM_NAMES = ["Kp", "Ti", "Td", "a", "b", "c"]
TD_INDEX = 2
LEVELS = np.array(list(itertools.product([-1.0, 0.0, 1.0], repeat=6)))   # 729 x 6
ALPHAS = np.logspace(-6, 2, 50)


def decode(center, coded, fx, fadd):
    theta = center * (1.0 + fx * coded)
    theta[..., TD_INDEX] = np.maximum(center[TD_INDEX] + fadd * coded[..., TD_INDEX], 0.0)
    return theta


COST = {"fn": J_FD}   # replaceable cost function (cfg -> float)


def evaluate(theta, base_cfg: LoopConfig):
    return COST["fn"](base_cfg.with_(D21=PIDCompact.from_vector(theta)))


def iterate(center, fx, fadd, base_cfg: LoopConfig):
    thetas = decode(center, LEVELS, fx, fadd)
    J = np.array([evaluate(th, base_cfg) for th in thetas])
    ok = np.isfinite(J)
    poly = PolynomialFeatures(degree=2, include_bias=True)
    Xp = poly.fit_transform(LEVELS[ok])
    model = RidgeCV(alphas=ALPHAS, fit_intercept=False, store_cv_results=True)
    model.fit(Xp, J[ok])
    pred = model.predict(Xp)
    ss_res = np.sum((J[ok] - pred) ** 2); ss_tot = np.sum((J[ok] - J[ok].mean()) ** 2)
    r2 = 1 - ss_res / ss_tot
    loo_rmse = float(np.sqrt(np.min(np.mean(model.cv_results_, axis=0))))
    f = lambda z: float(model.predict(poly.transform(z.reshape(1, -1)))[0])
    z0 = LEVELS[ok][np.argmin(J[ok])]
    res = minimize(f, z0, method="L-BFGS-B", bounds=[(-1, 1)] * 6)
    z_opt = res.x
    on_bound = [PARAM_NAMES[i] for i in range(6) if abs(abs(z_opt[i]) - 1) < 1e-6]
    new_center = decode(center, z_opt, fx, fadd)
    return {
        "center": center.tolist(), "J_center": float(J[np.all(LEVELS == 0, axis=1)][0]),
        "J_best_run": float(np.min(J[ok])), "n_unstable": int((~ok).sum()),
        "r2": float(r2), "loo_rmse": loo_rmse, "alpha": float(model.alpha_),
        "J_pred_opt": float(res.fun), "z_opt": z_opt.tolist(), "on_bound": on_bound,
        "new_center": new_center.tolist(), "fx": fx, "fadd": fadd,
    }


def run(center0, base_cfg: LoopConfig, fx0=0.4, fadd0=0.0005, shrink=0.5, fx_min=0.4/64,
        rel_tol=1e-3, max_iter=150, log_path=None, verbose=True, resume=True, time_budget=None):
    """Sequential RSM with an adaptive region (trust-region-like).

    The hypercube moves to the surrogate minimizer while the simulated cost at the new
    center improves on the best value found. When it does not, the search returns to the
    best center and the region is halved (fx and f_add). After two consecutive improvements
    with the minimizer on the region boundary, the region is doubled again (up to fx0).
    Stops when the relative improvement of an accepted step is below rel_tol, or when
    the region would become smaller than fx_min. With resume=True and an existing log,
    the run continues from the stored state.
    """
    import os
    state = None
    if resume and log_path and os.path.exists(log_path):
        history = json.load(open(log_path))
        if history and "state" in history[-1]:
            state = history[-1]["state"]
    else:
        history = []
    if state:
        center = np.asarray(state["center"]); fx = state["fx"]; fadd = state["fadd"]
        best_center = np.asarray(state["best_center"]); best_J = state["best_J"]
        revisit = state["revisit"]; streak = state["streak"]
        if state.get("done"):
            return history
    else:
        center = np.asarray(center0, float); fx, fadd = fx0, fadd0
        best_center, best_J = center.copy(), np.inf; revisit, streak = False, 0
    t_start = time.time()
    while len(history) < max_iter:
        t0 = time.time()
        it = iterate(center, fx, fadd, base_cfg)
        it["iteration"] = len(history) + 1
        Jc = it["J_center"]
        improved = (Jc < best_J) and not revisit
        rel = (best_J - Jc) / best_J if np.isfinite(best_J) and best_J > 0 else np.inf
        it["accepted"] = bool(improved or revisit); it["revisit"] = revisit
        done = False
        if revisit:
            revisit = False; center = np.asarray(it["new_center"])
        elif improved:
            best_J, best_center = Jc, center.copy()
            streak = streak + 1 if it["on_bound"] else 0
            if rel < rel_tol:
                done = True
            else:
                if streak >= 2 and fx < fx0:
                    fx = min(fx / shrink, fx0); fadd = min(fadd / shrink, fadd0); streak = 0
                center = np.asarray(it["new_center"])
        else:
            streak = 0
            if fx * shrink < fx_min:
                done = True
            else:
                fx *= shrink; fadd *= shrink
                center, revisit = best_center.copy(), True
        it["state"] = {"center": center.tolist(), "fx": fx, "fadd": fadd, "best_center": best_center.tolist(),
                       "best_J": best_J, "revisit": revisit, "streak": streak, "done": done}
        history.append(it)
        if verbose:
            print(f"it {it['iteration']:3d} fx={it['fx']:.4f} J_center={Jc:.6g} best={best_J:.6g} "
                  f"R2={it['r2']:.4f} acc={it['accepted']} bound={it['on_bound']} ({time.time()-t0:.1f}s)", flush=True)
        if log_path:
            json.dump(history, open(log_path, "w"), indent=1)
        if done or (time_budget and time.time() - t_start > time_budget):
            break
    return history


# ---------------------------------------------------------------------------
# General version: any number of parameters and a Box-Behnken-type design
# ---------------------------------------------------------------------------
def bb_design(k):
    """All pairs of factors at +/-1 with the others at 0, plus the center point."""
    rows = [np.zeros(k)]
    for i, j in itertools.combinations(range(k), 2):
        for si, sj in itertools.product([-1.0, 1.0], repeat=2):
            r = np.zeros(k); r[i], r[j] = si, sj; rows.append(r)
    return np.array(rows)


def run_general(center0, to_cfg, cost, td_idx, fx0=0.4, fadd0=None, shrink=0.5, fx_min=0.4 / 64,
                rel_tol=1e-3, max_iter=150, design="bb", verbose=True, log_path=None):
    center = np.asarray(center0, float); k = center.size
    levels = bb_design(k) if design == "bb" else np.array(list(itertools.product([-1.0, 0.0, 1.0], repeat=k)))
    fadd0 = np.asarray(fadd0 if fadd0 is not None else np.full(len(td_idx), 0.0005), float)
    fx, fadd = fx0, fadd0.copy()

    def decode_g(c, z, fx, fadd):
        th = c * (1.0 + fx * z)
        for m, i in enumerate(td_idx):
            th[..., i] = np.maximum(c[i] + fadd[m] * z[..., i], 0.0)
        return th

    best_center, best_J = center.copy(), np.inf
    history, revisit, shrunk = [], False, False
    poly = PolynomialFeatures(degree=2, include_bias=True)
    for it_n in range(1, max_iter + 1):
        t0 = time.time()
        thetas = decode_g(center, levels, fx, fadd)
        J = np.array([cost(to_cfg(th)) for th in thetas]); ok = np.isfinite(J)
        Xp = poly.fit_transform(levels[ok])
        model = RidgeCV(alphas=ALPHAS, fit_intercept=False, store_cv_results=True).fit(Xp, J[ok])
        pred = model.predict(Xp)
        r2 = 1 - np.sum((J[ok] - pred) ** 2) / np.sum((J[ok] - J[ok].mean()) ** 2)
        f = lambda z: float(model.predict(poly.transform(z.reshape(1, -1)))[0])
        res = minimize(f, levels[ok][np.argmin(J[ok])], method="L-BFGS-B", bounds=[(-1, 1)] * k)
        on_bound = [int(i) for i in range(k) if abs(abs(res.x[i]) - 1) < 1e-6]
        Jc = float(J[0])
        improved = (Jc < best_J) and not revisit
        rel = (best_J - Jc) / best_J if np.isfinite(best_J) and best_J > 0 else np.inf
        history.append({"iteration": it_n, "center": center.tolist(), "J_center": Jc, "J_best_run": float(np.min(J[ok])),
                        "n_unstable": int((~ok).sum()), "r2": float(r2), "alpha": float(model.alpha_),
                        "J_pred_opt": float(res.fun), "on_bound": on_bound, "fx": fx, "accepted": bool(improved or revisit),
                        "n_runs": int(len(levels))})
        if verbose:
            print(f"it {it_n:3d} fx={fx:.4f} J={Jc:.6g} best={min(best_J, Jc):.6g} R2={r2:.3f} "
                  f"unstable={int((~ok).sum())} ({time.time()-t0:.1f}s)", flush=True)
        if log_path:
            json.dump(history, open(log_path, "w"), indent=1)
        new_center = decode_g(center, res.x, fx, fadd)
        if revisit:
            revisit = False; center = new_center; continue
        if improved:
            best_J, best_center = Jc, center.copy()
            if rel < rel_tol:
                break
            if on_bound and fx < fx0 and not shrunk:
                fx = min(fx / shrink, fx0); fadd = np.minimum(fadd / shrink, fadd0)
            center = new_center
        else:
            if fx * shrink < fx_min:
                break
            fx *= shrink; fadd = fadd * shrink
            center, revisit, shrunk = best_center.copy(), True, True
    return history
