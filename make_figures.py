"""Generate every simulation figure and table of the paper from the Python model.

Requires: results/rsm_iae_history.json (run_rsm_iae.py), results/bench2_*.json (benchmarks2.py).
Outputs: figures/*.pdf and results/*.json
"""
import json, os
import numpy as np
from scipy.stats import qmc
import matplotlib.pyplot as plt

from safehx.model import (PIDCompact, LoopConfig, K11_NOMINAL, K22_NOMINAL, PLANT_TF,
                          OPERATING_POINT, simulate_step, closed_loop, max_real_pole, replace)
from safehx.metrics import X2BAR
from safehx.baselines import filtered_D21, static_D21, ideal_D21, perturbed_plant
from safehx.benchplants import iae_cross, setup as bench_setup, static_gains, ideal_dynamic, pade
from safehx.loadsim import cl_with_load
from safehx import plotstyle as ps
from scipy.linalg import expm

os.makedirs("figures", exist_ok=True); os.makedirs("results", exist_ok=True)
ps.apply()
X1BAR = OPERATING_POINT["x1"]; U1BAR = OPERATING_POINT["u1"]; U2BAR = OPERATING_POINT["u2"]
TF_FILTER = 0.1           # s
KP22_FACTOR = 103         # retuned flow controller with decoupling
KP22_ND = 7               # largest acceptable increase without decoupling
SAFEHX_IAE = {"H": 2000.0, "dt": 0.5}

hist = json.load(open("results/rsm_iae_history.json"))
best = min(hist, key=lambda h: h["J_center"])
PID_OPT = PIDCompact.from_vector(best["center"])
K22_FAST = replace(K22_NOMINAL, Kp=K22_NOMINAL.Kp * KP22_FACTOR)
K22_ND = replace(K22_NOMINAL, Kp=K22_NOMINAL.Kp * KP22_ND)

CFG = {"none": LoopConfig(), "static": LoopConfig(D21_tf=static_D21()),
       "filt": LoopConfig(D21_tf=filtered_D21(TF_FILTER)), "pid": LoopConfig(D21=PID_OPT)}


def save(fig, name):
    fig.savefig(f"figures/{name}.pdf"); fig.savefig(f"figures/{name}.png", dpi=200); plt.close(fig)


def iae_load_total(cfg, H, dt, scale=1.0):
    """Sum of IAE of both outputs for unit load steps at each plant input (separately)."""
    A, B, C, D = cl_with_load(cfg); n = A.shape[0]; tot = 0.0
    for ch in (2, 3):
        M = np.zeros((n + 1, n + 1)); M[:n, :n] = A * dt; M[:n, n] = B[:, ch] * dt
        E = expm(M); Ad, bd = E[:n, :n], E[:n, n]; x = np.zeros(n)
        for _ in range(int(H / dt)):
            y = C[:2] @ x + D[:2, ch]; tot += np.sum(np.abs(y)) * dt; x = Ad @ x + bd
    return tot


def load_response(cfg, ch, H, dt):
    A, B, C, D = cl_with_load(cfg); n = A.shape[0]
    M = np.zeros((n + 1, n + 1)); M[:n, :n] = A * dt; M[:n, n] = B[:, ch] * dt
    E = expm(M); Ad, bd = E[:n, :n], E[:n, n]; x = np.zeros(n); out = []
    for _ in range(int(H / dt) + 1):
        out.append(C @ x + D[:, ch]); x = Ad @ x + bd
    return np.arange(len(out)) * dt, np.array(out)


def gain_margin(cfg, hi=50.0):
    """Largest uniform multiplicative factor on all plant gains that keeps the loop stable."""
    def stable(k):
        pl = perturbed_plant({g: k for g in PLANT_TF}, {})
        return max_real_pole(cfg.with_(plant=pl)) < 0
    if stable(hi): return float(">%.0f" % hi) if False else np.inf
    lo, up = 1.0, hi
    for _ in range(40):
        m = 0.5 * (lo + up); lo, up = (m, up) if stable(m) else (lo, m)
    return lo


# ============================================================ 1. interaction without decoupling
t, Yp = simulate_step(CFG["none"], 0, 400, 0.1)
_, Yf = simulate_step(CFG["none"], 1, 400, 0.1)
fig, ax = plt.subplots(2, 2, figsize=(ps.TWO_COL_W, 3.3), sharex=True)
for j, (Y, title) in enumerate([(Yp, "Unit step on the pressure setpoint"), (Yf, "Unit step on the flow setpoint")]):
    a = ax[0, j]; a.plot(t, Y[:, 0], color=ps.INK, lw=1.4, label="Pressure")
    a.plot(t, (1 if j == 0 else 0) + 0 * t, **ps.REF, label="Setpoint")
    a.set_ylabel(r"$\Delta$ pressure [bar]"); a.set_title(title)
    a = ax[1, j]; a.plot(t, Y[:, 1], color=ps.C["none"]["color"], lw=1.4, label="Flow")
    a.plot(t, (1 if j == 1 else 0) + 0 * t, **ps.REF, label="Setpoint")
    a.set_ylabel(r"$\Delta$ flow [m$^3$/h]"); a.set_xlabel("Time [s]")
ax[0, 0].legend(loc="upper right"); ax[1, 1].legend(loc="lower right")
for a, s in zip(ax.flat, "abcd"): ps.panel_label(a, f"({s})")
fig.tight_layout(); save(fig, "fig_interaction_no_decoupler")

# ============================================================ 2. RSM history (SAFEHX)
it = np.array([h["iteration"] for h in hist]); Jc = np.array([h["J_center"] for h in hist])
I0 = iae_cross(CFG["none"], SAFEHX_IAE)
ov_it, t95_it = [], []
for h in hist:
    c = LoopConfig(D21=PIDCompact.from_vector(h["center"]))
    tt, Y = simulate_step(c, 0, 900, 0.5); ov_it.append(np.max(np.abs(Y[:, 1])) / X2BAR * 100)
    tt, Y = simulate_step(c, 1, 20000, 1.0); idx = np.where(Y[:, 1] >= 0.95)[0]
    t95_it.append(tt[idx[0]] if idx.size else np.nan)
fig, ax = plt.subplots(3, 1, figsize=(ps.COL_W, 4.2), sharex=True)
ax[0].semilogy(it, Jc / I0, color=ps.C["pid"]["color"], lw=0.9, marker="o", ms=2.5, label="Simulated at hypercube center")
ax[0].semilogy(it, np.minimum.accumulate(Jc) / I0, color=ps.INK, lw=1.2, ls="--", label="Best so far")
ax[0].set_ylabel(r"$IAE_{rel}$"); ax[0].legend(loc="lower left")
ax[1].semilogy(it, ov_it, color=ps.C["pid"]["color"], lw=1.2); ax[1].set_ylabel("Interaction\novershoot [%]")
ax[2].plot(it, t95_it, color=ps.C["pid"]["color"], lw=1.2); ax[2].set_ylabel(r"$T_{95}$ [s]"); ax[2].set_xlabel("RSM iteration")
for a, s in zip(ax, "abc"): ps.panel_label(a, f"({s})")
fig.tight_layout(); save(fig, "fig_rsm_history")
json.dump([{"iteration": int(i), "IAE": float(j), "IAE_rel": float(j / I0), "fx": h["fx"], "accepted": h["accepted"],
            "overshoot_%": float(o), "T95_s": float(t9), "r2": h["r2"], "alpha": h["alpha"], "loo_rmse": h["loo_rmse"],
            "J_pred_next": h["J_pred_opt"], "on_bound": h["on_bound"], "center": h["center"]}
           for i, j, o, t9, h in zip(it, Jc, ov_it, t95_it, hist)], open("results/rsm_history_table.json", "w"), indent=1)

# ============================================================ 3. SAFEHX decoupler comparison (nominal K22)
fig, ax = plt.subplots(2, 1, figsize=(ps.COL_W, 3.4), sharex=True)
for k in ["none", "static", "filt", "pid"]:
    t, Y = simulate_step(CFG[k], 0, 300, 0.05)
    ps.line(ax[0], t, Y[:, 1], k); ps.line(ax[1], t, Y[:, 3], k)
ax[0].set_ylabel(r"$\Delta$ flow [m$^3$/h]"); ax[1].set_ylabel(r"$\Delta$ VFD speed [%]"); ax[1].set_xlabel("Time [s]")
ax[0].legend(loc="lower right"); ax[0].set_title("Unit step on the pressure setpoint")
for a, s in zip(ax, "ab"): ps.panel_label(a, f"({s})")
fig.tight_layout(); save(fig, "fig_disturbance_comparison")
fig, ax = plt.subplots(1, 1, figsize=(ps.COL_W, 1.9))
for k in ["static", "filt", "pid"]:
    t, Y = simulate_step(CFG[k], 0, 300, 0.05); ps.line(ax, t, Y[:, 1], k)
ax.axhline(0, **ps.REF); ax.set_ylabel(r"$\Delta$ flow [m$^3$/h]"); ax.set_xlabel("Time [s]"); ax.legend(loc="lower right")
fig.tight_layout(); save(fig, "fig_disturbance_zoom")

# ============================================================ 4. flow setpoint and retuning (2x2)
fig, ax = plt.subplots(2, 2, figsize=(ps.TWO_COL_W, 3.8))
for key, cfg, lab in [("none", CFG["none"], r"No decoupler"), ("pid", CFG["pid"], r"PIDDec")]:
    t, Y = simulate_step(cfg, 1, 8000, 1.0); ps.line(ax[0, 0], t, Y[:, 1], key, label=lab)
ax[0, 0].axhline(0.95, **ps.REF, label="95% of step"); ax[0, 0].set_title(r"Flow setpoint step, nominal $K_{22}$")
ax[0, 0].legend(loc="lower right")
for key, cfg, lab in [("pid", CFG["pid"].with_(K22=K22_FAST), rf"PIDDec, $K_{{p,22}}\times{KP22_FACTOR}$"),
                      ("none", LoopConfig(K22=K22_ND), rf"No decoupler, $K_{{p,22}}\times{KP22_ND}$")]:
    t, Y = simulate_step(cfg, 1, 300, 0.05); ps.line(ax[0, 1], t, Y[:, 1], key, label=lab)
ax[0, 1].axhline(1, **ps.REF); ax[0, 1].set_title(r"Flow setpoint step, retuned $K_{22}$"); ax[0, 1].legend(loc="lower right")
t, Y = simulate_step(CFG["pid"], 0, 300, 0.05); ax[1, 0].plot(t, Y[:, 1], color=ps.C["pid"]["color"], ls="-", label=r"PIDDec, nominal $K_{22}$")
t, Y = simulate_step(CFG["pid"].with_(K22=K22_FAST), 0, 300, 0.05); ax[1, 0].plot(t, Y[:, 1], color=ps.INK, ls="--", label=rf"PIDDec, $K_{{p,22}}\times{KP22_FACTOR}$")
ax[1, 0].set_title("Pressure setpoint step"); ax[1, 0].legend(loc="lower right")
for key, cfg, lab in [("pid", CFG["pid"].with_(K22=K22_FAST), rf"PIDDec, $K_{{p,22}}\times{KP22_FACTOR}$"),
                      ("none", LoopConfig(K22=K22_ND), rf"No decoupler, $K_{{p,22}}\times{KP22_ND}$")]:
    t, Y = simulate_step(cfg, 0, 300, 0.05); ps.line(ax[1, 1], t, Y[:, 1], key, label=lab)
ax[1, 1].set_title(r"Pressure setpoint step, retuned $K_{22}$"); ax[1, 1].legend(loc="lower right")
for a, s in zip(ax.flat, "abcd"):
    ps.panel_label(a, f"({s})"); a.set_xlabel("Time [s]"); a.set_ylabel(r"$\Delta$ flow [m$^3$/h]")
fig.tight_layout(); save(fig, "fig_retuning")

# manipulated variables, absolute, final configuration, setpoint and load
fig, ax = plt.subplots(1, 3, figsize=(ps.TWO_COL_W, 1.9), sharey=True)
fin = CFG["pid"].with_(K22=K22_FAST)
for j, (title, fn) in enumerate([("Pressure setpoint step", lambda: simulate_step(fin, 0, 150, 0.05)),
                                 ("Flow setpoint step", lambda: simulate_step(fin, 1, 150, 0.05)),
                                 ("1% load step at the valve", lambda: load_response(fin, 2, 150, 0.05))]):
    t, Y = fn()
    ax[j].plot(t, U1BAR + Y[:, 2], color=ps.INK, lw=1.2, label="Valve opening")
    ax[j].plot(t, U2BAR + Y[:, 3], color=ps.C["pid"]["color"], lw=1.2, label="VFD speed")
    ax[j].set_title(title); ax[j].set_xlabel("Time [s]"); ax[j].set_ylim(0, 100); ps.panel_label(ax[j], f"({'abc'[j]})")
ax[0].set_ylabel("Actuator [%]"); ax[0].legend(loc="center right")
fig.tight_layout(); save(fig, "fig_manipulated_variables")

# ============================================================ 5. Bode: PIDDec vs filtered ideal and static
w = np.logspace(-4, 2, 600); s = 1j * w
fr = lambda num, den: np.polyval(num, s) / np.polyval(den, s)
K11y = fr(*K11_NOMINAL.tf("y"))
target = fr(*ideal_D21()) * K11y / (TF_FILTER * s + 1)
static_e = fr(*static_D21()) * K11y
pid_y = fr(*PID_OPT.tf("y"))
fig, ax = plt.subplots(2, 1, figsize=(ps.COL_W, 3.4), sharex=True)
for val, key, lab in [(target, "filt", r"$D_{21}^{f}(s)K_{11}(s)$"), (static_e, "static", r"$D_{21}(0)K_{11}(s)$"),
                      (pid_y, "pid", "PIDDec (RSM)")]:
    ax[0].loglog(w, np.abs(val), color=ps.C[key]["color"], ls=ps.C[key]["ls"], label=lab)
    ax[1].semilogx(w, np.unwrap(np.angle(val)) * 180 / np.pi, color=ps.C[key]["color"], ls=ps.C[key]["ls"])
ax[0].set_ylabel("Magnitude"); ax[1].set_ylabel("Phase [deg]"); ax[1].set_xlabel(r"$\omega$ [rad/s]")
ax[0].legend(loc="lower left")
for a, s_ in zip(ax, "ab"): ps.panel_label(a, f"({s_})")
fig.tight_layout(); save(fig, "fig_bode_compensator")

# ============================================================ 6. RGA
G = {k: fr(*PLANT_TF[k]) for k in PLANT_TF}
lam_plant = 1.0 / (1.0 - G["G12"] * G["G21"] / (G["G11"] * G["G22"]))
def closed_loop_fr(cfg):
    A, B, C, D = closed_loop(cfg); n = A.shape[0]
    return np.array([C[:2] @ np.linalg.solve(1j * wk * np.eye(n) - A, B) + D[:2] for wk in w])
fig, ax = plt.subplots(2, 1, figsize=(ps.COL_W, 3.4), sharex=True)
ax[0].semilogx(w, np.abs(lam_plant), color=ps.INK, lw=1.4)
ax[0].axhline(0.5, **ps.REF); ax[0].set_ylabel(r"$|\lambda_{11}|$ of $G(j\omega)$"); ax[0].set_title("Open-loop RGA (pairing)")
for k in ["none", "pid"]:
    T = closed_loop_fr(CFG[k]); lam12 = np.abs([Tk * np.linalg.inv(Tk).T for Tk in T])[:, 0, 1]
    ps.line(ax[1], w, lam12, k)
ax[1].set_ylabel(r"$|\lambda_{12}|$ of $T(j\omega)$"); ax[1].set_xlabel(r"$\omega$ [rad/s]")
ax[1].set_title("Closed-loop interaction"); ax[1].legend(loc="upper left")
for a, s_ in zip(ax, "ab"): ps.panel_label(a, f"({s_})")
fig.tight_layout(); save(fig, "fig_rga")
json.dump({"lambda11_static": float(np.real(lam_plant[0]))}, open("results/rga.json", "w"))

# ============================================================ 7. SAFEHX performance and stability table
L0_nom = iae_load_total(CFG["none"], 4000, 1.0)
rows = []
for name, cfg in [("No decoupler", CFG["none"]), ("Static decoupler", CFG["static"]),
                  ("Filtered ideal D21", CFG["filt"]), ("PIDDec (RSM)", CFG["pid"]),
                  ("No decoupler, Kp22 x7", LoopConfig(K22=K22_ND)),
                  ("Static decoupler, Kp22 x103", CFG["static"].with_(K22=K22_FAST)),
                  ("Filtered ideal D21, Kp22 x103", CFG["filt"].with_(K22=K22_FAST)),
                  ("PIDDec (RSM), Kp22 x103", CFG["pid"].with_(K22=K22_FAST))]:
    tt, Yd = simulate_step(cfg, 0, 900, 0.5); ts, Ys = simulate_step(cfg, 1, 20000, 1.0)
    idx = np.where(Ys[:, 1] >= 0.95)[0]
    ev = np.linalg.eigvals(closed_loop(cfg)[0]); ev = ev[np.abs(ev) > 1e-9]
    rows.append({"case": name, "IAE_rel_setpoint": iae_cross(cfg, SAFEHX_IAE) / I0,
                 "IAE_rel_load": iae_load_total(cfg, 4000, 1.0) / L0_nom,
                 "interaction_overshoot_%": float(np.max(np.abs(Yd[:, 1])) / X2BAR * 100),
                 "T95_s": float(ts[idx[0]]) if idx.size else None,
                 "slowest_pole": float(np.max(ev.real)), "gain_margin": gain_margin(cfg),
                 "u2_abs_range_%": [U2BAR + float(min(Yd[:, 3].min(), Ys[:, 3].min())), U2BAR + float(max(Yd[:, 3].max(), Ys[:, 3].max()))]})
    print(rows[-1], flush=True)
rows.append({"case": "No decoupler, Kp22 x103", "slowest_pole": float(max_real_pole(LoopConfig(K22=K22_FAST)))})
json.dump(rows, open("results/performance_table.json", "w"), indent=1)

# ============================================================ 8. Monte Carlo (SAFEHX, retuned K22)
N = 300
U = qmc.scale(qmc.LatinHypercube(d=8, seed=1).random(N), [0.8] * 8, [1.2] * 8)
keys = ["G11", "G12", "G21", "G22"]
mc = {k: {"iae": [], "ov": [], "stable": [], "T95": []} for k in ["none", "static", "filt", "pid"]}
for row in U:
    pl = perturbed_plant(dict(zip(keys, row[:4])), dict(zip(keys, row[4:])))
    for k in mc:
        cfg = CFG[k].with_(plant=pl, K22=K22_FAST if k != "none" else K22_ND)
        stable = max_real_pole(cfg) < 0; mc[k]["stable"].append(stable)
        if stable:
            t, Y = simulate_step(cfg, 0, 300, 0.1)
            mc[k]["ov"].append(float(np.max(np.abs(Y[:, 1])) / X2BAR * 100))
            mc[k]["iae"].append(float(np.sum(np.abs(Y[:, 1])) * 0.1))
            t, Y = simulate_step(cfg, 1, 600, 0.1); idx = np.where(Y[:, 1] >= 0.95)[0]
            mc[k]["T95"].append(float(t[idx[0]]) if idx.size else np.nan)
json.dump({k: {m: [bool(x) if m == "stable" else float(x) for x in v] for m, v in d.items()} for k, d in mc.items()},
          open("results/montecarlo.json", "w"))
order = ["none", "static", "filt", "pid"]
fig, ax = plt.subplots(1, 2, figsize=(ps.TWO_COL_W, 2.2))
for j, (m, lab) in enumerate([("ov", "Interaction overshoot [%]"), ("T95", r"$T_{95}$ [s]")]):
    data = [np.array(mc[k][m])[~np.isnan(mc[k][m])] for k in order]
    bp = ax[j].boxplot(data, patch_artist=True, widths=0.55,
                       flierprops={"marker": ".", "markersize": 2, "markeredgecolor": ps.MUTED},
                       medianprops={"color": ps.INK, "lw": 1.2})
    for patch, k in zip(bp["boxes"], order):
        patch.set_facecolor(ps.C[k]["color"]); patch.set_alpha(0.55); patch.set_edgecolor(ps.C[k]["color"])
    ax[j].set_xticks(range(1, 5)); ax[j].set_xticklabels([rf"None ($\times{KP22_ND}$)", "Static", r"$D_{21}^f$", "PIDDec"])
    ax[j].set_ylabel(lab); ps.panel_label(ax[j], f"({'ab'[j]})")
ax[0].set_yscale("log")
fig.tight_layout(); save(fig, "fig_montecarlo")

# ============================================================ 9. Benchmarks (Hylte, Wood-Berry)
B = {n: json.load(open(f"results/bench2_{n}.json")) for n in ["Hylte", "WoodBerry"]}
LOAD = json.load(open("results/bench2_load.json"))
bench_rows = []
for n in ["Hylte", "WoodBerry"]:
    b = B[n]
    bench_rows.append({"plant": n, "static_rule": b["rel_static_inverted"], "static_opt": b["rel_static_inverted_opt"],
                       "ideal": b["rel_ideal_dynamic"], "ideal_realizable": b["ideal_dynamic_realizable"],
                       "pid": b["rel_pid_rsm"], "iterations": b["iterations"],
                       "load_static": LOAD[n]["static_inv"], "load_pid": LOAD[n]["pid"], "load_ideal": LOAD[n].get("ideal")})
json.dump(bench_rows, open("results/benchmark_table.json", "w"), indent=1)

fig, ax = plt.subplots(2, 2, figsize=(ps.TWO_COL_W, 3.6))
for j, n in enumerate(["Hylte", "WoodBerry"]):
    base, P = bench_setup(n); KF1, KF2 = static_gains(n); d = 1 - KF1 * KF2
    cf = {"none": base,
          "static": base.with_(K11=replace(base.K11, Kp=base.K11.Kp / d), K22=replace(base.K22, Kp=base.K22.Kp / d),
                               D12_tf=(np.array([KF1]), np.array([1.0])), D21_tf=(np.array([KF2]), np.array([1.0]))),
          "pid": base.with_(D21=PIDCompact.from_vector(B[n]["pid"][:6]), D12=PIDCompact.from_vector(B[n]["pid"][6:]))}
    if B[n]["ideal_dynamic_realizable"]:
        dyn = {}
        for wch in ["D12", "D21"]:
            num, den, L = ideal_dynamic(n, wch); pn, pd = pade(L); dyn[wch + "_tf"] = (np.polymul(num, pn), np.polymul(den, pd))
        cf["filt"] = base.with_(**dyn)
    H = P["H"] * (0.5 if n == "Hylte" else 0.5)
    for k, c in cf.items():
        lab = {"none": "No decoupler", "static": "Static (Hägglund et al.)", "pid": "PIDDec (RSM)", "filt": "Ideal dynamic"}[k]
        t, Y = simulate_step(c, 0, H, P["dt"]); ps.line(ax[0, j], t, Y[:, 1], k, label=lab)
        t, Y = simulate_step(c, 1, H, P["dt"]); ps.line(ax[1, j], t, Y[:, 0], k, label=lab)
    unit = P["unit"]
    ax[0, j].set_title(f"{'Hylte pulp flow' if n == 'Hylte' else 'Wood–Berry column'}: step on $r_1$")
    ax[1, j].set_title(f"{'Hylte pulp flow' if n == 'Hylte' else 'Wood–Berry column'}: step on $r_2$")
    ax[0, j].set_ylabel(r"$y_2$"); ax[1, j].set_ylabel(r"$y_1$")
    for a in ax[:, j]: a.set_xlabel(f"Time [{unit}]")
h_, l_ = ax[0, 1].get_legend_handles_labels()
fig.legend(h_, l_, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.02))
for a, s_ in zip(ax.flat, "abcd"): ps.panel_label(a, f"({s_})")
fig.tight_layout(rect=(0, 0.06, 1, 1)); save(fig, "fig_benchmarks")

print("PIDDec optimized (IAE):", best["iteration"], best["J_center"] / I0, PID_OPT)
for k in mc:
    print(k, "stable", np.mean(mc[k]["stable"]), "median ov %.2f" % np.median(mc[k]["ov"]),
          "p95 ov %.2f" % np.percentile(mc[k]["ov"], 95), "median T95 %.1f" % np.nanmedian(mc[k]["T95"]))
