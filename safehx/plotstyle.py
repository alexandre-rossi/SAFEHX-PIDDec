"""Single figure style for every plot in the paper (Elsevier 5p two-column layout)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COL_W = 3.46     # in, single column (88 mm)
TWO_COL_W = 7.09 # in, double column (180 mm)

# Categorical order (validated for CVD separation); line styles give a second channel
C = {
    "none":     {"color": "#eb6834", "ls": "--", "label": "No decoupler"},
    "pid":      {"color": "#2a78d6", "ls": "-",  "label": "PIDDec (RSM)"},
    "filt":     {"color": "#1baf7a", "ls": "-.", "label": r"$D_{21}^{f}$ (filtered ideal)"},
    "static":   {"color": "#eda100", "ls": ":",  "label": "Static decoupler"},
}
REF = {"color": "#3d3d3a", "ls": (0, (1, 1.5)), "lw": 0.9}
INK = "#3d3d3a"; MUTED = "#8a8980"; GRID = "#e4e3dc"


def apply():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 1.4, "legend.frameon": False,
        "figure.dpi": 150, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    })


def line(ax, x, y, key, **kw):
    st = C[key]
    return ax.plot(x, y, color=st["color"], ls=st["ls"], label=kw.pop("label", st["label"]), **kw)


def panel_label(ax, s):
    ax.text(0.0, 1.02, s, transform=ax.transAxes, ha="left", va="bottom", fontsize=8, fontweight="bold", color=INK)
