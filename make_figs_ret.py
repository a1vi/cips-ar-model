import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cips_ar import Params, AnalyticWindow
from retention import RetentionWindow
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.linewidth": 0.8, "axes.labelsize": 9,
    "xtick.direction": "in", "ytick.direction": "in", "legend.frameon": False, "legend.fontsize": 7.5,
    "savefig.dpi": 300, "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.5})
OI = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
D = json.load(open("data/retention.json"))
def save(f, n):
    f.savefig("figures/%s.png" % n, bbox_inches="tight"); f.savefig("figures/%s.pdf" % n, bbox_inches="tight"); plt.close(f)
fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.7))
ax = axs[0]
for i, a in enumerate((0.1, 0.25, 0.5, 1.0, 2.0)):
    ax.semilogx(D["t"], D["time"][str(a)], color=OI[i], label="A$_R$ = %g" % a)
ax.set_xlabel("hold time (s)"); ax.set_ylabel("memory window (V)"); ax.legend(loc="upper right", fontsize=6.5)
ax.axvline(3.15e8, color="0.6", lw=0.8, ls=":"); ax.text(0.03, 0.04, "a", transform=ax.transAxes, fontweight="bold", va="bottom", ha="left") if False else ax.text(-0.18, 1.02, "a", transform=ax.transAxes, fontweight="bold")
AR = np.array(D["AR"])
ax = axs[1]
ax.semilogx(AR, [AnalyticWindow(Params(A_R=a, xp=6.0)).MW() for a in AR], "k--", lw=1, label="quasi-static")
for i, b in enumerate((10, 20, 40)):
    mw = [RetentionWindow(Params(A_R=a, xp=6.0), b=b).MW_t(3.15e8) for a in AR]
    ax.semilogx(AR, mw, color=OI[i], label="10 yr, b = %d" % b)
ax.set_xlabel("A$_R$"); ax.set_ylabel("memory window (V)"); ax.legend(fontsize=6.5, loc="lower left"); ax.text(-0.18, 1.02, "b", transform=ax.transAxes, fontweight="bold")
ax = axs[2]
for i, (xp, ls) in enumerate(((6.0, "-"), (10.0, "--"))):
    for j, b in enumerate((10, 40)):
        fr = [RetentionWindow(Params(A_R=a, xp=xp), b=b).MW_t(3.15e8) / max(AnalyticWindow(Params(A_R=a, xp=xp)).MW(), 1e-9) for a in AR]
        ax.semilogx(AR, fr, ls, color=OI[j * 2], label="x$_p$ = %g V, b = %d" % (xp, b))
ax.set_ylim(0, 1); ax.set_xlabel("A$_R$"); ax.set_ylabel("fraction of window left (10 yr)"); ax.legend(fontsize=6.3, loc="lower left"); ax.text(-0.18, 1.02, "c", transform=ax.transAxes, fontweight="bold")
fig.tight_layout(w_pad=1.0); save(fig, "fig14_retention")
# 3D surface
S = D["surf"]; X, Y = np.meshgrid(np.log10(S["AR"]), np.log10(S["t"])); Z = np.array(S["MW"])
f = plt.figure(figsize=(5.0, 3.5)); ax = f.add_subplot(111, projection="3d")
s = ax.plot_surface(X, Y, Z, cmap="viridis", edgecolor="none", rstride=1, cstride=1, vmin=0, vmax=1.75)
ax.set_xlabel("log$_{10}$ A$_R$", labelpad=1); ax.set_ylabel("log$_{10}$ t (s)", labelpad=1); ax.set_zlim(0, 1.8)
ax.view_init(26, -58); ax.tick_params(labelsize=7, pad=-2); ax.set_title("window (V), x$_p$ = 6 V, b = 20", fontsize=8.5, y=0.97); ax.text2D(0.02, 0.95, "d", transform=ax.transAxes, fontweight="bold", fontsize=10)
f.colorbar(s, ax=ax, shrink=0.55, pad=0.12, label="window (V)"); save(f, "fig15_3d_retention")
