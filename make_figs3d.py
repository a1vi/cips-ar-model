import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from cips_ar import Params, AnalyticWindow
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "savefig.dpi": 300})
R = json.load(open("data/ar_results.json"))
def save(fig, n):
    fig.savefig("figures/%s.png" % n, bbox_inches="tight"); fig.savefig("figures/%s.pdf" % n, bbox_inches="tight"); plt.close(fig)

FACES = []
def box(ax, x0, x1, y0, y1, z0, z1, color, alpha=1.0):
    v = [(x0,y0,z0),(x1,y0,z0),(x1,y1,z0),(x0,y1,z0),(x0,y0,z1),(x1,y0,z1),(x1,y1,z1),(x0,y1,z1)]
    f = [[0,1,2,3],[4,5,6,7],[0,1,5,4],[2,3,7,6],[1,2,6,5],[0,3,7,4]]
    import matplotlib.colors as mc
    rgba = mc.to_rgba(color, alpha)
    for q in f: FACES.append(([v[k] for k in q], rgba))

def flush(ax):
    global FACES
    pc = Poly3DCollection([f for f, _ in FACES], facecolors=[c for _, c in FACES], edgecolors="k", linewidths=0.5)
    ax.add_collection3d(pc); FACES = []

def stack(ax, AR):
    s = np.sqrt(AR); c = 0.5; h = 0.5 * s
    box(ax, -1.0, -0.5, -0.4, 1.4, 0, 0.14, "#dddddd")          # source
    box(ax, 1.5, 2.0, -0.4, 1.4, 0, 0.14, "#dddddd")            # drain
    box(ax, -0.5, 1.5, -0.4, 1.4, 0, 0.14, "#7fb2e5")           # channel
    box(ax, 0, 1, 0, 1, 0.14, 0.34, "#a9d8b8")                  # h-BN
    box(ax, 0, 1, 0, 1, 0.34, 0.42, "#e7c446")                  # floating gate
    box(ax, c - h, c + h, 0.5 - h, 0.5 + h, 0.42, 0.90, "#f4b4c4", 0.95)  # FE
    box(ax, c - h, c + h, 0.5 - h, 0.5 + h, 0.90, 0.98, "#bbbbbb")        # top gate
    flush(ax)
    ax.set_xlim(-1.0, 2.0); ax.set_ylim(-0.6, 1.6); ax.set_zlim(0, 1.1)
    ax.set_box_aspect((3.0, 2.2, 1.1)); ax.view_init(32, -52); ax.set_axis_off()
    ax.set_title("A$_R$ = %g" % AR, fontsize=9, y=0.9)

def fig11():
    fig = plt.figure(figsize=(7.2, 2.6))
    for i, a in enumerate([0.25, 1.0, 2.0]):
        ax = fig.add_subplot(1, 3, i + 1, projection="3d"); stack(ax, a)
    fig.legend(handles=[plt.Rectangle((0,0),1,1,fc=c,ec="k") for c in ["#bbbbbb","#f4b4c4","#e7c446","#a9d8b8","#7fb2e5"]],
               labels=["top gate","CuInP$_2$S$_6$","floating gate","h-BN","channel"], loc="lower center", ncol=5, frameon=False, fontsize=7.5)
    fig.subplots_adjust(wspace=-0.05, bottom=0.12)
    save(fig, "fig11_3d_stack")

def surf(ax, X, Y, Z, xl, yl, zl, cmap="viridis", elev=26, az=-58):
    s = ax.plot_surface(X, Y, Z, cmap=cmap, vmin=0, vmax=1.75, edgecolor="none", antialiased=True, rstride=1, cstride=1, alpha=0.97)
    ax.set_xlabel(xl, labelpad=1); ax.set_ylabel(yl, labelpad=1)
    ax.view_init(elev, az); ax.tick_params(labelsize=6.5, pad=-2); ax.set_zlim(0, 1.8)
    return s

def fig12():
    fig = plt.figure(figsize=(7.2, 3.1))
    H = R["heat"]; AR = np.log10(H["AR"]); ti = np.array(H["ti"]); M = np.array(H["MW"])
    X, Y = np.meshgrid(AR, ti); Z = M if M.shape == X.shape else M.T
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    s = surf(ax, X, Y, Z, "log$_{10}$ A$_R$", "t$_i$ (nm)", "window (V)"); ax.set_title("window (V), x$_p$ = 6 V", fontsize=8.5, y=0.97); ax.text2D(0.02, 0.95, "c", transform=ax.transAxes, fontweight="bold", fontsize=10)
    H2 = R["heat_xp"]; AR2 = np.log10(H2["AR"]); xp = np.array(H2["xp"]); M2 = np.array(H2["MW"])
    X2, Y2 = np.meshgrid(AR2, xp); Z2 = M2 if M2.shape == X2.shape else M2.T
    ax = fig.add_subplot(1, 2, 2, projection="3d")
    surf(ax, X2, Y2, Z2, "log$_{10}$ A$_R$", "x$_p$ (V)", "window (V)", az=-125); ax.set_title("window (V), t$_i$ = 10 nm", fontsize=8.5, y=0.97); ax.text2D(0.02, 0.95, "d", transform=ax.transAxes, fontweight="bold", fontsize=10)
    fig.colorbar(s, ax=fig.axes, shrink=0.6, pad=0.02, label="window (V)")
    save(fig, "fig12_3d_design_maps")

def fig13():
    ar = np.logspace(-1.3, 0.6, 34); tf = np.linspace(10, 100, 30)
    Z = np.zeros((len(tf), len(ar))); XS = np.zeros_like(Z)
    for i, t in enumerate(tf):
        for j, a in enumerate(ar):
            Z[i, j] = AnalyticWindow(Params(A_R=a, t_FE=t * 1e-7, xp=6.0)).MW()
    fig = plt.figure(figsize=(7.2, 3.1))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    X, Y = np.meshgrid(np.log10(ar), tf)
    s = ax.plot_surface(X, Y, Z, cmap="magma", edgecolor="none", rstride=1, cstride=1)
    ax.set_xlabel("log$_{10}$ A$_R$", labelpad=1); ax.set_ylabel("t$_{FE}$ (nm)", labelpad=1);     ax.view_init(26, -55); ax.tick_params(labelsize=6.5, pad=-2); ax.set_title("window (V) at x$_p$ = 6 V", fontsize=8.5, y=0.97); ax.text2D(0.02, 0.95, "d", transform=ax.transAxes, fontweight="bold", fontsize=10)
    # drive-efficiency surface: window per volt of drive vs (A_R, x_p)
    xps = np.linspace(2, 12, 30); E = np.zeros((len(xps), len(ar)))
    for i, x in enumerate(xps):
        for j, a in enumerate(ar):
            E[i, j] = AnalyticWindow(Params(A_R=a, xp=x)).MW() / x
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    X2, Y2 = np.meshgrid(np.log10(ar), xps)
    ax2.plot_surface(X2, Y2, E, cmap="cividis", edgecolor="none", rstride=1, cstride=1)
    ax2.set_xlabel("log$_{10}$ A$_R$", labelpad=1); ax2.set_ylabel("x$_p$ (V)", labelpad=1);     ax2.view_init(26, -125); ax2.tick_params(labelsize=6.5, pad=-2); ax2.set_title("window per volt of drive (V/V)", fontsize=8.5, y=0.97); ax2.text2D(0.02, 0.95, "e", transform=ax2.transAxes, fontweight="bold", fontsize=10)
    fig.subplots_adjust(wspace=0.05)
    save(fig, "fig13_3d_thickness_efficiency")

fig11(); fig12(); fig13()
