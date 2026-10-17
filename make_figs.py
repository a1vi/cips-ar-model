import json, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib import cm
from matplotlib.colors import Normalize

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.linewidth": 0.8, "axes.labelsize": 9,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "legend.frameon": False, "legend.fontsize": 7.5, "figure.dpi": 150, "savefig.dpi": 300,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.5,
})
OI = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]
R = json.load(open("data/ar_results.json"))
V = json.load(open("data/verify.json"))
os.makedirs("figures", exist_ok=True)


def save(fig, name):
    fig.savefig("figures/%s.png" % name, bbox_inches="tight")
    fig.savefig("figures/%s.pdf" % name, bbox_inches="tight")
    plt.close(fig)


def plainlog(ax, which="y", ticks=None):
    from matplotlib.ticker import ScalarFormatter, NullFormatter
    a = ax.yaxis if which == "y" else ax.xaxis
    if ticks is not None:
        (ax.set_yticks if which == "y" else ax.set_xticks)(ticks)
    a.set_major_formatter(ScalarFormatter()); a.set_minor_formatter(NullFormatter())


def panel(ax, s):
    ax.text(-0.14, 1.04, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


# ---------------- Fig 1: schematic ----------------
def fig1():
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw=dict(width_ratios=[1.25, 1]))
    ax.set_xlim(0, 10); ax.set_ylim(0, 8); ax.axis("off")
    def box(x, y, w, h, c, label, tc="k", fs=8):
        ax.add_patch(Rectangle((x, y), w, h, fc=c, ec="k", lw=0.7))
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=fs, color=tc)
    box(3.0, 6.6, 4.0, 0.7, "#BBBBBB", "top gate  (V$_G$)")
    box(3.0, 5.0, 4.0, 1.6, "#F4B6C2", "CuInP$_2$S$_6$   t$_{FE}$")
    box(0.8, 4.3, 8.4, 0.7, "#E8C547", "floating gate (neutral)")
    box(0.8, 3.1, 8.4, 1.2, "#A8D8B9", "h-BN interlayer   t$_i$")
    box(0.8, 2.3, 8.4, 0.8, "#7FB2E5", "channel (WSe$_2$)")
    box(0.0, 2.3, 0.8, 0.8, "#DDDDDD", "S", fs=7); box(9.2, 2.3, 0.8, 0.8, "#DDDDDD", "D", fs=7)
    ax.annotate("", xy=(3.0, 0.9), xytext=(7.0, 0.9), arrowprops=dict(arrowstyle="<->", lw=0.8))
    ax.text(5.0, 0.35, "ferroelectric capacitor area  =  A$_R$ × channel capacitor area", ha="center", fontsize=7.5)
    ax.annotate("", xy=(0.8, 1.6), xytext=(9.2, 1.6), arrowprops=dict(arrowstyle="<->", lw=0.8))
    ax.text(5.0, 1.65, "channel capacitor area", ha="center", va="bottom", fontsize=7.5)
    for x0 in (4.0, 5.0, 6.0):
        ax.add_patch(FancyArrowPatch((x0, 5.2), (x0, 6.4), arrowstyle="-|>", mutation_scale=7, color="#8B1A3A", lw=0.9))
    ax.text(7.25, 5.8, "P$_{sw}$", color="#8B1A3A", fontsize=8, va="center")
    ax.text(5.0, 7.85, "D$_{FE}$ = D / A$_R$  =  ε$_0$ε$_b$E$_{FE}$ + P$_{sw}$", ha="center", fontsize=8)
    ax2.axis("off"); ax2.set_xlim(0, 10); ax2.set_ylim(0, 8)
    txt = [
        "Stack equations (gate drive x > 0 accumulates holes)",
        "",
        "D = q (p − n) + Q$_{it}$",
        "D / A$_R$ = ε$_0$ε$_b$ E$_{FE}$ + P$_{sw}$",
        "x = t$_{FE}$ E$_{FE}$ + D (1/C$_i$ + 1/C$_{FG}$) + ψ$_s$",
        "",
        "Read at constant drain current fixes p, ψ$_s$:",
        "MW = t$_{FE}$ (ΔQ$_{it}$/A$_R$ − ΔP$_{sw}$)/(ε$_0$ε$_b$)",
        "       + ΔQ$_{it}$ (1/C$_i$ + 1/C$_{FG}$)",
        "",
        "A$_R$ enters ΔP$_{sw}$ through the program step:",
        "D$_p$ = A$_R$ D$_{FE,p}$ must be supplied by x$_p$",
    ]
    for i, t in enumerate(txt):
        ax2.text(0.0, 7.6 - 0.62 * i, t, fontsize=7.8, va="top", fontweight="bold" if i == 0 else "normal")
    save(fig, "fig1_schematic")


# ---------------- Fig 2: window vs A_R ----------------
def fig2():
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.55))
    ax = axs[0]
    E3 = R["E3"]
    sat = {s["xp"]: s["closed"] for s in R["sat"]}
    for c, xp in zip(OI, ("4.0", "6.0", "8.0", "10.0")):
        d = E3[xp]
        ax.plot(d["AR_an"], d["MW_an"], color=c, lw=1.4, label="x$_p$ = %g V" % float(xp))
        ax.plot(d["AR_solver"], d["MW_solver"], "o", ms=2.8, mfc="none", color=c, mew=0.7)
        if float(xp) in sat:
            ax.plot([sat[float(xp)]], [0.98 * R["MWmax"]], marker="v", color=c, ms=5)
    ax.axhline(R["MWmax"], color="0.4", ls="--", lw=0.8)
    ax.text(5.2, R["MWmax"] + 0.03, "M$_{max}$ = %.2f V" % R["MWmax"], ha="right", va="bottom", fontsize=7.5, color="0.3")
    ax.set_xscale("log"); ax.set_xlim(0.03, 6); ax.set_ylim(0, 1.95)
    ax.set_xlabel("floating-gate area ratio A$_R$"); ax.set_ylabel("memory window (V)"); plainlog(ax, "x", [0.1, 1])
    ax.legend(loc="lower left"); panel(ax, "a")
    ax = axs[1]
    d = E3["6.0"]
    ar = np.array(d["AR_solver"])
    ax.plot(ar, np.array(d["Pb"]) * 1e6, "o-", color=OI[0], ms=3, label="P$_b$ (backward)")
    ax.plot(ar, np.array(d["Pf"]) * 1e6, "s-", color=OI[1], ms=3, label="P$_f$ (forward)")
    ax.axhline(0, color="0.7", lw=0.6)
    ax.set_xscale("log"); ax.set_xlabel("A$_R$"); ax.set_ylabel("switched polarization at threshold (µC/cm$^2$)")
    ax.legend(loc="center left"); ax.set_title("x$_p$ = 6 V", fontsize=8); panel(ax, "b")
    ax = axs[2]
    d = E3["6.0"]
    ax.plot(ar, np.array(d["Eprog"]) / 3e5, "o-", color=OI[2], ms=3)
    ax.axhline(1.32, color="0.4", ls="--", lw=0.8); ax.text(0.035, 1.42, "E$_s$ = 1.32 E$_c$", fontsize=7, color="0.3")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("A$_R$"); ax.set_ylabel("E$_p$ / E$_c$ at program extreme")
    from matplotlib.ticker import ScalarFormatter, NullFormatter
    ax.set_yticks([1, 2, 3, 4, 6]); ax.yaxis.set_major_formatter(ScalarFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_title("x$_p$ = 6 V", fontsize=8); panel(ax, "c")
    fig.tight_layout(w_pad=1.1)
    save(fig, "fig2_window_vs_AR")


# ---------------- Fig 3: collapse, x_sat, A_R_sat(x_p) ----------------
def fig3():
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.55))
    ax = axs[0]
    E4 = R["E4"]
    for c, ti in zip(OI, ("2", "5", "10", "20", "40")):
        d = E4[ti]
        lam = np.array(d["AR"]) * float(ti) / 3.76
        ax.plot(lam, d["MW"], color=c, label="t$_i$ = %s nm" % ti)
    for r in R["collapse_check"]:
        ax.plot([r["lam"]], [r["MW"]], "k.", ms=3)
    ax.set_xscale("log"); ax.set_xlim(0.05, 20); ax.set_ylim(0, 1.95)
    ax.set_xlabel("Λ = A$_R$ t$_i$ / ε$_i$  (nm)"); ax.set_ylabel("memory window (V)")
    ax.legend(loc="lower left"); ax.set_title("x$_p$ = 6 V; dots: full solver", fontsize=8); panel(ax, "a")
    ax = axs[1]
    for c, ti in zip(OI, ("5", "10", "20")):
        d = R["E5_xsat"][ti]
        ar_ = np.array(d["AR"]); xs_ = np.array(d["xsat"]); k = xs_ < 60
        ax.plot(ar_[k], xs_[k], color=c, label="t$_i$ = %s nm" % ti)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("A$_R$"); ax.set_ylabel("x$_{sat}$ for 98% of M$_{max}$ (V)")
    plainlog(ax, "y", [2, 5, 10, 20, 50]); plainlog(ax, "x", [0.03, 0.1, 0.3, 1, 3])
    ax.legend(loc="upper left"); panel(ax, "b")
    ax = axs[2]
    xp = [s["xp"] for s in R["sat"]]
    ax.plot(xp, [s["closed"] for s in R["sat"]], "-", color=OI[0], label="closed form")
    sx = [s for s in R["sat"] if "solver" in s]
    ax.plot([s["xp"] for s in sx], [s["solver"] for s in sx], "o", mfc="none", color=OI[1], label="full solver")
    ax.set_xlabel("program drive x$_p$ (V)"); ax.set_ylabel("A$_R^{sat}$ (98% of M$_{max}$)")
    ax.legend(loc="upper left"); panel(ax, "c")
    fig.tight_layout(w_pad=1.1)
    save(fig, "fig3_collapse_and_drive")


# ---------------- Fig 4: design maps ----------------
def fig4():
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9))
    H = R["heat"]
    AR = np.array(H["AR"]); ti = np.array(H["ti"]); M = np.array(H["MW"])
    ax = axs[0]
    pc = ax.pcolormesh(AR, ti, M, shading="auto", cmap="viridis", vmin=0, vmax=1.75)
    ax.contour(AR, ti, M, levels=[0.98 * R["MWmax"]], colors="w", linewidths=1.2)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("A$_R$"); ax.set_ylabel("h-BN thickness t$_i$ (nm)"); plainlog(ax, "y", [2, 5, 10, 20, 40]); plainlog(ax, "x", [0.1, 1])
    ax.set_title("x$_p$ = 6 V, t$_{FE}$ = 30 nm (white: 98% of M$_{max}$)", fontsize=8); panel(ax, "a")
    cb = fig.colorbar(pc, ax=ax, pad=0.02); cb.set_label("memory window (V)")
    H2 = R["heat_xp"]
    AR2 = np.array(H2["AR"]); xp = np.array(H2["xp"]); M2 = np.array(H2["MW"])
    ax = axs[1]
    pc = ax.pcolormesh(AR2, xp, M2, shading="auto", cmap="viridis", vmin=0, vmax=1.75)
    ax.contour(AR2, xp, M2, levels=[0.98 * R["MWmax"]], colors="w", linewidths=1.2)
    ax.set_xscale("log"); plainlog(ax, "x", [0.1, 1]); ax.set_xlabel("A$_R$"); ax.set_ylabel("program drive x$_p$ (V)")
    ax.set_title("t$_i$ = 10 nm, t$_{FE}$ = 30 nm", fontsize=8); panel(ax, "b")
    cb = fig.colorbar(pc, ax=ax, pad=0.02); cb.set_label("memory window (V)")
    fig.tight_layout(w_pad=1.2)
    save(fig, "fig4_design_maps")


# ---------------- Fig 5: traps ----------------
def fig5():
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.0))
    E = R["E6_fixed_xp"]; AR = np.array(E["AR"])
    labs = {"0e+00": "0", "1e+11": "10$^{11}$", "1e+12": "10$^{12}$", "3e+12": "3×10$^{12}$", "1e+13": "10$^{13}$"}
    ax = axs[0, 0]
    cmap_ = {"0e+00": OI[0], "1e+11": OI[1], "1e+12": OI[2], "3e+12": OI[3], "1e+13": OI[4]}
    for k, v in E["data"].items():
        ax.plot(AR, v["MW"], "-o", ms=2.5, color=cmap_[k], label="D$_{it}$ = " + labs[k])
    ax.set_xscale("log"); ax.set_ylim(-0.4, 1.9); ax.set_xlabel("A$_R$"); ax.set_ylabel("memory window (V)")
    ax.set_title("fixed x$_p$ = 6 V", fontsize=8); ax.legend(loc="lower left", ncol=1); panel(ax, "a")
    ax = axs[0, 1]
    base = np.array(E["data"]["0e+00"]["MW"])
    for k, v in list(E["data"].items())[1:]:
        ax.plot(AR, 100 * (1 - np.array(v["MW"]) / base), "-o", ms=2.5, color={"1e+11": OI[1], "1e+12": OI[2], "3e+12": OI[3], "1e+13": OI[4]}[k], label=labs[k])
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("A$_R$"); ax.set_ylabel("window lost to traps (%)")
    ax.set_title("fixed x$_p$ = 6 V", fontsize=8); ax.legend(title="D$_{it}$ (cm$^{-2}$eV$^{-1}$)", loc="lower right", fontsize=6.5, title_fontsize=6.5); plainlog(ax, "y", [0.5, 1, 5, 10, 50, 100]); panel(ax, "b")
    ax = axs[1, 0]
    E2 = R["E6_scaled_xp"]; AR2 = np.array(E2["AR"])
    colmap = {"0e+00": OI[0], "1e+11": OI[1], "1e+12": OI[2], "3e+12": OI[3], "1e+13": OI[4]}
    for k, v in E2["data"].items():
        c = colmap[k]
        ax.plot(AR2, v["MW"], "-o", ms=2.5, color=c, label="D$_{it}$ = " + {"0e+00": "0", "1e+12": "10$^{12}$", "3e+12": "3×10$^{12}$", "1e+13": "10$^{13}$"}[k])
    ax.set_xscale("log"); ax.set_xlabel("A$_R$"); ax.set_ylabel("memory window (V)")
    ax.set_title("x$_p$ = x$_{sat}$(A$_R$): field held at 1.3 E$_c$", fontsize=8); ax.legend(loc="lower left"); panel(ax, "c")
    ax = axs[1, 1]
    ax.plot(AR, np.array(E["data"]["0e+00"]["Eprog"]) / 3e5, "-o", ms=2.5, color=OI[2], label="fixed x$_p$ = 6 V")
    ax.plot(AR2, np.array(E2["data"]["0e+00"]["Eprog"]) / 3e5, "-s", ms=2.5, color=OI[3], label="x$_p$ = x$_{sat}$")
    ax.axhline(3.3, color="0.4", ls="--", lw=0.8); ax.text(0.04, 3.55, "1 MV/cm", fontsize=7, color="0.3")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("A$_R$"); ax.set_ylabel("E$_p$ / E$_c$ at program extreme")
    plainlog(ax, "y", [1, 2, 3, 5]); ax.legend(loc="upper right"); panel(ax, "d")
    fig.tight_layout(h_pad=1.2, w_pad=1.2)
    save(fig, "fig5_traps")


# ---------------- Fig 6: finite-screening FG and ambipolar cycling ----------------
def fig6():
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7))
    ax = axs[0]
    for c, AR in zip(OI, (0.25, 1.0)):
        rows = [r for r in R["E7"] if r["AR"] == AR]
        base = [r for r in rows if r["C_FG"] is None][0]["MW"]
        cs = [r["C_FG"] * 1e6 for r in rows if r["C_FG"] is not None]
        mw = [100 * (1 - r["MW"] / base) for r in rows if r["C_FG"] is not None]
        ax.plot(cs, mw, "-o", ms=3, color=c, label="A$_R$ = %g" % AR)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_xaxis()
    ax.set_xlabel("floating-gate series capacitance C$_{FG}$ (µF/cm$^2$)"); ax.set_ylabel("window lost (%)")
    ax.axvline(0.333, color="0.4", ls=":", lw=0.8); ax.text(0.37, 0.05, "C$_i$", fontsize=7.5, color="0.3")
    ax.set_title("x$_p$ = 6 V", fontsize=8); ax.legend(loc="upper left"); panel(ax, "a")
    ax = axs[1]
    rows = [r for r in R["E8"]["rows"] if r["Eg"] == 1.65]
    for c, AR in zip(OI, (0.1, 0.25, 0.5, 1.0, 2.0)):
        rr = sorted([r for r in rows if r["AR"] == AR], key=lambda r: r["xe"])
        ax.plot([r["xe"] for r in rr], [r["MW"] for r in rr], "-o", ms=3, color=c, label="A$_R$ = %g" % AR)
        ax.plot([12.7], [R["E8"]["unipolar"][str(AR)]], ">", ms=4, color=c)
    ax.set_xlabel("erase extreme x$_e$ (V)"); ax.set_ylabel("steady-state window (V)")
    ax.set_title("ambipolar channel, E$_g$ = 1.65 eV, x$_p$ = 6 V", fontsize=8)
    ax.text(12.5, 0.15, "▶ first-cycle\n(unipolar) window", fontsize=7, ha="right", color="0.3")
    ax.set_xlim(5.5, 13); ax.legend(loc="lower left", fontsize=6.5, ncol=2); panel(ax, "b")
    fig.tight_layout(w_pad=1.2)
    save(fig, "fig6_screening_and_cycling")


# ---------------- Fig 7: verification ----------------
def fig7():
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5))
    cs = V["E1"]["cases"]
    ax = axs[0]
    mw = np.array([c["MW"] for c in cs]); eq = np.array([c["eq8"] for c in cs])
    ab = np.array([c["jump"] >= 0.2 for c in cs])
    ax.plot(mw[~ab], eq[~ab], ".", ms=3, color=OI[0], label="smooth crossing (n=%d)" % (~ab).sum())
    ax.plot(mw[ab], eq[ab], "x", ms=3.5, color=OI[1], label="abrupt crossing (n=%d)" % ab.sum())
    ax.plot([0, 2], [0, 2], "k-", lw=0.6)
    ax.set_xlabel("simulated window (V)"); ax.set_ylabel("Eq. (1) from P$_{sw}$, Q$_{it}$ (V)"); ax.legend(loc="lower right"); panel(ax, "a")
    ax = axs[1]
    res = np.array([abs(c["eq8"] - c["MW"]) / c["MW"] for c in cs])
    bins = np.logspace(-9, -1.5, 40)
    ax.hist(res[~ab], bins=bins, color=OI[0], alpha=0.85, label="smooth")
    ax.hist(res[ab], bins=bins, color=OI[1], alpha=0.85, label="abrupt")
    ax.set_xscale("log"); ax.set_xlabel("relative residual |Eq. (1) − simulation| / MW"); ax.set_ylabel("cases"); ax.legend(); panel(ax, "b")
    ax = axs[2]
    an = [(c["analytic"], c["MW"]) for c in cs if "analytic" in c]
    a = np.array(an)
    ax.plot(a[:, 1], a[:, 0], ".", ms=3, color=OI[2])
    ax.plot([0, 2], [0, 2], "k-", lw=0.6)
    ax.set_xlabel("simulated window (V)"); ax.set_ylabel("analytic window (V)"); panel(ax, "c")
    fig.tight_layout(w_pad=1.1)
    save(fig, "fig7_verification")


# ---------------- Fig 8: extraction ----------------
def fig8():
    X = json.load(open("data/extraction.json"))
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6))
    for ax, dev, c, lab in zip(axs[:2], ("Wang (MoS2, 87 nm CIPS)", "Lee (MoTe2, 86.5 nm CIPS)"), (OI[0], OI[1]), ("a", "b")):
        s = np.array(X[dev]["A_R_samples"])
        bins = np.geomspace(0.03, 10, 45)
        ax.hist(s, bins=bins, color=c, alpha=0.85)
        ax.set_xscale("log"); ax.set_xlabel("effective A$_R$"); ax.set_ylabel("posterior draws")
        r = X[dev]["A_R"]
        ax.axvline(r["median"], color="k", lw=1)
        ax.axvspan(r["p16"], r["p84"], color="k", alpha=0.12)
        ax.set_title(dev.split(",")[0] + ", " + dev.split(",")[1].split(")")[0].strip() if False else dev.replace("(", "(").replace(")", ")"), fontsize=7.5)
        ax.text(0.97, 0.95, "median %.2f\n68%%: %.2f–%.2f" % (r["median"], r["p16"], r["p84"]), transform=ax.transAxes, ha="right", va="top", fontsize=7)
        panel(ax, lab)
    ax = axs[2]
    syn = X["synthetic"]["rows"]
    t = np.array([s["true"] for s in syn]); m = np.array([s["med"] for s in syn])
    lo = np.array([s["lo68"] for s in syn]); hi = np.array([s["hi68"] for s in syn])
    ax.errorbar(t, m, yerr=[m - lo, hi - m], fmt="o", ms=2.5, color=OI[2], elinewidth=0.5, capsize=0)
    ax.plot([0.08, 2.5], [0.08, 2.5], "k-", lw=0.6)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("true A$_R$ (synthetic)"); ax.set_ylabel("recovered A$_R$ (median, 68%)")
    ax.set_title("synthetic recovery (n=%d): 68%% coverage %.2f" % (X["synthetic"]["n"], X["synthetic"]["cov68"]), fontsize=7.5)
    panel(ax, "c")
    fig.tight_layout(w_pad=1.1)
    save(fig, "fig8_extraction")


# ---------------- Fig 9: chip predictions ----------------
def fig9():
    fig, ax = plt.subplots(1, 1, figsize=(3.6, 2.8))
    for c, xp in zip(OI, (4.0, 6.0, 8.0)):
        rows = [r for r in R["chip"] if r["xp"] == xp]
        ax.plot([r["AR"] for r in rows], [r["MW"] for r in rows], "-o", ms=3, color=c, label="x$_p$ = %g V" % xp)
        ax.plot([r["AR"] for r in rows], [r["MW_trap1e12"] for r in rows], "--s", ms=2.5, color=c, alpha=0.8)
    ax.set_xscale("log"); ax.set_xlabel("A$_R$"); ax.set_ylabel("predicted memory window (V)")
    ax.text(0.97, 0.97, "solid: trap-free\ndashed: D$_{it}$ = 10$^{12}$", transform=ax.transAxes, ha="right", va="top", fontsize=7)
    ax.legend(loc="lower left")
    fig.tight_layout()
    save(fig, "fig9_chip_prediction")


def fig10():
    X = json.load(open("data/extra.json"))
    T = X["tfe_scaling"]; E = X["efficiency"]
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.7))
    ax = axs[0]; panel(ax, "a")
    t = [r["t_FE_nm"] for r in T]
    ax.plot(t, [r["MWmax"] for r in T], "-o", ms=3, color=OI[6], label="MW$_{max}$")
    ax.plot(t, [r["MW_AR0.25_xp6"] for r in T], "-s", ms=3, color=OI[0], label="A$_R$ = 0.25")
    ax.plot(t, [r["MW_AR1_xp6"] for r in T], "-^", ms=3, color=OI[1], label="A$_R$ = 1")
    ax.set_xlabel("t$_{FE}$ (nm)"); ax.set_ylabel("memory window (V)"); ax.legend(loc="upper left")
    ax.text(0.97, 0.05, "x$_p$ = 6 V", transform=ax.transAxes, ha="right", fontsize=7)
    ax = axs[1]; panel(ax, "b")
    ax.plot(t, [r["ARsat_xp6"] for r in T], "-o", ms=3, color=OI[2], label="x$_p$ = 6 V")
    ax.plot(t, [r["ARsat_xp10"] for r in T], "-s", ms=3, color=OI[3], label="x$_p$ = 10 V")
    ax.set_xlabel("t$_{FE}$ (nm)"); ax.set_ylabel("A$_R^{sat}$"); ax.legend(loc="upper right"); ax.set_ylim(0, 0.65)
    ax = axs[2]; panel(ax, "c")
    ax.loglog([r["AR"] for r in E], [r["eff"] for r in E], "-o", ms=3, color=OI[0])
    ax.set_xlabel("A$_R$"); ax.set_ylabel("MW / x$_{sat}$ (V/V)")
    plainlog(ax, "y", [0.03, 0.1, 0.3, 1]); plainlog(ax, "x", [0.03, 0.1, 0.3, 1, 3])
    ax.text(0.97, 0.95, "98% of MW$_{max}$,\n30 nm CIPS", transform=ax.transAxes, ha="right", va="top", fontsize=7)
    fig.tight_layout()
    save(fig, "fig10_thickness_efficiency")


if __name__ == "__main__":
    import sys
    which = sys.argv[1:] or ["1", "2", "3", "4", "5", "6", "7", "8", "9"]
    fns = {"1": fig1, "2": fig2, "3": fig3, "4": fig4, "5": fig5, "6": fig6, "7": fig7, "8": fig8, "9": fig9, "10": fig10}
    for w in which:
        fns[w]()
        print("fig", w, "ok")
