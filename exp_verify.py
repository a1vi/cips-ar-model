"""E1/E2/E10: verification of Eq. (8), analytic model, reference cross-checks, numerical robustness."""
import json, math, sys
import numpy as np
from multiprocessing import Pool
from dataclasses import replace
from cips_ar import Params, window, AnalyticWindow

X_REF = 6.0  # program drive that reproduces the reference paper's A_R = 0.25 / 1 / 4 windows


def case(args):
    AR, tFE, ti, Dit, CFG, xp = args
    prm = Params(A_R=AR, t_FE=tFE * 1e-7, t_i=ti * 1e-7, Dit=Dit, C_FG=CFG, xp=xp)
    r = window(prm, detail=True)
    if not r["ok"]:
        return None
    out = dict(AR=AR, tFE=tFE, ti=ti, Dit=Dit, CFG=CFG if not math.isinf(CFG) else None, xp=xp,
               MW=r["MW"], eq8=r["eq8"], dP=r["dP"], dQ=r["dQ"],
               jump=max(r["jump_f"], r["jump_b"]))
    if Dit == 0:
        out["analytic"] = AnalyticWindow(prm).MW()
    return out


def e1():
    rng = np.random.default_rng(7)
    cfgs = []
    while len(cfgs) < 300:
        AR = float(np.exp(rng.uniform(np.log(0.05), np.log(4.0))))
        tFE = float(rng.choice([15, 30, 45, 60]))
        ti = float(rng.choice([2, 5, 10, 20, 40]))
        Dit = float(rng.choice([0, 0, 1e11, 1e12, 3e12]))
        CFG = float(rng.choice([math.inf, math.inf, 20e-6, 5e-6, 2e-6]))
        xp = float(rng.uniform(4.0, 12.0))
        cfgs.append((AR, tFE, ti, Dit, CFG, xp))
    with Pool(2) as pool:
        res = pool.map(case, cfgs)
    res = [r for r in res if r is not None and r["MW"] > 0.1]
    resid = np.array([abs(r["eq8"] - r["MW"]) / r["MW"] for r in res])
    ana = [(r["analytic"], r["MW"]) for r in res if "analytic" in r]
    aerr = np.array([abs(a - m) / m for a, m in ana])
    smooth = np.array([abs(r["eq8"] - r["MW"]) / r["MW"] for r in res if r["jump"] < 0.2])
    abrupt = np.array([abs(r["eq8"] - r["MW"]) / r["MW"] for r in res if r["jump"] >= 0.2])
    out = dict(n=len(res), eq8_max=float(resid.max()), eq8_median=float(np.median(resid)),
               eq8_p95=float(np.percentile(resid, 95)),
               n_smooth=int(smooth.size), eq8_smooth_max=float(smooth.max()), eq8_smooth_p95=float(np.percentile(smooth, 95)),
               n_abrupt=int(abrupt.size), eq8_abrupt_max=float(abrupt.max()) if abrupt.size else 0.0,
               n_analytic=len(ana), analytic_max=float(aerr.max()), analytic_median=float(np.median(aerr)),
               analytic_p95=float(np.percentile(aerr, 95)), cases=res)
    return out


def e2():
    rows = {}
    for AR in (0.25, 0.5, 1.0, 2.0, 4.0):
        p = Params(A_R=AR, xp=X_REF)
        rows["AR=%g" % AR] = dict(MW=window(p)["MW"], analytic=AnalyticWindow(p).MW())
    rows["trap_bound"] = {}
    for D in (0, 1e10, 1e11, 1e12, 3e12, 1e13):
        r = window(Params(Dit=D, xp=X_REF), detail=True)
        rows["trap_bound"]["%.0e" % D] = dict(MW=r["MW"], dQ=r["dQ"], dP=r["dP"], eq8=r["eq8"])
    rows["mobility"] = {str(mu): window(Params(mu=mu, xp=X_REF))["MW"] for mu in (2.5, 25.0, 250.0)}
    # published-device stacks, reference-paper model values for comparison
    rows["wang_stack"] = {}
    for ti in (9, 4):
        for AR in (1.0, 0.5):
            p = Params(t_FE=87e-7, t_i=ti * 1e-7, A_R=AR, xp=4.0)
            rows["wang_stack"]["ti=%d,AR=%g" % (ti, AR)] = window(p)["MW"]
    rows["lee_stack"] = {}
    for xp in (3.0, 3.5, 4.0):
        p = Params(t_FE=86.5e-7, t_i=18.6e-7, A_R=1.0, xp=xp)
        rows["lee_stack"]["xp=%g,AR=1" % xp] = window(p)["MW"]
    rows["lee_stack"]["xp=4,AR=0.25"] = window(Params(t_FE=86.5e-7, t_i=18.6e-7, A_R=0.25, xp=4.0))["MW"]
    return rows


def rob_case(args):
    N, seed, AR = args
    p = Params(A_R=AR, xp=X_REF, N=N, sampling="random", seed=seed)
    return window(p)["MW"]


def e10():
    out = {}
    with Pool(2) as pool:
        for AR in (0.25, 1.0, 4.0):
            for N in (400, 1600, 6400):
                vals = pool.map(rob_case, [(N, s, AR) for s in range(20)])
                out["AR=%g,N=%d" % (AR, N)] = dict(mean=float(np.mean(vals)), std=float(np.std(vals, ddof=1)))
    # step size
    out["step"] = {}
    for AR in (0.25, 1.0, 4.0):
        p = Params(A_R=AR, xp=X_REF)
        out["step"]["AR=%g" % AR] = dict(fine_1mV=window(p, fine=1e-3)["MW"], fine_0p25mV=window(p, fine=2.5e-4)["MW"])
    return out


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "e1only":
        res = json.load(open("data/verify.json"))
        res["E1"] = e1()
    else:
        res = dict(E1=e1(), E2=e2(), E10=e10())
    json.dump(res, open("data/verify.json", "w"), indent=1)
    e1r = res["E1"]
    print("E1: n=%d eq8 max %.2e median %.2e p95 %.2e | analytic (n=%d) max %.2e median %.2e p95 %.2e" % (
        e1r["n"], e1r["eq8_max"], e1r["eq8_median"], e1r["eq8_p95"], e1r["n_analytic"], e1r["analytic_max"],
        e1r["analytic_median"], e1r["analytic_p95"]))
    print({k: v for k, v in e1r.items() if k != "cases"})
