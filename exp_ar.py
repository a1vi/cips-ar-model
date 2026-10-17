"""Main experiments: A_R sweeps, collapse, design maps, traps, finite-screening FG, ambipolar cycling, chip predictions."""
import json, math, sys
import numpy as np
from multiprocessing import Pool
from dataclasses import replace
from cips_ar import Params, window, AnalyticWindow, a_r_sat_closed_form, EPS0, Q

XREF = 6.0
GAMMA = 0.98


def mw_solver(args):
    kw, kind = args
    prm = Params(**kw)
    if kind == "amb":
        return window(prm, cycles=3, detail=True)
    return window(prm, detail=True)


def an_mw(**kw):
    return AnalyticWindow(Params(**kw)).MW()


def x_sat(AR, gamma=GAMMA, **kw):
    """Smallest program drive giving gamma * MW_max (analytic model, bisection on x_p)."""
    prm = Params(A_R=AR, **kw)
    target = gamma * AnalyticWindow(prm).MW_max()
    lo, hi = 0.3, 80.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if AnalyticWindow(replace(prm, xp=mid)).MW() >= target:
            hi = mid
        else:
            lo = mid
    return hi


def x_sat_closed(AR, gamma=GAMMA, **kw):
    """Closed-form x_sat: x_sat = t_FE E_s + c A_R D_FE,s + psi_p."""
    prm = Params(A_R=AR, **kw)
    cf = a_r_sat_closed_form(prm, gamma)
    aw = AnalyticWindow(prm)
    Es = cf["E_s"]
    Dfe = aw.e0eb * Es + prm.Pr * (2 * aw.Phi(Es) - 1)
    psi = aw.ch.u_of_p(max(AR * Dfe / Q, 1e-30))
    return prm.t_FE * Es + prm.c_series * AR * Dfe + psi


def a_r_sat_num(gamma=GAMMA, **kw):
    prm = Params(**kw)
    target = gamma * AnalyticWindow(prm).MW_max()
    lo, hi = 0.003, 50.0
    for _ in range(70):
        mid = math.sqrt(lo * hi)
        if AnalyticWindow(replace(prm, A_R=mid)).MW() > target:
            lo = mid
        else:
            hi = mid
    return lo


def a_r_sat_solver(gamma=GAMMA, **kw):
    prm = Params(**kw)
    target = gamma * AnalyticWindow(prm).MW_max()
    lo, hi = 0.02, 5.0
    for _ in range(14):
        mid = math.sqrt(lo * hi)
        if window(replace(prm, A_R=mid))["MW"] > target:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def main():
    out = {}
    pool = Pool(2)
    # ---------- E3: MW vs A_R for several drives -------------
    ARg = np.geomspace(0.03, 5.0, 28)
    ARf = np.geomspace(0.01, 8.0, 160)
    e3 = {}
    for xp in (4.0, 6.0, 8.0, 10.0):
        sol = pool.map(mw_solver, [(dict(A_R=float(a), xp=xp), "uni") for a in ARg])
        e3[str(xp)] = dict(AR_solver=ARg.tolist(), MW_solver=[s["MW"] for s in sol],
                           Pf=[s["Pf"] for s in sol], Pb=[s["Pb"] for s in sol],
                           Eprog=[s["E_prog"] for s in sol], uprog=[s["u_prog"] for s in sol],
                           AR_an=ARf.tolist(), MW_an=[an_mw(A_R=float(a), xp=xp) for a in ARf],
                           MWmax=AnalyticWindow(Params(xp=xp)).MW_max())
    out["E3"] = e3
    # ---------- A_R_sat vs xp: closed form vs analytic-numeric vs full solver ------
    xps = [3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0]
    sat = []
    for xp in xps:
        cf = a_r_sat_closed_form(Params(xp=xp), GAMMA)
        sat.append(dict(xp=xp, closed=cf["A_R_sat"], numeric=a_r_sat_num(xp=xp), Es=cf["E_s"], psi=cf["psi_p"]))
    for xp in (4.0, 6.0, 8.0, 12.0):
        r = [s for s in sat if s["xp"] == xp][0]
        r["solver"] = a_r_sat_solver(xp=xp)
    out["sat"] = sat
    out["MWmax"] = AnalyticWindow(Params()).MW_max()
    # E* (kappa)
    aw = AnalyticWindow(Params(A_R=1e9))
    out["kappa"] = aw.forward_field() / Params().Ec
    # ---------- E4: collapse ---------------------------------------------------------
    e4 = {}
    for ti in (2, 5, 10, 20, 40):
        ARs = np.geomspace(0.02, 8.0, 80)
        e4[str(ti)] = dict(AR=ARs.tolist(), MW=[an_mw(A_R=float(a), t_i=ti * 1e-7, xp=XREF) for a in ARs])
    # solver spot checks on the collapse
    coll = []
    for lam_ar in (0.25, 1.0):
        for ti in (2, 5, 10, 20, 40):
            AR = lam_ar * 10.0 / ti
            coll.append(dict(lam=lam_ar * 10 / 3.76, ti=ti, AR=AR, MW=window(Params(A_R=AR, t_i=ti * 1e-7, xp=XREF))["MW"]))
    out["E4"] = e4
    out["collapse_check"] = coll
    # ---------- x_sat vs A_R -----------------------------------------------------------
    e5 = {}
    for ti in (5, 10, 20):
        ARs = np.geomspace(0.02, 3.0, 40)
        e5[str(ti)] = dict(AR=ARs.tolist(),
                           xsat=[x_sat(float(a), t_i=ti * 1e-7) for a in ARs],
                           xsat_closed=[x_sat_closed(float(a), t_i=ti * 1e-7) for a in ARs])
    out["E5_xsat"] = e5
    # heat map MW(A_R, t_i) at XREF (analytic) and A_R_sat curve
    ARh = np.geomspace(0.03, 6.0, 70)
    tih = np.geomspace(2.0, 40.0, 60)
    H = [[an_mw(A_R=float(a), t_i=float(t) * 1e-7, xp=XREF) for a in ARh] for t in tih]
    out["heat"] = dict(AR=ARh.tolist(), ti=tih.tolist(), MW=H,
                       sat=[a_r_sat_num(t_i=float(t) * 1e-7, xp=XREF) for t in tih])
    # heat map MW(A_R, x_p) for t_i = 10 nm
    xph = np.linspace(2.0, 14.0, 49)
    H2 = [[an_mw(A_R=float(a), xp=float(x)) for a in ARh] for x in xph]
    out["heat_xp"] = dict(AR=ARh.tolist(), xp=xph.tolist(), MW=H2)
    # ---------- E6: traps ---------------------------------------------------------------
    ARt = np.geomspace(0.04, 4.0, 22)
    tr = {}
    for D in (0.0, 1e11, 1e12, 3e12, 1e13):
        sol = pool.map(mw_solver, [(dict(A_R=float(a), Dit=D, xp=XREF), "uni") for a in ARt])
        tr["%.0e" % D] = dict(MW=[s["MW"] for s in sol], Eprog=[s["E_prog"] for s in sol],
                              dQ=[s["dQ"] for s in sol], dP=[s["dP"] for s in sol],
                              uprog=[s["u_prog"] for s in sol])
    out["E6_fixed_xp"] = dict(AR=ARt.tolist(), data=tr)
    # scaled drive: x_p = x_sat(A_R) (the drive that just saturates the trap-free window)
    ARs2 = np.geomspace(0.04, 3.0, 20)
    xs2 = [x_sat(float(a)) for a in ARs2]
    tr2 = {}
    for D in (0.0, 1e12, 3e12, 1e13):
        sol = pool.map(mw_solver, [(dict(A_R=float(a), Dit=D, xp=float(x)), "uni") for a, x in zip(ARs2, xs2)])
        tr2["%.0e" % D] = dict(MW=[s["MW"] for s in sol], Eprog=[s["E_prog"] for s in sol])
    out["E6_scaled_xp"] = dict(AR=ARs2.tolist(), xp=xs2, data=tr2)
    # ---------- E7: finite-screening floating gate ---------------------------------------
    cfg = [math.inf, 100e-6, 50e-6, 20e-6, 10e-6, 5e-6, 2e-6, 1e-6, 0.5e-6]
    e7 = []
    for AR in (0.25, 1.0):
        sol = pool.map(mw_solver, [(dict(A_R=AR, C_FG=C, xp=XREF), "uni") for C in cfg])
        for C, s in zip(cfg, sol):
            e7.append(dict(AR=AR, C_FG=None if math.isinf(C) else C, MW=s["MW"], eq8=s["eq8"]))
    out["E7"] = e7
    # ---------- E8: ambipolar steady state ------------------------------------------------
    e8 = []
    jobs = []
    for AR in (0.1, 0.25, 0.5, 1.0, 2.0):
        for xe in (6.0, 8.0, 10.0, 12.0):
            jobs.append(dict(A_R=AR, xp=XREF, xe=xe, ambipolar=True, Eg=1.65))
    sol = pool.map(mw_solver, [(j, "amb") for j in jobs])
    for j, s in zip(jobs, sol):
        e8.append(dict(AR=j["A_R"], xe=j["xe"], Eg=j["Eg"], MW=s["MW"], Pf=s["Pf"], Pb=s["Pb"]))
    jobs = [dict(A_R=AR, xp=XREF, xe=10.0, ambipolar=True, Eg=Eg) for AR in (0.25, 1.0) for Eg in (1.0, 1.65, 2.0)]
    sol = pool.map(mw_solver, [(j, "amb") for j in jobs])
    for j, s in zip(jobs, sol):
        e8.append(dict(AR=j["A_R"], xe=j["xe"], Eg=j["Eg"], MW=s["MW"], Pf=s["Pf"], Pb=s["Pb"]))
    uni = {str(a): window(Params(A_R=a, xp=XREF))["MW"] for a in (0.1, 0.25, 0.5, 1.0, 2.0)}
    out["E8"] = dict(rows=e8, unipolar=uni)
    # ---------- E11: chip predictions -------------------------------------------------------
    chip = []
    for xp in (4.0, 6.0, 8.0):
        for AR in (0.1, 0.25, 0.5, 1.0, 2.0, 4.0):
            r0 = window(Params(A_R=AR, xp=xp))["MW"]
            r1 = window(Params(A_R=AR, xp=xp, Dit=1e12))["MW"]
            chip.append(dict(xp=xp, AR=AR, MW=r0, MW_trap1e12=r1))
    out["chip"] = chip
    # ---------- sensitivity of the A_R curve to ferroelectric parameters ----------------------
    sens = []
    base = dict(xp=XREF)
    for name, key in (("Ec", "Ec"), ("Pr", "Pr"), ("eps_b", "eps_b"), ("sigma", "sigma"), ("t_FE", "t_FE")):
        b0 = getattr(Params(), key)
        for f in (0.8, 1.2):
            kw = dict(base)
            kw[key] = b0 * f
            row = dict(param=name, factor=f,
                       MWmax=AnalyticWindow(Params(**kw)).MW_max(),
                       MW_AR025=an_mw(A_R=0.25, **kw), MW_AR1=an_mw(A_R=1.0, **kw),
                       ARsat=a_r_sat_closed_form(Params(**kw), GAMMA)["A_R_sat"])
            sens.append(row)
    # interlayer permittivity / thickness
    for name, kw2 in (("t_i x0.8", dict(t_i=8e-7)), ("t_i x1.2", dict(t_i=12e-7)),
                      ("eps_i x0.8", dict(eps_i=3.76 * 0.8)), ("eps_i x1.2", dict(eps_i=3.76 * 1.2))):
        kw = dict(base)
        kw.update(kw2)
        sens.append(dict(param=name, factor=None, MWmax=AnalyticWindow(Params(**kw)).MW_max(),
                         MW_AR025=an_mw(A_R=0.25, **kw), MW_AR1=an_mw(A_R=1.0, **kw),
                         ARsat=a_r_sat_closed_form(Params(**kw), GAMMA)["A_R_sat"]))
    sens.append(dict(param="baseline", factor=1.0, MWmax=AnalyticWindow(Params(**base)).MW_max(),
                     MW_AR025=an_mw(A_R=0.25, **base), MW_AR1=an_mw(A_R=1.0, **base),
                     ARsat=a_r_sat_closed_form(Params(**base), GAMMA)["A_R_sat"]))
    out["sens"] = sens
    json.dump(out, open("data/ar_results.json", "w"))
    print("done")


if __name__ == "__main__":
    main()
