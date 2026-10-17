"""Inverse extraction of the effective floating-gate area ratio of published CIPS MFMIS devices.

Monte-Carlo over the nuisance parameters; the model window (analytic model, unipolar first-cycle protocol
of the reference paper) is inverted for A_R.  A synthetic recovery test validates the procedure.
"""
import json, math, sys
import numpy as np
from multiprocessing import Pool
from cips_ar import Params, AnalyticWindow, window

PRIOR = dict(
    Ec=(240e3, 360e3),      # V/cm, uniform (Table 1 of the reference: 300 kV/cm +-20%)
    Pr=(4.0e-6, 5.5e-6),    # C/cm^2, uniform (Wang ~4, Lee ~5.4 uC/cm2)
    eps_b=(12.0, 40.0),     # uniform (range used in the reference paper)
    sigma=(0.15, 0.30),     # uniform, assumed
    eps_i=(3.0, 4.5),       # h-BN out-of-plane permittivity, uniform around 3.76
)
DEV = {
    "Wang (MoS2, 87 nm CIPS)": dict(tFE=(87.0, 2.0), ti=("u", 4.0, 9.0), xp=(4.0, 0.35), MW=(3.8, 0.3)),
    "Lee (MoTe2, 86.5 nm CIPS)": dict(tFE=(86.5, 2.0), ti=("n", 18.6, 0.5), xp=(4.0, 0.35), MW=(3.5, 0.4)),
}


def draw_nuis(rng, dev):
    d = DEV[dev]
    n = {}
    for k, (lo, hi) in PRIOR.items():
        n[k] = rng.uniform(lo, hi)
    n["tFE"] = rng.normal(*d["tFE"]) * 1e-7
    if d["ti"][0] == "u":
        n["ti"] = rng.uniform(d["ti"][1], d["ti"][2]) * 1e-7
    else:
        n["ti"] = rng.normal(d["ti"][1], d["ti"][2]) * 1e-7
    n["xp"] = rng.normal(*d["xp"])
    return n


def mk(n, AR):
    return Params(A_R=AR, Ec=n["Ec"], Pr=n["Pr"], eps_b=n["eps_b"], sigma=n["sigma"], eps_i=n["eps_i"],
                  t_FE=n["tFE"], t_i=n["ti"], xp=n["xp"])


def invert(n, MW_obs, lo=0.01, hi=30.0):
    """Return A_R with analytic MW(A_R) = MW_obs; None if above the bound, np.inf if below MW(hi)."""
    mwmax = AnalyticWindow(mk(n, 1.0)).MW_max()
    if MW_obs >= 0.999 * mwmax:
        return None
    if AnalyticWindow(mk(n, hi)).MW() > MW_obs:
        return float("inf")
    a, b = math.log(lo), math.log(hi)
    for _ in range(40):
        m = 0.5 * (a + b)
        if AnalyticWindow(mk(n, math.exp(m))).MW() > MW_obs:
            a = m
        else:
            b = m
    return math.exp(0.5 * (a + b))


def worker(args):
    dev, seed, n = args
    rng = np.random.default_rng(seed)
    res = []
    d = DEV[dev]
    for _ in range(n):
        nu = draw_nuis(rng, dev)
        mwo = rng.normal(*d["MW"])
        ar = invert(nu, mwo)
        res.append((ar, nu["ti"] * 1e7, nu["eps_i"], mwo,
                    AnalyticWindow(mk(nu, 1.0)).MW_max()))
    return res


def summarize(vals):
    v = np.array(vals)
    return dict(median=float(np.median(v)), p16=float(np.percentile(v, 16)), p84=float(np.percentile(v, 84)),
                p2_5=float(np.percentile(v, 2.5)), p97_5=float(np.percentile(v, 97.5)))


def run_extraction(pool, n_per=1500, nw=4):
    out = {}
    for dev in DEV:
        parts = pool.map(worker, [(dev, 1000 + i + (0 if dev.startswith("W") else 50), n_per) for i in range(nw)])
        rows = [r for p in parts for r in p]
        tot = len(rows)
        above = sum(1 for r in rows if r[0] is None)
        ok = [r for r in rows if r[0] is not None and np.isfinite(r[0])]
        ar = [r[0] for r in ok]
        lam = [r[0] * r[1] / r[2] for r in ok]  # A_R t_i/eps_i in nm
        out[dev] = dict(n_total=tot, frac_above_bound=above / tot, n_ok=len(ok),
                        A_R=summarize(ar), Lambda_nm=summarize(lam),
                        frac_AR_below_1=float(np.mean(np.array(ar) < 1.0)),
                        frac_AR_below_0p5=float(np.mean(np.array(ar) < 0.5)),
                        A_R_samples=ar[:4000], Lambda_samples=lam[:4000])
    # nominal-parameter values
    nom = {}
    for ti in (4.0, 6.5, 9.0):
        nu = dict(Ec=3e5, Pr=5.5e-6, eps_b=25.0, sigma=0.22, eps_i=3.76, tFE=87e-7, ti=ti * 1e-7, xp=4.0)
        nom["Wang ti=%g" % ti] = invert(nu, 3.8)
    nu = dict(Ec=3e5, Pr=5.5e-6, eps_b=25.0, sigma=0.22, eps_i=3.76, tFE=86.5e-7, ti=18.6e-7, xp=4.0)
    nom["Lee"] = invert(nu, 3.5)
    out["nominal"] = nom
    return out


# ------------------ synthetic recovery test ----------------------------------------
def synth_worker(args):
    seed, = args
    rng = np.random.default_rng(seed)
    dev = "Wang (MoS2, 87 nm CIPS)"
    # truth
    AR_true = math.exp(rng.uniform(math.log(0.1), math.log(2.0)))
    nu_true = draw_nuis(rng, dev)
    prm = mk(nu_true, AR_true)
    r = window(prm)
    if not r["ok"] or r["MW"] < 0.3:
        return None
    mw_obs = r["MW"] * (1.0 + rng.normal(0, 0.05))
    # inference with prior draws; measurement sd = 5% of the observation
    ars = []
    nf = 0
    ntot = 500
    for _ in range(ntot):
        nu = draw_nuis(rng, dev)
        mo = rng.normal(mw_obs, 0.05 * mw_obs)
        a = invert(nu, mo)
        if a is None:
            nf += 1
        elif np.isfinite(a):
            ars.append(a)
    if len(ars) < 50:
        return None
    ars = np.array(ars)
    return dict(true=AR_true, med=float(np.median(ars)), lo68=float(np.percentile(ars, 16)),
                hi68=float(np.percentile(ars, 84)), lo95=float(np.percentile(ars, 2.5)),
                hi95=float(np.percentile(ars, 97.5)), MW=r["MW"])


if __name__ == "__main__":
    pool = Pool(2)
    res = run_extraction(pool)
    syn = pool.map(synth_worker, [(5000 + i,) for i in range(70)])
    syn = [s for s in syn if s is not None]
    cov68 = float(np.mean([s["lo68"] <= s["true"] <= s["hi68"] for s in syn]))
    cov95 = float(np.mean([s["lo95"] <= s["true"] <= s["hi95"] for s in syn]))
    logerr = np.array([math.log(s["med"] / s["true"]) for s in syn])
    res["synthetic"] = dict(n=len(syn), cov68=cov68, cov95=cov95, median_log_ratio=float(np.median(logerr)),
                            iqr_log=[float(np.percentile(logerr, 25)), float(np.percentile(logerr, 75))],
                            rows=syn)
    json.dump(res, open("data/extraction.json", "w"))
    for k in DEV:
        r = res[k]
        print(k, "n=%d above-bound=%.3f" % (r["n_total"], r["frac_above_bound"]), "A_R", {a: round(b, 3) for a, b in r["A_R"].items()},
              "Lambda", {a: round(b, 2) for a, b in r["Lambda_nm"].items()},
              "P(A_R<1)=%.2f P(<0.5)=%.2f" % (r["frac_AR_below_1"], r["frac_AR_below_0p5"]))
    print("nominal", res["nominal"])
    print("synthetic", {k: v for k, v in res["synthetic"].items() if k != "rows"})
