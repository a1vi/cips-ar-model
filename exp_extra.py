"""Extra studies: thickness scaling, drive efficiency, ensemble/step convergence, coercive-field spread, h-BN scaling."""
import json, math, time
import numpy as np
from dataclasses import replace
from cips_ar import Params, window, AnalyticWindow, a_r_sat_closed_form
from exp_ar import x_sat, x_sat_closed

out = {}
# A. ferroelectric thickness scaling
rows = []
for tfe in (10, 15, 20, 30, 40, 60, 87, 100):
    kw = dict(t_FE=tfe * 1e-7)
    mwmax = AnalyticWindow(Params(A_R=1.0, **kw)).MW_max()
    r = {"t_FE_nm": tfe, "MWmax": mwmax}
    for xp in (6.0, 10.0):
        cf = a_r_sat_closed_form(Params(xp=xp, **kw))
        r["ARsat_xp%g" % xp] = cf["A_R_sat"]
    for AR in (0.25, 1.0):
        r["xsat_AR%g" % AR] = x_sat(AR, **kw)
        r["MW_AR%g_xp6" % AR] = AnalyticWindow(Params(A_R=AR, xp=6.0, **kw)).MW()
    rows.append(r)
out["tfe_scaling"] = rows
print("A done", flush=True)

# B. drive efficiency MW/x_sat versus A_R
rows = []
for AR in (0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 4.0):
    xs = x_sat(AR)
    prm = Params(A_R=AR, xp=xs)
    mw = AnalyticWindow(prm).MW()
    rows.append({"AR": AR, "xsat": xs, "MW": mw, "eff": mw / xs})
out["efficiency"] = rows
print("B done", flush=True)

# C. convergence: ensemble size and step (full solver, A_R=1 and 0.25, x_p=6)
rows = []
for N in (400, 1600, 6400, 25600):
    for AR in (0.25, 1.0):
        t = time.time()
        r = window(Params(A_R=AR, xp=6.0, N=N))
        rows.append({"N": N, "AR": AR, "MW": r["MW"], "sec": time.time() - t})
        print("N", N, AR, r["MW"], flush=True)
for fine in (4e-3, 2e-3, 1e-3, 0.5e-3):
    for AR in (0.25, 1.0):
        r = window(Params(A_R=AR, xp=6.0), fine=fine)
        rows.append({"fine_V": fine, "AR": AR, "MW": r["MW"]})
out["convergence"] = rows
print("C done", flush=True)

# D. coercive-field spread sigma (analytic): MW_max, ARsat, MW(A_R=1) at 6 V
rows = []
for s in (0.05, 0.10, 0.15, 0.22, 0.30, 0.40):
    prm = Params(sigma=s)
    cf = a_r_sat_closed_form(prm)
    rows.append({"sigma": s, "MWmax": cf["MW_max"], "ARsat": cf["A_R_sat"], "E_s_over_Ec": cf.get("E_s", float('nan')) / prm.Ec if cf.get("E_s") else None,
                 "MW_AR1": AnalyticWindow(Params(sigma=s, A_R=1.0)).MW()})
out["sigma_scan"] = rows
json.dump(out, open("data/extra.json", "w"), indent=1)
print(json.dumps(out, indent=1)[:6000])
