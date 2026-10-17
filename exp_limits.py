"""Limit checks (sharp-switching limit), retained-field stress versus A_R, and the large-A_R product MW*A_R.
Writes data/limits.json and data/large_AR.json."""
import json
from cips_ar import Params, AnalyticWindow

out = {"depol": [], "sigma_limit": []}
for a in (0.03, 0.05, 0.1, 0.2, 0.28, 0.5, 1.0, 2.0, 4.0):
    p = Params(A_R=a, xp=6.0)
    w = AnalyticWindow(p)
    r = w.retained()
    out["depol"].append(dict(AR=a, eb_over_Ec=r["eb"] / p.Ec, ef_over_Ec=r["ef"] / p.Ec, Ep_over_Ec=r["Ep"] / p.Ec,
                             atrisk=w.Phi(r["eb"]), MW=w.MW()))
for sg in (0.005, 0.01, 0.02, 0.05, 0.1, 0.22):
    p = Params(A_R=0.25, xp=6.0, sigma=sg)
    mm = AnalyticWindow(p).MW_max()
    out["sigma_limit"].append(dict(sigma=sg, MWmax=mm, ratio_to_2EctFE=mm / (2 * p.Ec * p.t_FE)))
json.dump(out, open("data/limits.json", "w"), indent=1)

big = {}
for xp in (4.0, 6.0, 8.0, 10.0, 12.0):
    rows = []
    for a in (1, 2, 4, 8, 16):
        mw = AnalyticWindow(Params(A_R=a, xp=xp)).MW()
        rows.append([a, round(mw, 4), round(mw * a, 3)])
    big[str(xp)] = rows
json.dump(big, open("data/large_AR.json", "w"))
print("done")
