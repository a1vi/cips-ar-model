import json, numpy as np
from cips_ar import Params, AnalyticWindow
from retention import RetentionWindow
AR = np.logspace(np.log10(0.02), np.log10(4.0), 40)
T = {"1 s": 1.0, "1e4 s": 1e4, "10 yr": 3.15e8}
out = {"AR": AR.tolist(), "scan": [], "time": {}, "ratio": {}}
for xp in (6.0, 10.0):
    for n in (1.0, 2.0):
        for b in (5.0, 10.0, 20.0, 40.0, 80.0):
            for tn, t in T.items():
                mw = np.array([RetentionWindow(Params(A_R=a, xp=xp), b=b, n=n).MW_t(t) for a in AR])
                m0 = np.array([AnalyticWindow(Params(A_R=a, xp=xp)).MW() for a in AR])
                k = int(np.argmax(mw)); plateau = AR[mw >= 0.99 * mw.max()]
                out["scan"].append(dict(xp=xp, n=n, b=b, t=tn, ARstar=float(AR[k]), MWmax=float(mw.max()),
                                        plateau_hi=float(plateau.max()), plateau_lo=float(plateau.min()),
                                        ret_frac_AR0p25=float(RetentionWindow(Params(A_R=0.25, xp=xp), b=b, n=n).MW_t(t) / AnalyticWindow(Params(A_R=0.25, xp=xp)).MW()),
                                        ret_frac_AR1=float(RetentionWindow(Params(A_R=1.0, xp=xp), b=b, n=n).MW_t(t) / AnalyticWindow(Params(A_R=1.0, xp=xp)).MW()),
                                        ret_frac_AR2=float(RetentionWindow(Params(A_R=2.0, xp=xp), b=b, n=n).MW_t(t) / AnalyticWindow(Params(A_R=2.0, xp=xp)).MW())))
# time curves, xp=6, b=20, n=1
ts = np.logspace(-6, 9, 46)
for a in (0.1, 0.25, 0.5, 1.0, 2.0):
    out["time"][str(a)] = [RetentionWindow(Params(A_R=a, xp=6.0), b=20.0).MW_t(t) for t in ts]
out["t"] = ts.tolist()
# surface MW(A_R, log t), b=20
ar2 = np.logspace(-1.3, 0.6, 30); ts2 = np.logspace(-6, 9, 30)
out["surf"] = dict(AR=ar2.tolist(), t=ts2.tolist(), MW=[[RetentionWindow(Params(A_R=a, xp=6.0), b=20.0).MW_t(t) for a in ar2] for t in ts2])
json.dump(out, open("data/retention.json", "w"))
for r in out["scan"]:
    if r["n"] == 1.0 and r["xp"] == 6.0 and r["t"] == "10 yr" or (r["xp"]==10.0 and r["n"]==1.0 and r["t"]=="10 yr"):
        print(r)
