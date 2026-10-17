"""Sensitivity of the extraction to the assumed program drive mean (source papers state only 'sweeping range of 4 V')."""
import json, copy
import numpy as np
from multiprocessing import Pool
import exp_extract as E
BASE_DEV = copy.deepcopy(E.DEV)
out = {}
for xp in (2.0, 3.0, 4.0, 5.0):
    for k, v in BASE_DEV.items():
        d = dict(v); d["xp"] = (xp, 0.35); E.DEV[k] = d
    res = {}
    with Pool(2) as pool:
        for dev in E.DEV:
            ps = pool.map(E.worker, [(dev, 9000 + i + (0 if dev.startswith("W") else 50), 600) for i in range(2)])
            rows = [r for p in ps for r in p]
            tot = len(rows)
            ar = np.array([r[0] for r in rows if r[0] is not None and np.isfinite(r[0])])
            res[dev] = dict(median=float(np.median(ar)), p16=float(np.percentile(ar, 16)), p84=float(np.percentile(ar, 84)),
                            P_lt_1=float(np.mean(ar < 1)), frac_above_bound=1 - len(ar) / tot, n=int(len(ar)))
    out[str(xp)] = res
    print(xp, res, flush=True)
json.dump(out, open("data/extraction_drive.json", "w"), indent=1)
