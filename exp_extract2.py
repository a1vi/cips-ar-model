"""Prior-sensitivity of the effective-A_R extraction."""
import json, copy
import numpy as np
from multiprocessing import Pool
import exp_extract as E

VARIANTS = {
    "base": {},
    "eps_b 20-30": {"PRIOR": {"eps_b": (20.0, 30.0)}},
    "E_c 200-400": {"PRIOR": {"Ec": (200e3, 400e3)}},
    "x_p sd 0.70 V": {"xp_sd": 0.70},
    "sigma 0.10-0.40": {"PRIOR": {"sigma": (0.10, 0.40)}},
}
BASE_PRIOR = copy.deepcopy(E.PRIOR); BASE_DEV = copy.deepcopy(E.DEV)

def run(pool, name, spec, n_per=600, nw=2):
    E.PRIOR.clear(); E.PRIOR.update(BASE_PRIOR); E.PRIOR.update(spec.get("PRIOR", {}))
    for k, v in BASE_DEV.items():
        d = dict(v)
        if "xp_sd" in spec: d["xp"] = (4.0, spec["xp_sd"])
        E.DEV[k] = d
    res = {}
    for dev in E.DEV:
        parts = pool.map(E.worker, [(dev, 7000 + i + (0 if dev.startswith("W") else 50), n_per) for i in range(nw)])
        rows = [r for p in parts for r in p]
        ar = np.array([r[0] for r in rows if r[0] is not None and np.isfinite(r[0])])
        res[dev] = {"median": float(np.median(ar)), "p16": float(np.percentile(ar, 16)), "p84": float(np.percentile(ar, 84)),
                    "P_lt_1": float(np.mean(ar < 1)), "n": int(len(ar))}
    return res

if __name__ == "__main__":
    # fork-based pool inherits modified module state only if created after modification; create a pool per variant
    out = {}
    for name, spec in VARIANTS.items():
        E.PRIOR.clear(); E.PRIOR.update(BASE_PRIOR); E.PRIOR.update(spec.get("PRIOR", {}))
        for k, v in BASE_DEV.items():
            d = dict(v)
            if "xp_sd" in spec: d["xp"] = (4.0, spec["xp_sd"])
            E.DEV[k] = d
        with Pool(2) as pool:
            parts = {}
            for dev in E.DEV:
                ps = pool.map(E.worker, [(dev, 7000 + i + (0 if dev.startswith("W") else 50), 600) for i in range(2)])
                rows = [r for p in ps for r in p]
                ar = np.array([r[0] for r in rows if r[0] is not None and np.isfinite(r[0])])
                parts[dev] = {"median": float(np.median(ar)), "p16": float(np.percentile(ar, 16)), "p84": float(np.percentile(ar, 84)),
                              "P_lt_1": float(np.mean(ar < 1)), "n": int(len(ar))}
        out[name] = parts
        print(name, parts, flush=True)
    json.dump(out, open("data/extraction_priors.json", "w"), indent=1)
