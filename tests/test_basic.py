"""Quick regression tests (a few seconds). Run from the repository root:  python -m pytest tests  or  python tests/test_basic.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from cips_ar import Params, AnalyticWindow, a_r_sat_closed_form
from retention import RetentionWindow


def test_window_plateau_and_decay():
    w = lambda a: AnalyticWindow(Params(A_R=a, xp=6.0)).MW()
    assert abs(w(0.1) - 1.742) < 0.01          # fully programmed plateau
    assert abs(w(1.0) - 1.443) < 0.01          # K/A_R branch
    assert w(4.0) < w(1.0) < w(0.25)


def test_saturation_ratio_closed_form():
    r = a_r_sat_closed_form(Params(xp=6.0))
    assert abs(r["A_R_sat"] - 0.282) < 0.01


def test_sharp_switching_limit():
    p = Params(A_R=0.25, xp=6.0, sigma=0.005)
    assert abs(AnalyticWindow(p).MW_max() / (2 * p.Ec * p.t_FE) - 1) < 0.002


def test_retention_reduces_to_quasistatic_at_tau0():
    p = Params(A_R=0.25, xp=6.0)
    assert abs(RetentionWindow(p).MW_t(1e-6) - AnalyticWindow(p).MW()) < 1e-9
    assert RetentionWindow(p, b=20).MW_t(3.15e8) < AnalyticWindow(p).MW()


if __name__ == "__main__":
    for n, f in list(globals().items()):
        if n.startswith("test_"):
            f(); print("ok", n)
