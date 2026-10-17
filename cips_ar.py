"""
cips_ar.py -- self-consistent MFMIS FeFET model for the floating-gate area-ratio study.

Independent re-implementation from the published model description of
Mahim & Rahman, Mater. Today Commun. 57 (2026) 116333 (their Sections 2-3), extended with
an explicit floating-gate area ratio A_R, a finite-screening floating gate, and a
semi-analytic programming-limited window model.

Units: cm, V, C, F.  Gate drive x = -(V_G - V_FB); x > 0 accumulates holes.

Stack equations (A_R = ferroelectric area / channel area):
    D      = q (p - n) + Q_it                      displacement on the channel side
    D_FE   = D / A_R = eps0 eps_b E_FE + P         floating gate neutral
    x      = t_FE E_FE + D (1/C_i + 1/C_FG) + psi_s(p, n)
"""
import math
from dataclasses import dataclass, replace

import numpy as np
from scipy.special import ndtr, ndtri

Q = 1.602176634e-19
KB = 1.380649e-23
HBAR = 1.054571817e-34
M0 = 9.1093837015e-31
EPS0 = 8.8541878128e-14  # F/cm
TEMP = 300.0
VT = KB * TEMP / Q


def dos_cm2(m_rel, g):
    """g m kT / (2 pi hbar^2) in cm^-2."""
    return g * m_rel * M0 * KB * TEMP / (2.0 * math.pi * HBAR ** 2) * 1e-4


def _sp(z):
    if z > 35.0:
        return z + math.exp(-z)
    return math.log1p(math.exp(z))


def _sig(z):
    if z >= 0.0:
        return 1.0 / (1.0 + math.exp(-z))
    e = math.exp(z)
    return e / (1.0 + e)


@dataclass(frozen=True)
class Params:
    # ferroelectric (Table 1 of the reference paper)
    t_FE: float = 30e-7          # cm
    eps_b: float = 25.0
    Pr: float = 5.5e-6           # C/cm^2
    Ec: float = 3.0e5            # V/cm
    sigma: float = 0.22          # relative Gaussian spread of coercive fields
    N: int = 6400                # hysterons
    # interlayer
    t_i: float = 10e-7           # cm
    eps_i: float = 3.76
    # geometry / floating gate
    A_R: float = 1.0
    C_FG: float = math.inf       # F/cm^2, series capacitance of an imperfectly screening FG
    # interface traps (hold-and-release bound)
    Dit: float = 0.0             # cm^-2 eV^-1
    # drive
    xp: float = 6.6              # program extreme (V)
    xe: float = 6.6              # erase extreme magnitude (V), used when ambipolar
    # channel (monolayer WSe2, 2-valley hole model)
    mK: float = 0.40
    gK: int = 2
    mG: float = 2.2
    gG: int = 1
    dEG: float = 0.157           # eV, Gamma below K
    ambipolar: bool = False
    Eg: float = 1.65             # eV
    me: float = 0.45
    ge: int = 2
    # read criterion
    I_per_W: float = 1e-7 / 1e-4  # A/cm  (100 nA/um)
    L: float = 2e-4              # cm
    Vds: float = 1.0
    mu: float = 25.0             # cm^2/Vs
    Rc_W: float = 50e3 * 1e-4    # ohm cm  (50 kOhm um)
    sampling: str = "quantile"
    seed: int = 0

    @property
    def b(self):
        return self.t_FE / (EPS0 * self.eps_b)

    @property
    def C_i(self):
        return EPS0 * self.eps_i / self.t_i

    @property
    def c_series(self):
        return 1.0 / self.C_i + (0.0 if math.isinf(self.C_FG) else 1.0 / self.C_FG)

    @property
    def a(self):
        return self.b / self.A_R + self.c_series

    @property
    def p_th(self):
        veff = self.Vds - self.I_per_W * self.Rc_W
        return self.I_per_W * self.L / (Q * self.mu * veff)


class Channel:
    def __init__(self, prm: Params):
        self.NK = dos_cm2(prm.mK, prm.gK)
        self.NG = dos_cm2(prm.mG, prm.gG)
        self.dG = prm.dEG / VT
        self.amb = prm.ambipolar
        self.Ne = dos_cm2(prm.me, prm.ge)
        self.eg = prm.Eg / VT

    def pn(self, u):
        z = u / VT
        p = self.NK * _sp(z) + self.NG * _sp(z - self.dG)
        dp = (self.NK * _sig(z) + self.NG * _sig(z - self.dG)) / VT
        n = dn = 0.0
        if self.amb:
            zn = (-u / VT) - self.eg
            n = self.Ne * _sp(zn)
            dn = -self.Ne * _sig(zn) / VT
        return p, dp, n, dn

    def u_of_p(self, p_target):
        """Surface potential for a given hole density (interpolated table, hole-only)."""
        key = (self.NK, self.NG, self.dG)
        tab = _UTAB.get(key)
        if tab is None:
            us = np.linspace(-1.2, 1.5, 20001)
            ps = np.array([self.NK * _sp(u / VT) + self.NG * _sp(u / VT - self.dG) for u in us])
            tab = (np.log(ps), us)
            _UTAB[key] = tab
        return float(np.interp(math.log(max(p_target, 1e-300)), tab[0], tab[1]))


_UTAB = {}


def coercive_fields(prm: Params):
    n = prm.N
    if prm.sampling == "quantile":
        z = ndtri((np.arange(n) + 0.5) / n)
    else:
        rng = np.random.default_rng(prm.seed)
        z = np.sort(rng.standard_normal(n))
    ec = prm.Ec * (1.0 + prm.sigma * z)
    return np.maximum(ec, 0.02 * prm.Ec)  # ascending


class Stack:
    """Quasi-static solver. State: hysteron signs, running max of surface potential."""

    def __init__(self, prm: Params):
        self.prm = prm
        self.ch = Channel(prm)
        self.Ec = coercive_fields(prm)
        self.N = prm.N
        self.s = -np.ones(prm.N, dtype=np.int8)
        self.nup = 0
        self.u = -1.0
        self.umax = -1e9
        self.psi_th = self.ch.u_of_p(prm.p_th)
        self.qd = Q * prm.Dit
        self.e0eb = EPS0 * prm.eps_b

    # ---- electrostatics -------------------------------------------------
    def _eval(self, u, T, umax_prev):
        prm = self.prm
        p, dp, n, dn = self.ch.pn(u)
        D = Q * (p - n)
        dD = Q * (dp - dn)
        if self.qd > 0.0:
            um = u if u > umax_prev else umax_prev
            D += self.qd * max(0.0, um - self.psi_th)
            if u > umax_prev and u > self.psi_th:
                dD += self.qd
        return prm.a * D + u - T, prm.a * dD + 1.0, D, p

    def solve(self, T, u0, umax_prev):
        f, df, D, p = self._eval(u0, T, umax_prev)
        if f == 0.0:
            return u0, D, p
        if f > 0:
            hi, fhi = u0, f
            step = 0.5
            lo = hi - step
            flo = self._eval(lo, T, umax_prev)[0]
            while flo > 0:
                hi, fhi = lo, flo
                step *= 2
                lo = hi - step
                flo = self._eval(lo, T, umax_prev)[0]
        else:
            lo, flo = u0, f
            step = 0.5
            hi = lo + step
            fhi = self._eval(hi, T, umax_prev)[0]
            while fhi < 0:
                lo, flo = hi, fhi
                step *= 2
                hi = lo + step
                fhi = self._eval(hi, T, umax_prev)[0]
        u = u0 if lo < u0 < hi else 0.5 * (lo + hi)
        for _ in range(200):
            f, df, D, p = self._eval(u, T, umax_prev)
            if f > 0:
                hi = u
            else:
                lo = u
            un = u - f / df
            if not (lo < un < hi):
                un = 0.5 * (lo + hi)
            if abs(un - u) < 1e-13 or (hi - lo) < 1e-13:
                u = un
                break
            u = un
        f, df, D, p = self._eval(u, T, umax_prev)
        return u, D, p

    def _P(self, nup):
        return self.prm.Pr * (2.0 * nup / self.N - 1.0)

    def _field(self, D, P):
        return (D / self.prm.A_R - P) / self.e0eb

    # ---- hysteron relaxation at fixed gate drive -------------------------
    def settle(self, x):
        prm = self.prm
        b = prm.b
        Pr = prm.Pr
        for _ in range(200):
            P = self._P(self.nup)
            self.u, D, p = self.solve(x + b * P, self.u, self.umax)
            E = self._field(D, P)
            if E > 0:
                cand = np.flatnonzero(self.s < 0)
                if cand.size == 0 or self.Ec[cand[0]] >= E:
                    break
                lo, hi = 1, cand.size
                while lo < hi:
                    mid = (lo + hi + 1) // 2
                    Pm = self._P(self.nup + mid - 1)
                    um, Dm, _ = self.solve(x + b * Pm, self.u, self.umax)
                    if self._field(Dm, Pm) > self.Ec[cand[mid - 1]]:
                        lo = mid
                    else:
                        hi = mid - 1
                self.s[cand[:lo]] = 1
                self.nup += lo
            elif E < 0:
                cand = np.flatnonzero(self.s > 0)
                if cand.size == 0 or self.Ec[cand[0]] >= -E:
                    break
                lo, hi = 1, cand.size
                while lo < hi:
                    mid = (lo + hi + 1) // 2
                    Pm = self._P(self.nup - (mid - 1))
                    um, Dm, _ = self.solve(x + b * Pm, self.u, self.umax)
                    if self._field(Dm, Pm) < -self.Ec[cand[mid - 1]]:
                        lo = mid
                    else:
                        hi = mid - 1
                self.s[cand[:lo]] = -1
                self.nup -= lo
            else:
                break
        P = self._P(self.nup)
        self.u, D, p = self.solve(x + b * P, self.u, self.umax)
        if self.u > self.umax:
            self.umax = self.u
        # trapped charge after the update of the running maximum
        Qit = self.qd * max(0.0, self.umax - self.psi_th)
        n = self.ch.pn(self.u)[2]
        return dict(x=x, u=self.u, p=p, n=n, P=P, D=D, Qit=Qit, E=self._field(D, P))

    def reset_traps(self):
        """Full depletion: traps emptied."""
        self.umax = -1e9


def _interp_cross(rec_a, rec_b, lnp_th):
    la, lb = math.log(rec_a["p"]), math.log(rec_b["p"])
    w = (lnp_th - la) / (lb - la)
    out = {}
    for k in ("x", "P", "Qit", "u", "D"):
        out[k] = rec_a[k] + w * (rec_b[k] - rec_a[k])
    out["jump_lnp"] = abs(lb - la)      # change of ln p across the bracketing step
    out["jump_P"] = abs(rec_b["P"] - rec_a["P"])
    return out


def _ramp(stack, x0, x1, lnp_th, want, fine=1e-3, coarse=0.05, band=1.2, cross_dir=None):
    """Ramp x0 -> x1; return crossing record of ln p through lnp_th in direction cross_dir."""
    x = x0
    prev = stack.settle(x)
    cross = None
    sgn = 1.0 if x1 > x0 else -1.0
    while (x1 - x) * sgn > 1e-12:
        lp = math.log(max(prev["p"], 1e-300))
        dx = fine if abs(lp - lnp_th) < band else coarse
        x = x + sgn * min(dx, abs(x1 - x))
        cur = stack.settle(x)
        if want and cross is None:
            lc = math.log(max(cur["p"], 1e-300))
            if cross_dir == "up" and lp < lnp_th <= lc:
                cross = _interp_cross(prev, cur, lnp_th)
            if cross_dir == "down" and lp > lnp_th >= lc:
                cross = _interp_cross(prev, cur, lnp_th)
        prev = cur
    return cross, prev


def window(prm: Params, cycles=None, detail=False, fine=1e-3, coarse=0.05):
    """Constant-current memory window.

    Unipolar protocol of the reference paper (cycles=None): start from the erased state (all hysterons
    down) in full depletion, ramp to +xp and back to -xp; forward crossing on the up-ramp, backward crossing
    on the down-ramp.  Ambipolar steady-state protocol (cycles=n): sweep -xe -> +xp -> -xe n times and use the
    last loop.
    """
    st = Stack(prm)
    lnp = math.log(prm.p_th)
    if cycles is None:
        x_start = -prm.xp
        st.u = -1.0
        st.settle(x_start)
        cf, last = _ramp(st, x_start, prm.xp, lnp, True, fine=fine, coarse=coarse, cross_dir="up")
        prog = dict(last)
        cb, end = _ramp(st, prm.xp, x_start, lnp, True, fine=fine, coarse=coarse, cross_dir="down")
    else:
        x_start = -prm.xe
        st.settle(x_start)
        for c in range(cycles):
            st.reset_traps()
            cf, last = _ramp(st, x_start, prm.xp, lnp, c == cycles - 1, fine=fine, coarse=coarse, cross_dir="up")
            prog = dict(last)
            cb, end = _ramp(st, prm.xp, x_start, lnp, c == cycles - 1, fine=fine, coarse=coarse, cross_dir="down")
    if cf is None or cb is None:
        res = dict(MW=float("nan"), ok=False)
        return res
    MW = cf["x"] - cb["x"]
    if not detail:
        return dict(MW=MW, ok=True)
    dP = cf["P"] - cb["P"]
    dQ = cf["Qit"] - cb["Qit"]
    # Eq. (8): interlayer term uses the total series compliance c_series = 1/C_i + 1/C_FG
    eq8 = prm.t_FE * (dQ / prm.A_R - dP) / (EPS0 * prm.eps_b) + dQ * prm.c_series
    return dict(MW=MW, ok=True, Pf=cf["P"], Pb=cb["P"], dP=dP, dQ=dQ, eq8=eq8,
                xf=cf["x"], xb=cb["x"], jump_f=cf["jump_lnp"], jump_b=cb["jump_lnp"], jumpP_f=cf["jump_P"], jumpP_b=cb["jump_P"], uf=cf["u"], ub=cb["u"],
                P_prog=prog["P"], u_prog=prog["u"], D_prog=prog["D"], E_prog=prog["E"],
                nup_prog_frac=None, psi_th=st.psi_th, umax=st.umax)


# ---------------------------------------------------------------------------------
# Semi-analytic programming-limited model
# ---------------------------------------------------------------------------------
class AnalyticWindow:
    """Closed-form-style window for the unipolar first-cycle protocol.

    Programming balance at the extreme:
        xp = t_FE E_p + D_p c + A_R-free:  D_p = A_R (eps0 eps_b E_p + P_p),   P_p = Pr (2 Phi(E_p) - 1)
    Relaxation at the read threshold (D_th ~ q p_th + Q_it):
        backward:  eps0 eps_b e_b = P_b - D_th/A_R,   P_b = Pr (2 [Phi(E_p) - Phi(e_b)] - 1)
        forward :  eps0 eps_b e_f = D_th/A_R - P_f,   P_f = Pr (2 Phi(e_f) - 1)
    MW = t_FE (P_b - P_f) / (eps0 eps_b)     (Eq. (8) with dQ_it = 0)
    Phi(E) = normal CDF of (E/Ec - 1)/sigma.
    """

    def __init__(self, prm: Params):
        self.prm = prm
        self.ch = Channel(prm)
        self.e0eb = EPS0 * prm.eps_b
        self.D_th = Q * prm.p_th

    def Phi(self, E):
        z = (E / self.prm.Ec - 1.0) / self.prm.sigma
        return 0.5 * math.erfc(-z / math.sqrt(2.0))

    def _bisect(self, f, lo, hi, n=64):
        flo = f(lo)
        for _ in range(n):
            mid = 0.5 * (lo + hi)
            fm = f(mid)
            if (fm > 0) == (flo > 0):
                lo, flo = mid, fm
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def forward_field(self):
        """Depolarization-limited field at the forward crossing (erased start, D = D_th)."""
        prm = self.prm

        def ff(e):
            Pf = prm.Pr * (2.0 * self.Phi(e) - 1.0)
            return self.e0eb * e - (self.D_th / prm.A_R - Pf)

        return self._bisect(ff, 0.0, 3.0e6)

    def program_field(self):
        """E_p from the programming balance (hole-only channel); never below the forward relaxed field,
        because the up-ramp already passes through E_f."""
        prm = self.prm
        Ef = self.forward_field()

        def g(E):
            Pp = prm.Pr * (2.0 * self.Phi(E) - 1.0)
            D = prm.A_R * (self.e0eb * E + Pp)
            p = max(D / Q, 1e-30)
            u = self.ch.u_of_p(p)
            return prm.xp - (prm.t_FE * E + D * prm.c_series + u)

        if g(Ef) <= 0:
            return Ef
        return self._bisect(g, Ef, 3.0e6)

    def retained(self):
        prm = self.prm
        Ep = self.program_field()
        Php = self.Phi(Ep)
        Pp = prm.Pr * (2.0 * Php - 1.0)
        Ef = self.forward_field()
        Pf = prm.Pr * (2.0 * self.Phi(Ef) - 1.0)
        if Pp <= self.D_th / prm.A_R:
            # depolarization field does not push against the programmed state: no back-switching
            return dict(Ep=Ep, Pb=Pp, Pf=Pf, eb=0.0, ef=Ef)

        def fb(e):
            Pb = prm.Pr * (2.0 * (Php - self.Phi(e)) - 1.0)
            return self.e0eb * e - (Pb - self.D_th / prm.A_R)

        eb = self._bisect(fb, 0.0, Ep)
        Pb = prm.Pr * (2.0 * (Php - self.Phi(eb)) - 1.0)
        return dict(Ep=Ep, Pb=Pb, Pf=Pf, eb=eb, ef=Ef)

    def MW(self):
        r = self.retained()
        return self.prm.t_FE * (r["Pb"] - r["Pf"]) / self.e0eb

    def MW_max(self):
        """Fully programmed limit: E_p -> infinity."""
        prm = self.prm

        def fb(e):
            Pb = prm.Pr * (2.0 * (1.0 - self.Phi(e)) - 1.0)
            return self.e0eb * e - (Pb - self.D_th / prm.A_R)

        eb = self._bisect(fb, 0.0, 3.0e6)
        Pb = prm.Pr * (2.0 * (1.0 - self.Phi(eb)) - 1.0)

        def ff(e):
            Pf = prm.Pr * (2.0 * self.Phi(e) - 1.0)
            return self.e0eb * e - (self.D_th / prm.A_R - Pf)

        ef = self._bisect(ff, 0.0, 3.0e6)
        Pf = prm.Pr * (2.0 * self.Phi(ef) - 1.0)
        return prm.t_FE * (Pb - Pf) / self.e0eb


def a_r_sat_closed_form(prm: Params, gamma=0.98):
    """Closed-form estimate of the saturation area ratio.

    The window reaches gamma of its fully-programmed limit when the programming field E_s makes the
    switched fraction Phi(E_s) high enough; E_s is found from the relaxation condition, and the
    programming balance then gives A_R directly (no root finding in A_R):
        A_R_sat = (xp - t_FE E_s - psi_p) / [ c (eps0 eps_b E_s + Pr (2 Phi(E_s) - 1)) ]
    where c = 1/C_i + 1/C_FG and psi_p is the channel surface potential at the programmed displacement
    (evaluated self-consistently, psi_p ~ 0.1 V for the stack studied).
    """
    aw = AnalyticWindow(prm)
    mw_max = aw.MW_max()
    # find E_s such that the window equals gamma * MW_max (use a generic A_R-independent relaxation model)
    prm1 = replace(prm, A_R=1.0)
    e0eb = EPS0 * prm.eps_b

    def mw_of_Ep(Ep):
        Php = aw.Phi(Ep)

        def fb(e):
            Pb = prm.Pr * (2.0 * (Php - aw.Phi(e)) - 1.0)
            return e0eb * e - (Pb - aw.D_th / prm.A_R)

        eb = aw._bisect(fb, 0.0, 3.0e6)
        Pb = prm.Pr * (2.0 * (Php - aw.Phi(eb)) - 1.0)

        def ff(e):
            Pf = prm.Pr * (2.0 * aw.Phi(e) - 1.0)
            return e0eb * e - (aw.D_th / prm.A_R - Pf)

        ef = aw._bisect(ff, 0.0, 3.0e6)
        Pf = prm.Pr * (2.0 * aw.Phi(ef) - 1.0)
        return prm.t_FE * (Pb - Pf) / e0eb

    Es = aw._bisect(lambda E: mw_of_Ep(E) - gamma * mw_max, 0.3 * prm.Ec, 3.0e6)
    Pp = prm.Pr * (2.0 * aw.Phi(Es) - 1.0)
    Dfe = e0eb * Es + Pp
    # self-consistent psi_p: D_p = A_R * Dfe, with A_R unknown -> iterate
    psi = 0.1
    for _ in range(50):
        AR = (prm.xp - prm.t_FE * Es - psi) / (prm.c_series * Dfe)
        D = AR * Dfe
        psi = aw.ch.u_of_p(max(D / Q, 1e-30))
    return dict(A_R_sat=AR, E_s=Es, psi_p=psi, MW_max=mw_max)


if __name__ == "__main__":
    import time

    t0 = time.time()
    prm = Params()
    r = window(prm, detail=True)
    print("baseline", r, "t=%.2fs" % (time.time() - t0))
    print("analytic", AnalyticWindow(prm).MW(), AnalyticWindow(prm).MW_max())
