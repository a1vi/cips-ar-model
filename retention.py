"""Retention-aware extension of the analytic window (thermally activated, Merz-type back-switching).

A hysteron with coercive field Ec_i held under a reverse field e switches after
    tau_i(e) = tau0 * exp[ b * ((Ec_i/e)^n - 1) ],
so that tau_i = tau0 at e = Ec_i (the quasi-static model of the main text corresponds to t = tau0).
After a hold time t the hysterons with Ec_i <= g e have switched, with
    g(t) = [1 + ln(t/tau0)/b]^(1/n)  >= 1.
Replacing Phi(e) by Phi(g e) in the relaxation balances gives the window at time t.
"""
import math
from cips_ar import Params, AnalyticWindow


class RetentionWindow(AnalyticWindow):
    def __init__(self, prm, tau0=1e-6, b=20.0, n=1.0):
        super().__init__(prm)
        self.tau0, self.b, self.n = tau0, b, n

    def g(self, t):
        if t <= self.tau0:
            return 1.0
        return (1.0 + math.log(t / self.tau0) / self.b) ** (1.0 / self.n)

    def retained_t(self, t):
        prm, g = self.prm, self.g(t)
        Ep = self.program_field()
        Php = self.Phi(Ep)
        Pp = prm.Pr * (2.0 * Php - 1.0)

        def Pf_of(e):
            return prm.Pr * (2.0 * self.Phi(g * e) - 1.0)

        def ff(e):
            return self.e0eb * e - (self.D_th / prm.A_R - Pf_of(e))
        ef = self._bisect(ff, 0.0, 3.0e6)
        Pf = Pf_of(ef)
        if Pp <= self.D_th / prm.A_R:
            return dict(Ep=Ep, Pb=Pp, Pf=Pf, eb=0.0, ef=ef)

        def Pb_of(e):
            return prm.Pr * (2.0 * (Php - self.Phi(min(g * e, Ep))) - 1.0)

        def fb(e):
            return self.e0eb * e - (Pb_of(e) - self.D_th / prm.A_R)
        eb = self._bisect(fb, 0.0, Ep)
        return dict(Ep=Ep, Pb=Pb_of(eb), Pf=Pf, eb=eb, ef=ef)

    def MW_t(self, t):
        r = self.retained_t(t)
        return max(self.prm.t_FE * (r["Pb"] - r["Pf"]) / self.e0eb, 0.0)
