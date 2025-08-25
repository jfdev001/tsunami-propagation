from numpy import (linspace, cos, exp, abs, pi, sqrt, sin, zeros_like, inf,
                   atleast_1d)
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import j0, ellipk, ellipe
from unittest import TestCase, main, skip

from pdb import set_trace


def integrand_superposition_simple_axisymmetric_wave(rho, r, t):
    """Integrand in Equation (7) from Carrier2005."""
    return rho*j0(rho*r)*cos(rho*t)*exp((-rho**2)/4)


def integrand_superposition_with_Gdot_axisymmetric_wave(rho, r, t):
    """Integrand in Equation (8) from Carrier2005."""
    Gdot = piecewise_Gdot
    return 2*exp(-rho**2)*Gdot(rho, r, t)


def piecewise_Gdot(rho, r, t):
    rho = atleast_1d(rho)
    Gdot = zeros_like(rho, dtype=float)

    # Case 1: t > r + rho
    mask1 = t > (r + rho)
    rho1 = rho[mask1]
    Gdot[mask1] = Gdot1(rho1, r, t)

    # Case 2: |r - rho| < t < r + rho
    mask2 = (abs(r - rho) < t) & (t < (r + rho))
    rho2 = rho[mask2]
    Gdot[mask2] = Gdot2(rho2, r, t)

    # Case 3: Gdot stays zero
    return Gdot if Gdot.size > 1 else Gdot.item()


def Gdot1(rho, r, t):
    """Equation (10) derivative with respect to `t` of the first piece.

    Case 1: t > r + rho
    """
    return -2*rho*t*(-(-4*r*rho/(t**2 - (r - rho)**2) + 1)*ellipk(4*r*rho/(t**2 - (r - rho)**2)) + ellipe(4*r*rho/(t**2 - (r - rho)**2)))/(pi*(t**2 - (r - rho)**2)**(3/2)*(-4*r*rho/(t**2 - (r - rho)**2) + 1)) - 2*rho*t*ellipk(4*r*rho/(t**2 - (r - rho)**2))/(pi*(t**2 - (r - rho)**2)**(3/2))


def Gdot2(rho, r, t):
    """Equation (10) second piece.

    Case 2: |r - rho| < t < r + rho
    """
    return t*sqrt(rho/r)*(-(1 - (t**2 - (r - rho)**2)/(4*r*rho))*ellipk((t**2 - (r - rho)**2)/(4*r*rho)) + ellipe((t**2 - (r - rho)**2)/(4*r*rho)))/(pi*(1 - (t**2 - (r - rho)**2)/(4*r*rho))*(t**2 - (r - rho)**2))


def integrand_superposition_no_Gdot_axisymmetric_wave(rho, r, t):
    """Integrand in Equation (9) from Carrier2005."""
    G = piecewise_G
    integrand = 2*exp(-rho**2)*G(rho, r, t)
    return integrand


def piecewise_G(rho, r, t):
    """Equation (10) from Carrier (2005)"""
    rho = atleast_1d(rho)
    G = zeros_like(rho, dtype=float)

    # define k parameter
    denom = t**2 - (r - rho)**2
    k = (4*r*rho) / denom

    # Case 1: t > r + rho
    mask1 = t > (r + rho)
    G[mask1] = (2*rho[mask1])/(pi*sqrt(denom[mask1])) * ellipk(k[mask1])

    # Case 2: |r - rho| < t < r + rho
    mask2 = (abs(r - rho) < t) & (t < (r + rho))
    G[mask2] = (1/pi)*sqrt(rho[mask2]/r) * ellipk(1/k[mask2])

    # Case 3: G stays zero
    return G if G.size > 1 else G.item()


def integrate_rho_for_ts_and_rs(
        integrand, rho_start, rho_stop, ts, rs):
    """ 
    TODO: should probably do this with symbolic library?
    i.e., once integral is computed, then use d/dt ??...
    could do something like integrate Gdot1 to inf and
    Gdot2 to inf, get these closed forms, then use these
    close forms in branches based on the cases defining G
    TODO: Other approaches: analyze G_t to see
    where discontinuities would appear and integrate around
    those... or analyze plots of integrand no gdot
    and visually identify discontinuities...
    """
    t_to_r_to_wave_displacement = []
    for t in ts:
        r_to_wave_displacement = []
        for r in rs:
            wave_displacement, err = quad(
                integrand,
                rho_start,
                rho_stop,
                args=(r, t)
            )
            r_to_wave_displacement.append(wave_displacement)
        t_to_r_to_wave_displacement.append(r_to_wave_displacement)
    return t_to_r_to_wave_displacement


class TestAxisymmetricWaves(TestCase):
    def test_plot_fig2_carrier2002(self):
        """Reproduce figure 2 from Carrier2002, Green's function and discontin.

        G(b, sigma, lambda) == G(rho, r, t)
        """
        ts = lambdas = [0.1, 0.2, 0.3, 0.4, 0.5]
        r = sigma = 0.05
        rhos = bs = linspace(0, 1, 100)

        def singularity_at_rho(t, r):
            # not including 1/2 because t = 1/2 lambda it seems...
            return t - r

        G = piecewise_G
        Gs = [G(rhos, r, t) for t in ts]
        fig, ax = plt.subplots()
        for ix, t in enumerate(ts):
            singularity = singularity_at_rho(t, r)
            ax.vlines(
                singularity, ymin=0, ymax=max(Gs[ix]), color="red",
                label="singularity" if ix == 0 else "")
            ax.plot(rhos, Gs[ix], label=f"t={t}")
        ax.set_xlabel("rho")
        ax.set_ylabel("G")
        ax.legend()
        plt.show()
        return

    def test_plot_fig2(self):
        """Reproduce figure (2) from Carrier 2005."""
        t = 102
        r = 100
        fig, axs = plt.subplots(3, 1, figsize=(6, 10))

        # 2af
        rho = linspace(0, 5, 840)
        fig2a_out = integrand_superposition_simple_axisymmetric_wave(
            rho, r, t)
        axs[0].plot(rho, fig2a_out, alpha=0.75)

        # 2b
        rho = linspace(0, 3, 500)
        fig2b_out = integrand_superposition_with_Gdot_axisymmetric_wave(
            rho, r, t)
        axs[1].plot(rho, fig2b_out)
        axs[1].set_ylim(-0.06, 0.06)

        # 2c
        fig2c_out = integrand_superposition_no_Gdot_axisymmetric_wave(
            rho, r, t)
        axs[2].plot(rho, fig2c_out)

        axs[2].set_xlabel("rho")
        # plt.show()
        return

    @skip("not correct yet")
    def test_plot_fig3(self):
        """Numerically solve equation (9).

        TODO: Then take derivative??

        TODO: The Green function has a known singularity!! 
        Carrier2002 "Tsunami Run-up and Draw-Down on a Plane Beach"....
        could either (a) split integral to perform integration? or (b)
        try and see about ISML software? since this what's recommneded in the
        paper...
        """
        rs = linspace(0, 10, 100)
        t_a = [0.5, 0.75, 1.0, 1.5, 2.0, 5.0]
        fig, axs = plt.subplots(2, 1, figsize=(15, 8))

        # plot fig3a
        t_a_to_r_to_wave_displacement = integrate_rho_for_ts_and_rs(
            integrand_superposition_no_Gdot_axisymmetric_wave,
            rho_start=0, rho_stop=3, ts=t_a, rs=rs)

        for tix in range(len(t_a)):
            t = t_a[tix]
            r_to_wave_displacement = t_a_to_r_to_wave_displacement[tix]
            axs[0].plot(rs, r_to_wave_displacement, label=f"t={t}")

        axs[0].legend()

        t_b = [1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0]
        t_b_to_r_to_wave_displacement = integrate_rho_for_ts_and_rs(
            integrand_superposition_no_Gdot_axisymmetric_wave,
            rho_start=0, rho_stop=3, ts=t_b, rs=rs)

        for tix in range(len(t_b)):
            t = t_b[tix]
            r_to_wave_displacement = t_b_to_r_to_wave_displacement[tix]
            axs[1].plot(rs, r_to_wave_displacement, label=f"t={t}")

        axs[1].legend()

        plt.show()
        return


if __name__ == "__main__":
    main()
