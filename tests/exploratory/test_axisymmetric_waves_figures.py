from argparse import ArgumentParser
import matplotlib.pyplot as plt
from mpmath import quad as mpquad
from mpmath import exp as mpexp
from numpy import (linspace, cos, exp, abs, pi, sqrt, sin, zeros_like, inf,
                   atleast_1d, isclose)
from numpy import all as np_all
from scipy.integrate import quad, IntegrationWarning
from scipy.special import j0, ellipk, ellipe
from unittest import TestCase, main, skip
from typing import Callable
from warnings import simplefilter

from pdb import set_trace

simplefilter("always", IntegrationWarning)


def integrand_superposition_simple_axisymmetric_wave(rho, r, t):
    """Integrand in Equation (7) from Carrier2005."""
    return rho*j0(rho*r)*cos(rho*t)*exp((-rho**2)/4)


# TODO: could have adaptive delta... i.e., proportion of singularity or try
# until no integration warning???
def split_integrate_with_Gdot_axisymmetric_wave(
        rho_start: float, rho_stop: float,
        r: float, t: float,
        use_manual_splitting: bool,
        singularity_delta: float,
        **quad_kwargs):
    rho_singularity = singularity_at_rho(t, r)
    # use quad builtin in splitting
    if not use_manual_splitting:
        integral, err = quad(
            integrand_superposition_with_Gdot_axisymmetric_wave,
            rho_start,
            rho_stop,
            args=(r, t),
            points=[rho_singularity],
            **quad_kwargs)
    # manually split integral
    else:
        integral_one, err_one = quad(
            integrand_superposition_with_Gdot_axisymmetric_wave,
            rho_start,
            rho_singularity - singularity_delta,
            args=(r, t),
            **quad_kwargs)
        integral_two, err_two = quad(
            integrand_superposition_with_Gdot_axisymmetric_wave,
            rho_singularity + singularity_delta,
            rho_stop,
            args=(r, t),
            **quad_kwargs)
        integral = integral_one + integral_two
    return integral


def integrand_superposition_with_Gdot_axisymmetric_wave(rho, r, t):
    """Integrand in Equation (8) from Carrier2005."""
    Gdot = greens_function_dt
    return 2*exp(-rho**2)*Gdot(rho, r, t)


def greens_function_dt(rho, r, t):
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


def piecewise_mpmath_Gdot(rho, r, t):
    raise


# TODO: could have args for the callables to either np/scipy/mpath funcs
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
    G = greens_function
    integrand = 2*exp(-rho**2)*G(rho, r, t)
    return integrand


def greens_function(rho, r, t):
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
                args=(r, t),
                limit=1000,
            )
            r_to_wave_displacement.append(wave_displacement)
        t_to_r_to_wave_displacement.append(r_to_wave_displacement)
    return t_to_r_to_wave_displacement


def compute_integral_of_rho_for_ts_and_rs(
        compute_integral, rho_start, rho_stop, ts, rs):
    """
    Dirty copy with one change
    """
    t_to_r_to_wave_displacement = []
    for t in ts:
        r_to_wave_displacement = []
        for r in rs:
            integral = compute_integral(rho_start, rho_stop, t, r)
            r_to_wave_displacement.append(integral)
        t_to_r_to_wave_displacement.append(r_to_wave_displacement)
    return t_to_r_to_wave_displacement


def plot_singularities_in_greens_function(ax, rhos, r, ts):
    G = greens_function
    Gs = [G(rhos, r, t) for t in ts]
    for ix, t in enumerate(ts):
        singularity = singularity_at_rho(t, r)
        ax.vlines(
            singularity, ymin=0, ymax=max(Gs[ix]), color="red",
            label="singularity" if ix == 0 else "")
        ax.plot(rhos, Gs[ix], label=f"t={t}")
    return


def singularity_at_rho(t, r):
    # not including 1/2 because t = 1/2 lambda it seems...
    return t - r


def plot_sanity_checked_integrand(
        compute_integrand: Callable,
        suptitle: str,
        include_singularities: bool,
        rho_start: float = 0,
        rho_stop: float = 3,
        n_rho_points: int = 500,
        rs: list[float] = [0, 2, 4],
        ts: list[float] = [0.5, 0.75, 1.0],
        ylims: tuple[float] = tuple(),
        **split_integrate_kwargs) -> tuple:
    """
    args:
        compute_integrand: Func with params for rho, r, t
    """
    rhos = linspace(rho_start, rho_stop, n_rho_points)
    if len(split_integrate_kwargs) == 0:
        split_integrate_kwargs = dict(
            use_manual_splitting=True,
            singularity_delta=1E-5)

    fig, axs = plt.subplots(2, 1, figsize=(6, 10))
    singularity_labeled_in_legend = False
    for r in rs:
        for t in ts:
            integrand = compute_integrand(rhos, r, t)
            if not np_all(isclose(integrand, 0)):
                singularity = singularity_at_rho(t, r)
                if singularity >= 0:
                    if include_singularities:
                        axs[0].vlines(
                            singularity,
                            min(integrand), max(integrand),
                            color="red",
                            label="singularity"
                            if not singularity_labeled_in_legend else "")
                        singularity_labeled_in_legend = True

                    # Compute and plot split integral
                    integral = split_integrate_with_Gdot_axisymmetric_wave(
                        rho_start, rho_stop, r, t,
                        **split_integrate_kwargs)
                else:
                    integral, err = quad(
                        integrand_superposition_with_Gdot_axisymmetric_wave,
                        rho_start, rho_stop, args=(r, t))

                round_digits_to_n_places = 2
                integral = round(integral, round_digits_to_n_places)

                # plot non zero entire integrands
                axs[0].plot(
                    rhos,
                    integrand, label=f"t={t}, r={r}, I={integral}")

                # plot only positive part
                positive_integrand_mask = integrand > 0
                positive_integrand = integrand[positive_integrand_mask]
                parallel_positive_rhos = rhos[positive_integrand_mask]
                if positive_integrand.shape[0] > 0:
                    axs[1].plot(
                        parallel_positive_rhos,
                        positive_integrand,
                        label=f"t={t}, r={r}, I={integral}")

    axs[0].set_title("entire domain of integrand")
    axs[0].set_ylim(*ylims)
    axs[0].legend()

    axs[1].set_title("positive parts of integrands")
    axs[1].legend()

    fig.supxlabel("rho")
    fig.suptitle(suptitle)
    return fig, axs


def plot_hlines_through_origin(axs) -> None:
    for ax in axs:
        ax.axhline(0, color="black")
    return


class TestAxisymmetricWaves(TestCase):
    def test_plot_fig2_greens_func_and_discontinuities(self):
        """Reproduce figure 2 from Carrier2002, Green's function and discontin.

        G(b, sigma, lambda) == G(rho, r, t)
        """

        fig, axs = plt.subplots(2, 1, figsize=(6, 10))
        fig.suptitle("Singularities of Green's Function")

        G = greens_function

        # From carrier2002
        ts = lambdas = [0.1, 0.2, 0.3, 0.4, 0.5]
        r = sigma = 0.05
        rhos = bs = linspace(0, 1, 100)
        plot_singularities_in_greens_function(axs[0], rhos, r, ts)
        axs[0].set_ylabel("G")
        axs[0].set_title(f"r = {r}")
        axs[0].legend()

        # Analysis for fig2a of carrier2005
        ts = [102]
        r = 100
        rhos = linspace(0, 3, 500)
        plot_singularities_in_greens_function(axs[1], rhos, r, ts)
        axs[1].set_xlabel("rho")
        axs[1].set_ylabel("G")
        axs[1].set_title(f"r = {r}")
        axs[1].legend()

        fig.suptitle("Supplement (2): Green's Function and its Singularities")

        return

    def test_plot_different_integrands_at_fixed_t_and_r(self):
        """Reproduce figure (2) from Carrier 2005."""
        t = 102
        r = 100
        fig, axs = plt.subplots(3, 1, figsize=(6, 10))

        # 2a
        rho = linspace(0, 5, 840)
        fig2a_out = integrand_superposition_simple_axisymmetric_wave(
            rho, r, t)
        axs[0].plot(rho, fig2a_out, alpha=0.75)
        axs[0].set_title("Equation (7)")

        # 2b
        rho = linspace(0, 3, 500)
        fig2b_out = integrand_superposition_with_Gdot_axisymmetric_wave(
            rho, r, t)
        axs[1].plot(rho, fig2b_out)
        axs[1].set_ylim(-0.06, 0.06)
        axs[1].set_title("Equation (8)")

        # 2c
        fig2c_out = integrand_superposition_no_Gdot_axisymmetric_wave(
            rho, r, t)
        axs[2].plot(rho, fig2c_out)

        axs[2].set_xlabel("rho")
        axs[2].set_title("Equation (9)")

        fig.suptitle("Figure (2): Different Integrands for t=102, r=100")

        plot_hlines_through_origin(axs)
        return

    def test_plot_sanity_check_bessel_integrand_and_computed_integrals(self):
        suptitle = (
            "Supplement (3.2):\n"
            r"$\text{integrand}(\rho, r, t) = \rho\ J_0(\rho r)\ \cos(\rho t)\ \exp(\frac{-\rho^2}{4})$"
        )
        include_singularities = False
        compute_integrand = integrand_superposition_simple_axisymmetric_wave
        fig, axs = plot_sanity_checked_integrand(
            compute_integrand, suptitle, include_singularities)
        plot_hlines_through_origin(axs)
        return

    def test_plot_sanity_check_Gdot_integrand_and_computed_integrals(self):
        """
        Plot the integrand with Gdot if values are not 0 for ts and rs.
        Mark the computed singularity with a vertical. Then compute the the
        integral to sanity check values.
        """

        suptitle = (
            "Supplement (3.1):\n"
            r"$\text{integrand}(\rho, r, t) = 2 \exp(-\rho^2) G_t(\rho, r, t)$"
            " with Singularities"
        )
        include_singularities = True
        ylims = (-100, 100)
        compute_integrand = integrand_superposition_with_Gdot_axisymmetric_wave
        fig, axs = plot_sanity_checked_integrand(compute_integrand, suptitle,
                                                 include_singularities, ylims=ylims)
        plot_hlines_through_origin(axs)
        return

    def test_plot_integral_of_axisymmetric_wave_integrand(self):
        """Reproduce figure 3 from Carrier 2005 

        TODO: The Green function has a known singularity!!
        Carrier2002 "Tsunami Run-up and Draw-Down on a Plane Beach"....
        could either (a) split integral to perform integration? or (b)
        try and see about ISML software? since this what's recommneded in the
        paper...

        TODO: maybe move this also into another function so that you can
        plug and play with different integrands...
        """
        integrand = integrand_superposition_simple_axisymmetric_wave

        rs = linspace(0, 10, 100)
        t_a = [0.5, 0.75, 1.0, 1.5, 2.0, 5.0]
        fig, axs = plt.subplots(2, 1, figsize=(6, 8))
        rho_start = 0
        rho_stop = 3

        # plot fig3a
        # t_a_to_r_to_wave_displacement = compute_integral_of_rho_for_ts_and_rs(
        # split_integrate_with_Gdot_axisymmetric_wave,
        # rho_start,
        # rho_stop,
        # t_a,
        # rs)
        t_a_to_r_to_wave_displacement = integrate_rho_for_ts_and_rs(
            integrand, rho_start=rho_start, rho_stop=rho_stop, ts=t_a, rs=rs)
        for tix in range(len(t_a)):
            t = t_a[tix]
            r_to_wave_displacement = t_a_to_r_to_wave_displacement[tix]
            axs[0].plot(rs, r_to_wave_displacement, label=f"t={t}")

        axs[0].legend()

        # plot fig3b
        t_b = [1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0]
        t_b_to_r_to_wave_displacement = integrate_rho_for_ts_and_rs(
            integrand,
            rho_start=rho_start, rho_stop=rho_stop, ts=t_b, rs=rs)

        for tix in range(len(t_b)):
            t = t_b[tix]
            r_to_wave_displacement = t_b_to_r_to_wave_displacement[tix]
            axs[1].plot(rs, r_to_wave_displacement, label=f"t={t}")

        axs[1].legend()
        axs[1].set_xlabel("r")

        title = r"$\int_0^3 2\exp(-\rho^2)"
        if "with_Gdot" in integrand.__name__:
            title += r"G_t(\rho, r, t)$"
        elif "no_Gdot" in integrand.__name__:
            title += r"G(\rho, r, t)$"
        elif "simple" in integrand.__name__:
            title = r"$\int_0^3 \rho\ J_0(\rho r)\ \cos(\rho t)\ \exp(\frac{-\rho^2}{4})$"
        else:
            raise ValueError(f"unrecognized integrand: {integrand.__name__}")
        fig.suptitle(f"Figure (3): {title}")

        plot_hlines_through_origin(axs)
        return

    @classmethod
    def tearDownClass(cls):
        plt.show()
        return


if __name__ == "__main__":
    # parser = ArgumentParser("additional args for testing")
    # parser.add_argument("flag", required=True)
    # args = parser.parse_args()
    main()
