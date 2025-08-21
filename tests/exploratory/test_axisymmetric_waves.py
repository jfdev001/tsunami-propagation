from numpy import linspace, cos, exp, abs, pi, sqrt, sin
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import j0
from unittest import TestCase, main


def integrand_superposition_simple_axisymmetric_wave(rho, r, t):
    """Integrand in Equation (7) from Carrier2005."""
    return rho*j0(rho*r)*cos(rho*t)*exp((-rho**2)/4)


def integrand_superposition_with_dt_axisymmetric_wave(rho, r, t):
    """Integrand in Equation (8) from Carrier2005."""
    raise NotImplementedError("requires d/dt[G]")


def integrand_superposition_no_dt_axisymmetric_wave(rho, r, t):
    """Integrand in Equation (9) from Carrier2005."""
    G = piecewise_G
    return 2*exp(-rho**2)*G(rho, r, t)


def piecewise_G(rho, r, t):
    """Equation (10) from Carrier2005"""
    G = None
    K = complete_elliptic_integral_of_first_kind
    k = (4*r*rho)/(t**2 - (r - rho)**2)

    raise NotImplementedError(
        "ambiguity with true values,, array vs single val")
    if t > (r + rho):
        G = (2*rho)/(pi*sqrt(t**2 - (r - rho)**2))*K(k)
    elif abs(r - rho) < t and t < (r + rho):
        G = (1/pi)*(sqrt(rho/r))*K(1/k)
    elif t < abs(r - rho):
        G = 0
    else:
        raise ValueError("unexpected rho, r, t args")

    return G


def complete_elliptic_integral_of_first_kind(k):
    """Defined after equation (10) in Carrier2005"""
    def integrand(v, k):
        return (1/sqrt(1 - k*(sin(v))**2))

    return quad(integrand, 0, pi/2, k)


class TestAxisymmetricWaves(TestCase):
    def test_plot_fig2(self):
        """Reproduce figure (2) from Carrier 2005."""
        rho = linspace(0, 5, 840)
        t = 102
        r = 100
        fig, axs = plt.subplots(3, 1, figsize=(8, 15))

        fig2a_out = integrand_superposition_simple_axisymmetric_wave(
            rho, r, t)
        axs[0].plot(rho, fig2a_out, alpha=0.75)

        pass

        # fig2c_out = integrand_superposition_no_dt_axisymmetric_wave(
        # rho, r, t)
        # axs[2].plot(rho, fig2c_out)

        axs[2].set_xlabel("rho")
        plt.show()
        return


if __name__ == "__main__":
    main()
