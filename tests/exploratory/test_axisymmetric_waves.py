from numpy import linspace, cos, exp, abs, pi, sqrt, sin, zeros_like
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import j0, ellipk
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
    """Equation (10) from Carrier (2005)"""
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

        pass  # 2b

        rho = linspace(0, 3, 500)
        fig2c_out = integrand_superposition_no_dt_axisymmetric_wave(
            rho, r, t)
        axs[2].plot(rho, fig2c_out)

        axs[2].set_xlabel("rho")
        plt.show()
        return


if __name__ == "__main__":
    main()
