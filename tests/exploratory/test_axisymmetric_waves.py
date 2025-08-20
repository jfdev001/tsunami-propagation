from numpy import meshgrid, cos, exp
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import jv
from unittest import TestCase, main


def integrand_strip_source_axisymmetric_wave(rho, r, t):
    """Integrand in Equation (7) from Carrier2005."""
    order = 0
    return rho*jv(order, rho*r)*cos(rho*t*exp((-rho**2)/4))


class TestAxisymmetricWaves(TestCase):
    pass


if __name__ == "__main__":
    main()
