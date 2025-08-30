import matplotlib.pyplot as plt
from numpy import linspace, trapezoid
from pdb import set_trace
from scipy.integrate import quad, IntegrationWarning
from unittest import TestCase, main
from warnings import simplefilter

from tests.exploratory.test_axisymmetric_waves_figures import (
    integrand_superposition_with_Gdot_axisymmetric_wave,
    singularity_at_rho)

simplefilter("always", IntegrationWarning)


class TestIntegrationWarning(TestCase):
    def test_integrate_Gdot(self):
        """
        Adaptive integration shouldn't really be necessary because 
        you can see from the function that the integral shape is always the 
        same... should likely be sufficient just to do split interval integration
        using singularity - delta where delta is 1E-5 since you already know
        from other tests that this reproduces the expected integral

        References:
            * Should be able to make integration points relatively imprecise, 
            i.e., singularity +- some numbers
            https://stackoverflow.com/questions/64992360/how-does-points-work-in-integrate-quad
        """
        ts = [0.1, 0.2, 0.3, 0.4, 0.5]
        r = 0.05
        rho_start = 0
        rho_stop = 1
        rhos = linspace(rho_start, rho_stop, 500)
        fig, ax = plt.subplots()
        delta = 0.5
        for t in ts:
            singularity = singularity_at_rho(t, r)
            singularity_left = singularity - singularity*delta
            singularity_right = singularity + singularity*delta
            points = [singularity_left,
                      singularity_right] if singularity > 0 else None
            print(points)
            integral, err, info, message = quad(
                integrand_superposition_with_Gdot_axisymmetric_wave,
                rho_start,
                rho_stop,
                args=(r, t),
                points=points,
                limit=500,
                full_output=True)
            print(t, message)
            Gs = integrand_superposition_with_Gdot_axisymmetric_wave(
                rhos, r, t)
            ax.plot(rhos, Gs, label=f"t={t}, I={integral:.2f}")
        ax.legend()
        plt.show()
        return
