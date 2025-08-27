#!/usr/bin/env python
""" 

References:
* https://stackoverflow.com/questions/38841587/how-do-i-define-a-conditional-function-using-sympy
* https://docs.sympy.org/latest/modules/functions/elementary.html#sympy.functions.elementary.piecewise.Piecewise
"""
from unittest import main, TestCase, skip
from sympy import (
    symbols, log, Piecewise, init_printing, pprint, elliptic_k, pi, sqrt, Abs,
    integrate, diff, exp, simplify)
from sympy import oo as inf

init_printing(use_unicode=True)
rho, r, t, v = symbols("rho r t v")


def wave_integrand():

    f = (2*rho)/(pi*sqrt(t**2 - (r - rho)**2))
    k1 = (4*r*rho)/(t**2 - (r - rho)**2)
    K1 = elliptic_k(k1)
    G1 = 2*exp(-rho**2)*f*K1
    case1 = t > (r + rho)

    g = (1/pi)*sqrt(rho/r)
    k2 = 1/k1
    K2 = elliptic_k(k2)
    G2 = 2*exp(-rho**2)*g*K2
    case2 = (Abs(r - rho) < t) & (t < (r + rho))

    G3 = 0
    case3 = t < Abs(r - rho)

    p = Piecewise((G1, case1), (G2, case2), (G3, case3))
    return p


class TestPiecewise(TestCase):
    @skip
    def test_simple_piecewise(self):
        x, y = symbols("x y")
        f = x**2
        g = log(x)
        p = Piecewise((0, x < -1), (f, x <= 1), (g, True))
        pprint(p)
        return

    def test_piecewise_wave(self):
        print("---------------")
        print("wave integrand")
        print("---------------")
        n = wave_integrand()
        pprint(n)

        print("---------------")
        print("wave integral")
        print("---------------")
        # TODO: problem here since integration is taking way too long...
        return
        integral = n.integrate((rho, 0, 3))
        pprint(integral)

        return


if __name__ == "__main__":
    main()
