"""Use symbolic manipulation to get d/dt[G] in equation (8) of Carrier 2005"""
from unittest import main, TestCase
from sympy import (
    symbols, sin, sqrt, pi, diff, integrate, init_printing, pprint, elliptic_k
)

init_printing(use_unicode=True)


class TestAxisymmetricSymbolic(TestCase):
    def test_symbolic_G(self):

        rho, r, t, v = symbols("rho r t v")

        f = (2*rho)/(pi*sqrt(t**2 - (r - rho)**2))
        k1 = (4*r*rho)/(t**2 - (r - rho)**2)
        K1 = elliptic_k(k1)
        G1 = f*K1
        G1dot = diff(G1, t)

        g = (1/pi)*sqrt(rho/r)
        k2 = 1/k1
        K2 = elliptic_k(k2)
        G2 = g*K2
        G2dot = diff(G2, t)

        print("G1 = ")
        pprint(G1)
        print()
        print(G1)
        print()
        print("G1dot = ")
        pprint(G1dot)
        print()
        print(G1dot)
        print("-----------------------------------------------")
        print("G2 = ")
        pprint(G2)
        print()
        print(G2)
        print()
        print("G2dot =")
        pprint(G2dot)
        print()
        print(G2dot)

        return


if __name__ == "__main__":
    main()
