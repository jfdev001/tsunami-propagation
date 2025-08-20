% Reproduce figure 2 from Carrier 2005 

r = 100
t = 102
rho = linspace(0, 5, 100)
n = axisymmetric_wave(r, t, rho)

function n = axisymmetric_wave(r, t, rho)
    order = 0;
    n = rho*besselj(order, rho*r)*cos(rho*t*exp((-rho*rho)/4));
end

