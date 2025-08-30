# tsunami-propagation

A Python and PETSc implementation of the shallow water equations as described 
in [Carrier 2005](https://www.techscience.com/CMES/v10n2/24866). The [method
of manufactured solutions](https://mooseframework.inl.gov/python/mms.html) could
be used to verify a numerical solution, however an exact solution to the shallow water equations 
is provided in Carrier 2005 so that will be used to verify a backward time central space
discretization implementation in PETSc---i.e., an implicit time scheme with a 
central finite difference spatial discretization.

The governing equation of interest is given by

$$
g h \Delta \eta - \eta_{tt}= 0,
$$

where $\Delta$ is the Laplace operator in dimensional coordinates, $h$ is the
depth, $t$ is the time, $g$ is the gravitational acceleration, and the letter
subscript idenotes partial differentiation.

For figures 2 and 3, the below equations are relevant:

Exact integral representation of wave displacement,

$$
\eta(r, t) = \int_0^\infty \rho J_0(\rho r)\ \cos(\rho t)\ \exp(\frac{-\rho^2}{4})\ d\rho.  \tag{eq. 7}
$$

The above equation using the Bessel function of the first kind of order 0 can be
rewritten as below:

Wave displacement *with* time derivative in the integrand,

$$
\eta(r, t) = \int_{0}^{\infty} 2 \exp(-\rho^2) G_t(\rho, r, t)\ d\rho, \tag{eq. 8}
$$

where 

$$
G(\rho, r, t) = \begin{cases}
    \frac{2\rho}{\pi \sqrt{t^2 - (r-\rho)^2}} K(\frac {4 r \rho} {t^2 - (r - \rho)^2}) ,\ \text{for}\ t > r + \rho, \\
    \frac{1}{\pi} \sqrt{\frac{\rho}{r}}K(\frac{t^2 - (r - \rho)^2}{4 r \rho}),\ \text{for}\ |r - \rho| < t < r + \rho, \\
    0,\ \text{for}\ t < |r - \rho|,
\end{cases}
$$

and $K$ is the Complete Elliptical Integral of the first kind, 

$$
K(k) = \int_{0}^{\pi / 2} \frac{dv}{\sqrt{1 - k \sin^2 v}}.
$$

The wave displacement *without* time derivative in the integrand, 

$$
\eta(r, t) = \frac{\partial}{\partial t}\int_{0}^{\infty} 2 \exp(-\rho^2) G(\rho, r, t)\ d\rho. \tag{eq. 9}
$$

# References


