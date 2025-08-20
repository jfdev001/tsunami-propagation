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
