# CaTchDes
MATLAB codes for Caratheodory-Tchakaloff Near-Optimal Regression Designs 

The software package CaTchDes contains two main Matlab functions for the computation of near-optimal sampling sets and weights 
(designs) for polynomial regression on discrete design spaces (for example grid discretizations of planar, surface and solid 
domains). This topic has strong connections with computational statistics and approximation theory.

The main functions are:

NORD: computes a Near Optimal Regression Design with a given G-efficiency on a discrete set X in 2d or 3d, by a basic 
multiplicative algorithm and Caratheodory-Tchakaloff discrete measure concentration

CTDC: computes Caratheodory-Tchakaloff discrete measure concentration in 2d and 3d, for example probability measures (designs) 
but also quadrature formulas; the moments are invariant (close to machine precision) up to degree deg; adapts to the 
(numerical) dimension of the polynomial space on the support; works satisfactorily for low/moderate degrees


Auxiliary function: 

chebvand2d3d: computes by recurrence the Chebyshev-Vandermonde matrix on a 2d or 3d arbitrarily located mesh, in the 
total-degree product Chebyshev basis of a given rectangle with the graded lexicographical order
