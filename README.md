README
================
Yves Deville
05/09/2016

Welcome to the **NSGEV** package! The main function of the package is
`TVGEV` which creates an object with class `"TVGEV"` representing a
Time-Varying model with GEV margins. This kind of model is especially
useful to study block maxima, usually annual maxima. A popular use is
assessing the impact of global warming using series of annual maxima of
the daily maximal temperatures.

# INSTALLATION

This package needs compilation. Hence if you are using a MS Windows
system you need to have the Rtools installed in order to install
**NSGEV** from its sources.

# NEWS

## 2022-05

The probability functions (density, distribution and quantile) are now
computed differently for small shape *ξ* ≈ 0, using a Taylor
approximation at *ξ* = 0. These approximations are used when
\|*ξ*\| \< *ϵ* where *ϵ* is can be chosen before compiling the package.
We use a 2-nd order approximation for the function, a 1-st order
approximation for their gradient and a zero-order (constant)
approximation for the Hessian. This leads to very good approximations
that are smooth and consistent.

The computations required are quite tedious. We used
[Maxima](https://maxima.sourceforge.io/) to compute the derivatives and
their Taylor approximation, along with the LaTeX package
[maxiplot](https://maxima.sourceforge.io/contrib/maxiplot/maxiplot_en.pdf).
However the formulas given by Maxima need some work to make them usable
in the numeric computations (in C). We need to express the derivatives
with a small number of well-chosen auxiliary variables that are computed
only once.

This feature seems to be unique among the available R packages devoted
to extreme-value. Although very small values of *ξ* are not frequently
used, it turns out that the exact derivatives are difficult to evaluate
for small *ξ* which may have some consequences in optimisations. This is
especially the case for the 2-nd order derivatives and for the quantile
function.
