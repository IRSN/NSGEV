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

## INSTALLATION

This package needs compilation. Hence if you are using a MS Windows
system you need to have the Rtools installed in order to install
**NSGEV** from its sources.

## Using the *devtools* package

Note that if you are using Windows, you need to have the
[Rtools](https://cran.r-project.org/bin/windows/Rtools) installed.
Provided that the **devtools** package is installed you can then in an R
session use

``` r
library(devtools)
install_github("yvesdeville/NSGEV", dependencies = TRUE, auth_token = myToken)
```

where `myToken` stands for *your* token. This should install the package
and make it ready to use.

You can also select a specific branch or a specific commit by using the
suitable syntax for `install_github`, see the **devtools** package
documentation.

## Clone, build and install

### Cloning the repository

If you do not have yet a local `NSGEV` repository, use `git clone` to
clone the `NSGEV` repository

``` bash
git clone https://github.com/yvesdeville/NSGEV
```

This will create a `NSGEV` sub-directory of the current directory,
i.e. the directory from which the git command was issued. Of course this
can work only if you have the authorisation to clone.

### Installation on Unix and MacOs systems

With these systems you can install a package from its source. Move to
the parent directory of your cloned repository and use the following
command from a terminal to create a tarball source file

``` bash
R CMD build NSGEV
```

This will produce a source tarball `NSGEV_x.y.z` where `x`, `y` and `z`
stand for the major, minor and patch version numbers. Then you can
install from a command line

``` bash
R CMD INSTALL NSGEV_x.y.z
```

Note that you must also have all the packages required by **NSGEV**
installed.

## NEWS

### 2022-05

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

This feature is quite rare among the available R packages devoted to
extreme-value, although the
[revdbayes](https://github.com/paulnorthrop/revdbayes) package by Paul
Northorp uses the same principle in R implementations.

Even if very small values of *ξ* are not frequently used, it turns out
that the exact derivatives are difficult to evaluate for small *ξ*,
which may have some bad consequences in optimisations. This is
especially the case for the 2-nd order derivatives or for the
log-density function.
