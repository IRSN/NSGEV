README
================
Yves Deville
05/09/2016

Welcome to the **NSGEV** package! The main function of the package is
`TVGEV` which creates an object with class `"TVGEV"` representing a
Time-Varying model with GEV margins. This kind of model is especially
useful to study *block maxima*, usually annual maxima. A popular use is
assessing the impact of global warming using series of annual maxima of
the daily maximal temperatures.

# NEWS

See the file [NEWS.md](NEWS.md). This file will now on be used in place
of the `ChangeLog` file.

The **nieve** package used to provide the GEV distribution functions
also has its own NEWS file.

# INSTALLATION

## Using the *remotes* package

In an R session use

``` r
library(remotes)
install_github("IRSN/NSGEV", dependencies = TRUE)
```

This should install the package and make it ready to use.

Mind that by default this does not build the vignette shipped with the
package (long-form documentation). To build the vignette, use instead

``` r
install_github("IRSN/NSGEV", dependencies = TRUE, build_vignettes = TRUE)
```

The installation will then take a longer time but the vignette will be
accessible from the help of the package (link above the “Help Pages”
section).

You can also select a specific branch or a specific commit by using the
suitable syntax for `install_github`. For instance to install the branch
`develop` use

``` r
install_github("IRSN/NSGEV@develop", dependencies = TRUE)
```

See the **remotes** package documentation for more details.

## Clone, build and install

### Cloning the repository

You can also clone the repository to install the package. If you do not
have yet a local `NSGEV` repository, use `git clone` to clone the
`NSGEV` repository

``` bash
git clone https://github.com/IRSN/NSGEV
```

This will create a `NSGEV` sub-directory of the current directory,
i.e. the directory from which the git command was issued. Of course this
can work only if you have the authorisation to clone.

### Installation

Move to the parent directory of your cloned repository and use the
following command from a terminal to create a tarball source file

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

You can also use the **RStudio** IDE to install the package.
