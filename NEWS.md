
**NSGEV** Package News
===========================

# New in version 0.1.9 (pending)

## Enhancements

- The `TVGEV` function did not work when some (non temporal)
  covariates add missing values. 

- The functions `qMax.TVGEV` and `pMax.TVGEV` now compute the gradient
  and the Hessian w.r.t. the parameters, opening the road to profile
  likelihood inference on the quantile of the maximum. New function
  `dMax.TVGEV`.

## Bug fixes

- In `quantMax.TVGEV` the check on `prob` was misleading.

- The `quantMax.TVGEV` method did not work for `TVGEV` object
  describing a stationary model.


# New in version 0.1.8

## Enhancements

- For a `TVGEV` object, the `negLogLikFun` function now optionally
  returns the Hessian.
  
- New method `quantMax` for the `"TVGEV"` class and subsequent method
  `autoplot` for the related `quantMax.TVGEV` objects. Using these
  methods is illustrated in the main vignette.

- New methods `quantMaxFun` and `cdfMaxFun` for the `"TVGEV"`
  class. The return functions (closures).
  
- An `autoplot` method has been added for some classes:`"TVGEV"` and
  more.

- The vignettes of the **NSGEVVal** package are now included in
  **NSGEV**.  The required technical functions from **NSGEVVal** are
  now included in **NSGEV** as "check" functions. The bibliography
  has been somewhat cleaned.

## Technical changes (non user-visible)
 
- The `NAMESPACE` file is now generated with **roxygen2** making the
  package easier to maintain.

# News in version 0.1.7

## Bug fixes

- The `TVGEV` function did not work when a constant location was
  specified along with a non-constant scale and/or shape. Thanks to
  Jesper Rydén.

# News in version 0.1.6

## Changes

- The GEV probability functions `dGEV`, `pGEV`, `qGEV` and `rGEV` have
  been moved (with some changes) to the **nieve** package, available
  both [on CRAN](N.R-project.org/package=nieve)
  and [on GitHub](https://github.com/yvesdeville/nieve/). So **NSGEV** 
  now imports these functions from **nieve**. As a consequence the
  compiled code used by this functions has been discarded, making the
  package easier to install at least for Windows users.
  
  

