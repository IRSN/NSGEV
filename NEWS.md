
**NSGEV** Package News
===========================

# News in version 0.1.7

## Bug fixes

- The `TVGEV` function did not work when a constant location 
  was specified along with a non-constant scale and/or shape. 
  Thanks to Jesper Rydén.

# News in version 0.1.6

## Changes

- The GEV probability functions `dGEV`, `pGEV`, `qGEV` and `rGEV` have
  been moved (with some changes) to the **nieve** package, available
  both [on CRAN](N.R-project.org/package=nieve)
  and [on GitHub](https://github.com/yvesdeville/nieve/). So **NSGEV** 
  now imports these functions from ** nieve**. As a consequence the
  compiled code used by this functions has been discarded, making the
  package easier to install at least for Windows users.
  
  

