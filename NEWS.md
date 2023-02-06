
**NSGEV** Package News
===========================

# News in version 0.1.6

## Changes

- The GEV probability functions `dGEV`, `pGEV`, `qGEV` and `rGEV` have
  been moved (with some changes) to the **nieve** package, available
  both [on CRAN](https://cran.r-project.org/web/packages/nieve/index.html) 
  and [on GitHub](https://github.com/yvesdeville/nieve/). So **NSGEV** 
  now imports these functions from ** nieve**. As a consequence the
  compiled code used by this functions has been discarded, making the
  package easier to install at least for Windows users.
  
  

