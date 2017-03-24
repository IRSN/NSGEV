#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#define NODEBUG

SEXP Call_dGEV(SEXP x,            /*  double                          */
	       SEXP loc,          /*  double                          */
	       SEXP scale,        /*  double                          */
	       SEXP shape,        /*  double                          */
	       SEXP logFlag,      /*  integer                         */
	       SEXP derivFlag) {  /*  integer                         */
 
  int n, nx, nloc, nscale, nshape, i, ix, iloc, iscale, ishape;
  
  double eps = 1e-6, z, emz, V, xi;
  
  SEXP val;
  
  x = coerceVector(x, REALSXP);
  loc = coerceVector(loc, REALSXP);
  scale = coerceVector(scale, REALSXP);
  shape = coerceVector(shape, REALSXP);

  double *rx = REAL(x), *rloc = REAL(loc), *rscale = REAL(scale), *rshape = REAL(shape);

  nx = LENGTH(x);						
  nloc = LENGTH(loc);						
  nscale = LENGTH(scale);		
  nshape = LENGTH(shape);
  
  if ((nx == 0) || (nloc == 0) || (nscale == 0) || (nshape == 0)) 			
    return(allocVector(REALSXP, 0));				
  
  n = nx;							
  if (n < nloc) n = nloc;						
  if (n < nscale) n = nscale;
  if (n < nshape) n = nshape;
  
  PROTECT(val = allocVector(REALSXP, n));
  double *rval = REAL(val);
  
  if (INTEGER(derivFlag)[0]) {

    double U, W;
    SEXP grad, attrNm;

    PROTECT(grad = allocVector(REALSXP, n * 3));
    double *rgrad = REAL(grad);

    PROTECT(attrNm = NEW_CHARACTER(1)); 
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));

    for (i = ix = iloc = iscale = ishape = 0;  i < n; 
	 ix = (++ix == nx) ? 0 : ix, 	       
	 iloc = (++iloc == nloc) ? 0 : iloc, 
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {
      
      if (ISNA(rx[ix]) || (rscale[iscale] <= 0.0)) {
	
	rval[i] = NA_REAL;
	
	rgrad[i] = NA_REAL;
	rgrad[i + n] = NA_REAL;
	rgrad[i + 2 * n] = NA_REAL; 
	

      } else {
	
	z = (rx[ix] - rloc[iloc]) / rscale[iscale];
	xi = rshape[ishape];

	if (fabs(xi) < eps) {
	  
	  emz = exp(-z);
	  rval[i] = -log(rscale[iscale]) - z - emz;
	  
	  rgrad[i] = (1.0 - emz) / rscale[iscale];
	  rgrad[i + n] =  (-1.0 + z * (1.0 - emz)) / rscale[iscale];
	  rgrad[i + 2 * n] = z * z * (1 - emz) / 2.0 - z;
	  
	} else {
	  
	  V = 1.0 + xi * z;
	  // Rprintf("%d, %d, %6.3f, %6.3f\n", i, ishape, xi, V);

	  if (V > 0.0) {
	
	    rval[i] = -log(rscale[iscale]) - pow(V, - 1.0 / xi) -
	      (1.0 / xi + 1.0) * log(V);	  
	    
	    W = pow(V, -1.0 / xi);
	    U = (1.0 + xi - W) / V / rscale[iscale];
	    
	    rgrad[i] = U;
	    rgrad[i + n] = -1.0 / rscale[iscale] + z * U;
	    rgrad[i + 2 * n] = log(V) * (1.0 - W) / xi / xi - 
	      z * U * rscale[iscale] / xi;
	  
	  } else {

	    rval[i] = R_NegInf;
	    rgrad[i] = 0.0;
	    rgrad[i + n] = 0.0;
	    rgrad[i + 2 * n] = 0.0;

	  }

	} /* non-Gumbel case */
	
	if (!INTEGER(logFlag)[0]) {
	  rval[i] = exp(rval[i]);
	  rgrad[i] *= rval[i];
	  rgrad[i + n] *= rval[i];
	  rgrad[i + 2 * n] *= rval[i];
	}
	
      }   /* non-NA case     */
      
    }
  
    SET_ATTR(val, attrNm, grad);
    UNPROTECT(3);
    return(val);
    
  } else {

    for (i = ix = iloc = iscale = ishape = 0;  i < n; 
	 ix = (++ix == nx) ? 0 : ix, 	       
	 iloc = (++iloc == nloc) ? 0 : iloc, 
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {
      
      if (ISNA(rx[ix]) || (rscale[iscale] <= 0.0)) {
	
	rval[i] = NA_REAL;
		
      } else {
	
	z = (rx[ix] - rloc[iloc]) / rscale[iscale];
	xi = rshape[ishape];

	if (fabs(xi) < eps) {
	  
	  emz = exp(-z);
	  rval[i] = -log(rscale[iscale]) - z - emz;	  
	  
	} else {
	  
	  V = 1.0 + xi * z;
	  // Rprintf("%d, %d, %6.3f, %6.3f\n", i, ishape, xi, V);

	  if (V > 0.0) {
	    rval[i] = -log(rscale[iscale]) - pow(V, - 1.0 / xi) -
	      (1.0 / xi + 1.0) * log(V);
	  } else {
	    rval[i] = R_NegInf;
	  }
	  
	  
	} /* non-Gumbel case */
	
	if (!INTEGER(logFlag)[0]) {
	  rval[i] = exp(rval[i]);
	}
	
      }   /* non-NA case     */
      
    }
  
    UNPROTECT(1);
    return(val);
    
  }

}
