#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP Call_dGEV(SEXP x,
	       SEXP loc, SEXP scale, SEXP shape,
	       SEXP logFlag, SEXP derivFlag, SEXP hessianFlag);
  
SEXP Call_pGEV(SEXP q, 
	       SEXP loc, SEXP scale, SEXP shape, 
	       SEXP lowerTailFlag, SEXP derivFlag);
 
SEXP Call_qGEV(SEXP p,
	       SEXP loc, SEXP scale, SEXP shape,
	       SEXP lowerTailFlag, SEXP derivFlag, SEXP hessianFlag); 
 
