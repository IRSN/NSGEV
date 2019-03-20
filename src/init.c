#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP Call_dGEV(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Call_pGEV(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Call_qGEV(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"Call_dGEV", (DL_FUNC) &Call_dGEV, 7},
    {"Call_pGEV", (DL_FUNC) &Call_pGEV, 6},
    {"Call_qGEV", (DL_FUNC) &Call_qGEV, 7},
    {NULL, NULL, 0}
};

void R_init_NSGEV(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
