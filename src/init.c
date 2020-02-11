#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void mcmcloop(void *, void *,void *,void *,void *,void *,void *,void *,void *,
                     void *, void *,void *,void *,void *,void *,void *,void *,void *,
                     void *, void *,void *,void *,void *,void *,void *,void *,void *,
                     void *, void *,void *,void *,void*);

static const R_CMethodDef CEntries[] = {
    {"mcmcloop", (DL_FUNC) &mcmcloop, 32},
    {NULL, NULL, 0}
};

void R_init_HDCurves(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
