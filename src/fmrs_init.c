#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void FMR_Norm_Surv_CwTuneParSel(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void FMR_Norm_Surv_EM_MLE(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void FMR_Norm_Surv_EM_VarSel(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void FMR_Weibl_Surv_CwTuneParSel(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void FMR_Weibl_Surv_EM_MLE(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void FMR_Weibl_Surv_EM_VarSel(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"FMR_Norm_Surv_CwTuneParSel",  (DL_FUNC) &FMR_Norm_Surv_CwTuneParSel,  19},
    {"FMR_Norm_Surv_EM_MLE",        (DL_FUNC) &FMR_Norm_Surv_EM_MLE,        31},
    {"FMR_Norm_Surv_EM_VarSel",     (DL_FUNC) &FMR_Norm_Surv_EM_VarSel,     35},
    {"FMR_Weibl_Surv_CwTuneParSel", (DL_FUNC) &FMR_Weibl_Surv_CwTuneParSel, 20},
    {"FMR_Weibl_Surv_EM_MLE",       (DL_FUNC) &FMR_Weibl_Surv_EM_MLE,       31},
    {"FMR_Weibl_Surv_EM_VarSel",    (DL_FUNC) &FMR_Weibl_Surv_EM_VarSel,    37},
    {NULL, NULL, 0}
};

void R_init_fmrs(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
