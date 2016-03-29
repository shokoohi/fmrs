#' @rdname weights-methods
#' @rdname fmrsfit-class
#' @aliases weights,fmrsfit-method
setMethod("weights", signature = "fmrsfit", weights.fmrsfit)

#' @rdname residuals-methods
#' @rdname fmrsfit-class
#' @aliases residuals,fmrsfit-method
setMethod("residuals", signature = "fmrsfit", residuals.fmrsfit)

#' @rdname nobs-methods
#' @rdname fmrsfit-class
#' @aliases nobs,fmrsfit-method
setMethod("nobs", signature = "fmrsfit", nobs.fmrsfit)

#' @rdname ncov-methods
#' @rdname fmrsfit-class
#' @aliases ncov,fmrsfit-method
setMethod("ncov", signature = "fmrsfit", ncov.fmrsfit)

#' @rdname ncomp-methods
#' @rdname fmrsfit-class
#' @aliases ncomp,fmrsfit-method
setMethod("ncomp", signature = "fmrsfit", ncomp.fmrsfit)

#' @rdname mixProp-methods
#' @rdname fmrsfit-class
#' @aliases mixProp,fmrsfit-method
setMethod("mixProp", signature = "fmrsfit", mixProp.fmrsfit)

#' @rdname logLik-methods
#' @rdname fmrsfit-class
#' @aliases logLik,fmrsfit-method
setMethod("logLik", signature = "fmrsfit", logLik.fmrsfit)

#' @rdname fitted-methods
#' @rdname fmrsfit-class
#' @aliases fitted,fmrsfit-method
setMethod("fitted", signature = "fmrsfit", fitted.fmrsfit)

#' @rdname deviance-methods
#' @rdname fmrsfit-class
#' @aliases deviance,fmrsfit-method
setMethod("deviance", signature = "fmrsfit", deviance.fmrsfit)

#' @rdname coefficients-methods
#' @rdname fmrsfit-class
#' @aliases coefficients,fmrsfit-method
setMethod("coefficients", signature = "fmrsfit", coefficients.fmrsfit)

#' @rdname BIC-methods
#' @rdname fmrsfit-class
#' @aliases BIC,fmrsfit-method
setMethod("BIC", signature = "fmrsfit", BIC.fmrsfit)

#' @rdname summary-methods
#' @rdname fmrsfit-class
#' @aliases summary,fmrsfit-method
setMethod("summary", signature = "fmrsfit", summary.fmrsfit)

#' @rdname show-methods
#' @rdname fmrsfit-class
#' @aliases show,fmrsfit-method
setMethod(f = "show", signature = "fmrsfit", show.fmrsfit)

