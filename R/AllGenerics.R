#' @title Method nobs
#' @description nobs method
#' @name nobs
#' @rdname nobs-methods
#' @param object An \code{\link{fmrsfit-class}} object
#' @param ... other possible arguments
#' @return An integer value
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' deviance = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDeviance = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' nobs(res.mle)
#' @exportMethod nobs
setGeneric("nobs", function(object, ...) standardGeneric("nobs"))
nobs.fmrsfit <- function(object, ...) {object@nobs}

#' @title Method ncov
#' @description ncov method
#' @name ncov
#' @rdname ncov-methods
#' @param object An \code{\link{fmrsfit-class}} object
#' @param ... other possible arguments
#' @return An integer value
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' deviance = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDeviance = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' ncov(res.mle)
#' @exportMethod ncov
setGeneric("ncov", function(object, ...) standardGeneric("ncov"))
ncov.fmrsfit <- function(object, ...) {object@ncov}

#' @title Method ncomp
#' @description ncomp method
#' @name ncomp
#' @rdname ncomp-methods
#' @param object An \code{\link{fmrsfit-class}} object
#' @param ... other possible arguments
#' @return An integer value
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' deviance = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDeviance = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' ncomp(res.mle)
#' @exportMethod ncomp
setGeneric("ncomp", function(object, ...) standardGeneric("ncomp"))
ncomp.fmrsfit <- function(object, ...) {object@ncomp}

#' @title Method coefficients
#' @description coefficients method
#' @name coefficients
#' @rdname coefficients-methods
#' @param object An \code{\link{fmrsfit-class}} object
#' @param ... other possible arguments
#' @return A numeric array
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' deviance = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDeviance = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' coefficients(res.mle)
#' @exportMethod coefficients
setGeneric("coefficients",
           function(object, ...) standardGeneric("coefficients"))
coefficients.fmrsfit <- function(object, ...) {object@coefficients}

#' @title Method deviance
#' @description deviance method
#' @name deviance
#' @rdname deviance-methods
#' @param object An \code{\link{fmrsfit-class}} object
#' @param ... other possible arguments
#' @return A numeric array
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' deviance = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDeviance = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' deviance(res.mle)
#' @exportMethod deviance
setGeneric("deviance", function(object, ...) standardGeneric("deviance"))
deviance.fmrsfit <- function(object, ...) {object@deviance}

#' @title Method mixProp
#' @description mixProp method
#' @name mixProp
#' @rdname mixProp-methods
#' @param object An \code{\link{fmrsfit-class}} object
#' @param ... other possible arguments
#' @return A numeric array
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' deviance = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDeviance = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' mixProp(res.mle)
#' @exportMethod mixProp
setGeneric("mixProp", function(object, ...) standardGeneric("mixProp"))
mixProp.fmrsfit <- function(object, ...) {object@mixProp}

#' @title Method fitted
#' @description fitted method
#' @name fitted
#' @rdname fitted-methods
#' @param object An \code{\link{fmrsfit-class}} object
#' @param ... other possible arguments
#' @return A numeric array
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' deviance = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDeviance = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' head(fitted(res.mle))
#' @exportMethod fitted
setGeneric("fitted", function(object, ...) standardGeneric("fitted"))
fitted.fmrsfit <- function(object, ...) {object@fitted}

#' @title Method residuals
#' @description residuals method
#' @name residuals
#' @rdname residuals-methods
#' @param object An \code{\link{fmrsfit-class}} object
#' @param ... other possible arguments
#' @return A numeric array
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' deviance = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDeviance = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' head(residuals(res.mle))
#' @exportMethod residuals
setGeneric("residuals", function(object, ...) standardGeneric("residuals"))
residuals.fmrsfit <- function(object, ...) {object@residuals}

#' @title Method weights
#' @description weights method
#' @name weights
#' @rdname weights-methods
#' @param object An \code{\link{fmrsfit-class}} object
#' @param ... other possible arguments
#' @return A numeric array
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' deviance = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDeviance = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' head(weights(res.mle))
#' @exportMethod weights
setGeneric("weights", function(object, ...) standardGeneric("weights"))
weights.fmrsfit <- function(object, ...) {object@weights}

#' @title Method logLik
#' @description logLik method
#' @name logLik
#' @rdname logLik-methods
#' @param object An \code{\link{fmrsfit-class}} object
#' @param ... other possible arguments
#' @return A numeric value
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' deviance = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDeviance = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' logLik(res.mle)
#' @exportMethod logLik
setGeneric("logLik", function(object, ...) standardGeneric("logLik"))
logLik.fmrsfit <- function(object, ...) {object@logLik}

#' @title Method BIC
#' @description BIC method
#' @name BIC
#' @rdname BIC-methods
#' @param object An \code{\link{fmrsfit-class}} object
#' @param ... other possible arguments
#' @return A numeric value
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' deviance = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDeviance = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' BIC(res.mle)
#' @exportMethod BIC
setGeneric("BIC", function(object, ...) standardGeneric("BIC"))
BIC.fmrsfit <- function(object, ...) {object@BIC}

#' @title Method summary
#' @description summary method
#' @name summary
#' @rdname summary-methods
#' @param object An \code{\link{fmrsfit-class}} object
#' @param ... other possible arguments
#' @return A fitted fmrs model
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' deviance = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDeviance = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' summary(res.mle)
#' @exportMethod summary
setGeneric("summary", function(object, ...) standardGeneric("summary"))
summary.fmrsfit <- function(object, ...) {
  if(object@model == "FMR") {
    modelfmr = "Finite Mixture of Regression Models"
  }else if(object@disFamily == "lnorm"){
    modelfmr = "Finite Mixture of Accelerated Failure Time Regression
    Models \n  Log-Normal Sub-Distributions"
  }else{
    modelfmr = "Finite Mixture of Accelerated Failure Time Regression
    Models \n  Weibull Sub-Distributions"
  }
  cat("-------------------------------------------\n")
  cat("Fitted Model: \n")
  cat("-------------------------------------------\n")
  cat(" ", modelfmr, "\n")
  cat("  ", object@ncomp, " Components; ",
      object@ncov," Covariates; ", object@nobs,
      " samples.", sep = "")
  cat("\n\nCoefficients:\n")
  print.default(object@coefficients)
  cat("\nDeviances:\n")
  print.default(object@deviance)
  cat("\nMixing Proportions:\n")
  print.default(object@mixProp)
  cat("\nLogLik: ", object@logLik, "; BIC: ", object@BIC, sep="")
  }

#' @title Method show
#' @description show method
#' @name show
#' @rdname show-methods
#' @param object An \code{\link{fmrsfit-class}} object
#' @return A fitted fmrs model
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' deviance = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDeviance = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' show(res.mle)
#' @exportMethod show
setGeneric("show")

show.fmrsfit <- function(object) {
  if(object@model == "FMR") {
    modelfmr = "Finite Mixture of Regression Models"
  }else if(object@disFamily == "lnorm"){
    modelfmr = "Finite Mixture of Accelerated Failure Time Regression Models
    Log-Normal Sub-Distributions"
  }else{
    modelfmr = "Finite Mixture of Accelerated Failure Time Regression Models
    Weibull Sub-Distributions"
  }
  cat("An object of class '", class(object), "'\n", sep = "")
  cat(" ", modelfmr, "\n")
  cat("  ", object@ncomp, " Components; ",
      object@ncov," Covariates; ", object@nobs,
      " samples.\n", sep = "")
  }





validity.fmrsfit <- function(object) {
  if(!all(sapply(object@data, is.numeric))) {
    return("data not numeric")
  } else return(TRUE)
}



