#' @title An S4 class to represent a fit of an FMRs model
#'
#' @description fmrs.fit is an S4 class to represent a fit of FMRs models resulted from running \code{\link{fmrs.mle}} or \code{\link{fmrs.varsel}}.
#' @name fmrs.fit-class
#' @import methods
#' @slot dims A length-three numeric vector represents number of observations (\code{n}), number of covariates (\code{nCov}) and the order of the mixture model (\code{nComp})
#' @slot coefficients A dimension-\code{nCov}-\code{nComp} numeric matrix
#' @slot deviance A length-\code{nComp} numeric vector
#' @slot pi A length-\code{nComp} numeric vector
#' @slot logLik A length-one numeric vector
#' @slot BIC A length-one numeric vector
#' @slot nIterEMconv A length-one numeric vector
#' @slot disFamily A length-one character vector
#' @slot penFamily A length-one character vector
#' @slot lambPen A length-\code{nComp} numeric vector
#' @slot lamRidge A length-one numeric vector
#' @slot method A length-one character vector
#' @slot fitted A dimension-\code{n}-\code{nComp} numeric matrix
#' @slot residuals A dimension-\code{n}-\code{nComp} numeric matrix
#' @slot weights A dimension-\code{n}-\code{nComp} numeric matrix
#' @slot data A list including \code{y}, \code{x} and \code{delta}
#' @docType class
#' @exportClass fmrs.fit
frms.fit <- setClass("fmrs.fit",
                     representation(dims = "vector",
                                    coefficients = "matrix",
                                    deviance = "vector",
                                    pi = "vector",
                                    logLik = "numeric",
                                    BIC = "numeric",
                                    nIterEMconv = "numeric",
                                    disFamily = "character",
                                    penFamily = "character",
                                    lambPen = "vector",
                                    lamRidge = "numeric",
                                    method = "character",
                                    fitted = "matrix",
                                    residuals = "matrix",
                                    weights = "matrix",
                                    data = "list"
                     )

)

