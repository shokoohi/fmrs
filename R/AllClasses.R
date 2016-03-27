#' @title An S4 class to represent a fit of an FMRs model
#'
#' @description fmrs is an S4 class represents a fit of FMRs models
#'     resulted from running \code{\link{fmrsmle}}
#'     or \code{\link{fmrsvarsel}}
#' @name fmrs-class
#' @import methods
#' @slot dims A length-three numeric vector represents number of observations
#'     (\code{n}), number of covariates (\code{nCov}) and the order of the
#'     mixture model (\code{nComp})
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
#' @exportClass fmrs
frms <- setClass("fmrs",
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

#' @title An S4 class to represent estimated optimal lambdas
#'
#' @description An S4 class to represent estimated optimal lambdas resulted
#'     from runnig \code{\link{fmrstunsel}}.
#' @name tunepar-class
#' @docType class
#' @exportClass tunepar
#' @import methods
#' @slot lamPen A length-\code{nComp} numeric vector
#' @slot disFamily A length-one character vector
#' @slot penFamily A length-one character vector
#' @slot lamRidge A length-one numeric vector
#' @slot method A length-one character vector
#' @slot data A list including \code{y}, \code{x} and \code{delta}
#' @export
tunepar <- setClass("tunepar",
                        representation(lamPen = "vector",
                                       disFamily = "character",
                                       penFamily = "character",
                                       lamRidge = "numeric",
                                       method = "character",
                                       data = "list"
                        )
)
