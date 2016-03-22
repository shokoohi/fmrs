#' @title An S4 class to represent a fit from an FMRs model.
#'
#' @description An S4 class to represent a fit form an FMRs model
#' @name fmrs.fit
#' @import methods
#' @slot dims A length-three numeric vector
#' @slot coefficients A dimension-\code{nCov}-\code{nComp} numeric matrix
#' @slot sigma A length-\code{nComp} numeric vector
#' @slot pi A length-\code{nComp} numeric vector
#' @slot logLik A length-one numeric vector
#' @slot BIC A length-one numeric vector
#' @slot nIterEMconv A length-one numeric vector
#' @slot disFamily A length-one character vector
#' @slot penFamily A length-one character vector
#' @slot lambPen A length-\code{nComp} numeric vector
#' @slot lamRidge A length-one numeric vector
#' @slot method A length-one character vecotr
#' @slot fitted A length-\code{n} numeric vecotr
#' @slot residuals A length-\code{n} numeric vector
#' @slot data A list including \code{y}, \code{x} and \code{delta}
frms.fit <- setClass("fmrs.fit",
                     slots = list(dims = "numeric",
                                  coefficients = "matrix",
                                  sigma = "numeric",
                                  pi = "numeric",
                                  logLik = "numeric",
                                  BIC = "numeric",
                                  nIterEMconv = "numeric",
                                  disFamily = "character",
                                  penFamily = "character",
                                  lambPen = "numeric",
                                  lamRidge = "numeric",
                                  method = "character",
                                  fitted = "numeric",
                                  residuals = "numeric",
                                  data = "list"
                     )
)
