#' @title An S4 class to represent estimated optimal lambdas
#'
#' @description An S4 class to represent estimated optimal lambdas resulted from runnig \code{\link{fmrs.tunsel}}.
#' @name fmrs.lambda-class
#' @docType class
#' @exportClass fmrs.lambda
#' @import methods
#' @slot lamPen A length-\code{nComp} numeric vector
#' @slot disFamily A length-one character vector
#' @slot penFamily A length-one character vector
#' @slot lamRidge A length-one numeric vector
#' @slot method A length-one character vector
#' @slot data A list including \code{y}, \code{x} and \code{delta}
#' @export
frms.lambda <- setClass("fmrs.lambda",
                        representation(lamPen = "vector",
                                       disFamily = "character",
                                       penFamily = "character",
                                       lamRidge = "numeric",
                                       method = "character",
                                       data = "list"
                        )
)
