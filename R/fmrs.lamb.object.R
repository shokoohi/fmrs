#' @title An S4 class to represent estimated optimal lambdas.
#'
#' @description An S4 class to represent estimated optimal lambdas.
#' @name fmrs.lambda
#' @import methods
#' @slot lamPen A length-\code{nComp} numeric vector
#' @slot disFamily A length-one character vector
#' @slot penFamily A length-one character vector
#' @slot lamRidge A length-one numeric vector
#' @slot method A length-one character vecotr
#' @slot data A list including \code{y}, \code{x} and \code{delta}
frms.lambda <- setClass("fmrs.lambda",
                     slots = list(lamPen = "numeric",
                                  disFamily = "character",
                                  penFamily = "character",
                                  lamRidge = "numeric",
                                  method = "character",
                                  data = "list"
                     )
)
