#' @useDynLib fmrs

.onAttach <- function(lib, pkg){
  packageStartupMessage(
    paste(' fmrs package, version 1.0-1, Released 2016-03-24 \n fmrs provides parameter estimation and variable selection in Finite Mixture of Accelerated Failure Time Regression models and Finite Mixture of Regression models')
  )
}

.onUnload <- function(libpath)
  library.dynam.unload("fmrs", libpath)
