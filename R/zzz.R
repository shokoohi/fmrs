#' @useDynLib fmrs

.onAttach <- function(lib, pkg){
  fmrsversion <- utils::packageVersion("fmrs")
  fmrsdate <- utils::packageDescription("fmrs")$Date
  fmrsdescrip <- utils::packageDescription("fmrs")$Description
  packageStartupMessage(
    paste('fmrs package, Version ', fmrsversion,', Released ',fmrsdate,' \n',fmrsdescrip, sep = "")
  )
}

.onUnload <- function(libpath)
  library.dynam.unload("fmrs", libpath)
