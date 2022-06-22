.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to simmethods")
  suppressPackageStartupMessages({
    require(powsimR)
    require(rtracklayer)
  })
}
