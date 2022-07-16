.onAttach <- function(libname, pkgname) {
  if("roxygen2" %in% (.packages())){
    suppressWarnings(detach("package:roxygen2", unload = FALSE, force = TRUE))
  }
  packageStartupMessage("Welcome to simmethods")
}
