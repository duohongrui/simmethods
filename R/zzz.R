.onAttach <- function(libname, pkgname) {
  suppressWarnings(detach("package:roxygen2", unload = FALSE, force = TRUE))
  packageStartupMessage("Welcome to simmethods")
}
