#' Get Information of Simple
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @importFrom splatter newSimpleParams
#' @export
#'
#' @examples
#' Simple_method_definition <- Simple_method_definition()

Simple_method_definition <- function(...){
  Simple_parameters <- parameter_sets(
    param_reference(
      id = "input_refernece",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      description = "The reference data to go through estimation and evaluation process. Usually, no reference data is also OK."
    ),
    param_others(
      id = "SimpleParams",
      type = "SimpleParams",
      default = "splatter::newSimpleParams()",
      process = "estimation",
      description = "Usually it is default by splatter::newSimpleParams function. Users can change the parameters by splatter::setParam function."
    ),
    param_integer(
      id = "nGenes",
      default = 10000L,
      lower = 1L,
      description = "The number of genes to simulate."
    ),
    param_integer(
      id = "nCells",
      default = 100L,
      lower = 1L,
      description = "The number of cells to simulate."
    ),
    param_integer(
      id = "seed",
      default = 687680L,
      description = "Seed to use for generating random numbers."
    ),
    param_numeric(
      id = "mean.shape",
      default = 0.6,
      lower = 0,
      border = FALSE,
      description = "The shape parameter for the mean gamma distribution."
    ),
    param_numeric(
      id = "mean.rate",
      default = 0.3,
      lower = 0,
      border = FALSE,
      description = "The rate parameter for the mean gamma distribution."
    ),
    param_numeric(
      id = "count.disp",
      default = 0.1,
      description = "The dispersion parameter for the counts negative binomial distribution."
    )
  )

  Simple_method <- method_definition(
    method = "Simple",
    programming = "R",
    url = "https://bioconductor.org/packages/release/bioc/html/splatter.html",
    authors = authors_definition(
      first = "Luke",
      last = "Zappia",
      email = NULL,
      github = "https://github.com/Oshlack/splatter",
      orcid = "0000-0001-7744-8565"
    ),
    manuscript = manuscript_definition(
      title = "Splatter: simulation of single-cell RNA sequencing data",
      doi = "10.1186/s13059-017-1305-0",
      journal = "Genome Biology",
      date = "2017",
      peer_review = TRUE
    ),
    description = "Splatter is a package for the simulation of single-cell RNA sequencing count data")

  list(Simple_method = Simple_method,
       Simple_parameters = Simple_parameters)
}