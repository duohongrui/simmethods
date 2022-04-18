#' Get Information of Splat
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' splat_method_definition <- splat_method_definition()
splat_method_definition <- function(...){

  splat_parameters <- parameter_sets(
    param_reference(
      id = "input_refernece",
      type = c("matrix", "SingleCellExperiment"),
      process = "estimation",
      description = "The reference data to go through estimation and evaluation process. Usually, no reference data is also OK."
    ),
    param_others(
      id = "splatParams",
      type = "SplatParams",
      process = "estimation",
      description = "Usually it is default by splatter::newSplatParams function. Users can change the parameters by splatter::setParam function."
    ),
    param_integer(
      id = "nBatches",
      default = 1L,
      lower = 1L
    ),
    param_integer(
      id = "batchCells",
      default = 100L,
      lower = 1L
    ),
    param_numeric(
      id = "batch.facLoc",
      default = 0.1,
      lower = 0,
      border = FALSE,
      upper = 1
    ),
    param_numeric(
      id = "batch.facScale",
      default = 0.1,
      lower = 0,
      border = FALSE,
      upper = 1
    ),
    param_Boolean(
      id = "batch.rmEffect",
      default = FALSE
    ),
    param_numeric(
      id = "mean.shape",
      default = 0.6,
      lower = 0,
      border = FALSE,
    ),
    param_numeric(
      id = "mean.rate",
      default = 0.3,
      lower = 0,
      border = FALSE,
    ),
    param_integer(
      id = "lib.loc",
      default = 11L,
      lower = 1L
    ),
    param_numeric(
      id = "lib.scale",
      default = 0.2,
      lower = 0,
      border = FALSE
    ),
    param_Boolean(
      id = "lib.norm",
      default = FALSE
    ),
    param_numeric(
      id = "out.prob",
      default = 0.05,
      lower = 0,
      upper = 1,
      border = FALSE
    ),
    param_numeric(
      id = "out.facLoc",
      default = 4L,
      lower = 0,
      border = FALSE
    ),
    param_numeric(
      id = "out.facScale",
      default = 0.5,
      lower = 0,
      border = FALSE,
      upper = 1
    ),
    param_integer(
      id = "nGroups",
      default = 1L,
      lower = 1L
    ),
    param_numeric(
      id = "group.prob",
      default = 1,
      lower = 0,
      border = TRUE,
      upper = 1
    ),
    param_numeric(
      id = "de.prob",
      default = 0.1,
      lower = 0,
      border = TRUE,
      upper = 1
    ),
    param_numeric(
      id = "de.downProb",
      default = 0.5,
      lower = 0,
      border = TRUE,
      upper = 1
    ),
    param_numeric(
      id = "de.facLoc",
      default = 0.1,
      lower = 0,
      border = FALSE,
    ),
    param_numeric(
      id = "de.facScale",
      default = 0.4,
      lower = 0,
      border = FALSE,
    ),
    param_numeric(
      id = "bcv.common",
      default = 0.1,
      lower = 0,
      border = FALSE,
    ),
    param_integer(
      id = "bcv.df",
      default = 60L,
      lower = 1L
    ),
    param_character(
      id = "dropout.type",
      default = "none",
      alternatives = c("none", "experiment", "batch", "group", "cell")
    ),
    param_numeric(
      id = "dropout.mid",
      default = 0,
    ),
    param_numeric(
      id = "dropout.shape",
      default = -1,
    ),
    param_vector(
      id = "path.from",
      default = 0
    ),
    param_integer(
      id = "path.nSteps",
      default = 100L,
      lower = 1L
    ),
    param_numeric(
      id = "path.skew",
      default = 0.5,
      lower = 0,
      upper = 1
    ),
    param_numeric(
      id = "path.nonlinearProb",
      default = 0.1,
      lower = 0,
      upper = 1
    ),
    param_numeric(
      id = "path.sigmaFac",
      default = 0.8,
      lower = 0,
      upper = 1
    ),
    param_integer(
      id = "nGenes",
      default = 10000L,
      lower = 1L
    ),
    param_integer(
      id = "nCells",
      default = 100L,
      lower = 1L
    ),
    param_integer(
      id = "seed",
      default = 687680L
    )
  )

  splat_method <- method_definition(
    method = "splat",
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

  list(splat_method = splat_method,
       splat_parameters = splat_parameters)
}
