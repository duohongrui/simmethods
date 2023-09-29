#' Get Information of phenopath
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @importFrom splatter newPhenoParams
#' @export
#'
#' @examples
#' phenopath_method_definition <- phenopath_method_definition()
#'
phenopath_method_definition <- function(...){

  phenopath_parameters <- parameter_sets(
    param_reference(
      id = "input_refernece",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      description = "The reference data to go through estimation and evaluation process. Usually, no reference data is also OK.",
      function_name = "phenoEstimate"
    ),
    param_others(
      id = "PhenoParams",
      type = "PhenoParams",
      default = "splatter::newPhenoParams()",
      process = "estimation",
      description = "Usually it is default by splatter::newPhenoParams function. Users can change the parameters by splatter::setParam function.",
      function_name = "phenoEstimate"
    ),
    param_integer(
      id = "n.de",
      default = 2500L,
      lower = 0L,
      description = "Number of genes to simulate from the differential expression regime.",
      function_name = "phenoSimulate"
    ),
    param_integer(
      id = "n.pst",
      default = 2500L,
      lower = 0L,
      description = "Number of genes to simulate from the pseudotime regime.",
      function_name = "phenoSimulate"
    ),
    param_integer(
      id = "n.pst.beta",
      default = 2500L,
      description = "Number of genes to simulate from the pseudotime + beta interactions regime.",
      function_name = "phenoSimulate"
    ),
    param_integer(
      id = "n.de.pst.beta",
      default = 2500L,
      lower = 0L,
      description = "Number of genes to simulate from the differential expression + pseudotime + interactions regime.",
      function_name = "phenoSimulate"
    ),
    param_integer(
      id = "nGenes",
      default = 10000L,
      lower = 1L,
      description = "The number of genes to simulate.",
      function_name = "phenoSimulate"
    ),
    param_integer(
      id = "nCells",
      default = 100L,
      lower = 1L,
      description = "The number of cells to simulate.",
      function_name = "phenoSimulate"
    ),
    param_integer(
      id = "seed",
      default = 484583L,
      description = "Seed to use for generating random numbers.",
      function_name = "phenoSimulate"
    )
  )

  phenopath_method <- method_definition(
    method = "phenopath",
    programming = "R",
    url = "https://bioconductor.org/packages/release/bioc/html/phenopath.html",
    authors = authors_definition(
      first = "Kieran",
      last = "Campbell",
      email = "c.yau@bham.ac.uk",
      github = "https://github.com/kieranrcampbell/phenopath",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "Uncovering genomic trajectories with heterogeneous genetic and environmental backgrounds across single-cells and populations",
      doi = "10.1101/159913",
      journal = "bioRxiv",
      date = "2017",
      peer_review = FALSE
    ),
    description = "Genomic trajectories (pseudotimes) in the presence of heterogenous environmental and genetic backgrounds.",
    vignette = "http://47.254.148.113/software/Simsite/references/methods/37-phenopath/")

  list(phenopath_method = phenopath_method,
       phenopath_parameters = phenopath_parameters)
}
