#' Get Information of CancerInSilico
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' CancerInSilico_method_definition <- CancerInSilico_method_definition()
CancerInSilico_method_definition <- function(...){

  CancerInSilico_parameters <- parameter_sets(
    param_reference(
      id = "ref_data",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      description = "A counts.",
      function_name = "calibratePathway"
    ),
    param_others(
      id = "model",
      type = "CellModel",
      description = "A CellModel object.",
      function_name = "inSilicoGeneExpression"
    ),
    param_others(
      id = "pathways",
      type = "list",
      description = "List of genes pathways.",
      function_name = "inSilicoGeneExpression"
    ),
    param_others(
      id = "params",
      type = "GeneExpressionParams",
      description = "GeneExpressionParams object.",
      function_name = "inSilicoGeneExpression"
    ))

  CancerInSilico_method <- method_definition(
    method = "CancerInSilico",
    programming = "R",
    url = "https://www.bioconductor.org/packages/release/bioc/html/CancerInSilico.html",
    authors = authors_definition(
      first = "Thomas D.",
      last = "Sherman",
      email = NULL,
      github = "https://github.com/FertigLab/CancerInSilico",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "CancerInSilico: An R/Bioconductor package for combining mathematical and statistical modeling to simulate time course bulk and single cell gene expression data in cancer",
      doi = "10.1371/journal.pcbi.1006935",
      journal = "PLoS Computational Biology",
      date = "2019",
      peer_review = TRUE
    ),
    description = "The CancerInSilico package provides an R interface for running mathematical models of tumor progresson and generating gene expression data from the results.",
    vignette = "http://47.254.148.113/software/Simsite/references/methods/25-cancerinsilico/")

  list(CancerInSilico_method = CancerInSilico_method,
       CancerInSilico_parameters = CancerInSilico_parameters)
}
