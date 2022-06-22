#' Get Information of scDD
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' scDD_method_definition <- scDD_method_definition()
#'
scDD_method_definition <- function(...){

  scDD_parameters <- parameter_sets(
    param_reference(
      id = "counts",
      type = c("matrix", "SingleCellExperiment"),
      default = NULL,
      process = "estimation",
      force = TRUE,
      description = "Either a counts matrix or a SingleCellExperiment object containing count data to estimate parameters from.",
      function_name = "scDDEstimate"
    ),
    param_others(
      id = "params",
      type = "newSCDDParams",
      default = "splatter::newSCDDParams()",
      process = "estimation",
      description = "SCDDParams object to store estimated values in.",
      function_name = "scDDEstimate"
    ),
    param_Boolean(
      id = "verbose",
      default = TRUE,
      process = "estimation",
      description = "Logical. Whether to show progress messages.",
      function_name = c("scDDEstimate", "scDDSimulate")
    ),
    param_others(
      id = "BPPARAM",
      type = "BiocParallelParam",
      default = "SerialParam()",
      process = "estimation",
      description = "A BiocParallelParam instance giving the parallel back-end to be used. Default is SerialParam which uses a single core.",
      function_name = c("scDDEstimate", "scDDSimulate")
    ),
    param_vector(
      id = "condition",
      force = TRUE,
      process = "estimation",
      description = "String giving the column that represents biological group of interest.",
      function_name = "scDDEstimate"
    ),
    param_others(
      id = "params",
      type = "SCDDParams",
      force = TRUE,
      process = "simulation",
      description = "SCDDParams object containing simulation parameters.",
      function_name = "scDDSimulate"
    ),
    param_Boolean(
      id = "plots",
      default = FALSE,
      process = "simulation",
      description = "logical. whether to generate scDD fold change and validation plots.",
      function_name = "scDDSimulate"
    ),
    param_Boolean(
      id = "sparsify",
      default = TRUE,
      process = "simulation",
      description = "logical. Whether to automatically convert assays to sparse matrices if there will be a size reduction.",
      function_name = "scDDSimulate"
    )
  )

  scDD_method <- method_definition(
    method = "scDD",
    programming = "R",
    url = "https://www.bioconductor.org/packages/release/bioc/html/scDD.html",
    authors = authors_definition(
      first = "Keegan",
      last = "Korthauer",
      email = NULL,
      github = "https://github.com/kdkorthauer/scDD",
      orcid = "0000-0002-4565-1654"
    ),
    manuscript = manuscript_definition(
      title = "A statistical approach for identifying differential distributions in single-cell RNA-seq experiments",
      doi = "10.1186/s13059-016-1077-y",
      journal = "Genome Biology",
      date = "2016",
      peer_review = TRUE
    ),
    description = "This package implements a method to analyze single-cell RNA- seq Data utilizing flexible Dirichlet Process mixture models. Genes with differential distributions of expression are classified into several interesting patterns of differences between two conditions. The package also includes functions for simulating data with these patterns from negative binomial distributions.")

  list(scDD_method = scDD_method,
       scDD_parameters = scDD_parameters)
}
