#' Get Information of zinbwave
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @importFrom splatter newSimpleParams
#' @export
#'
#' @examples
#' zinbwave_method_definition <- zinbwave_method_definition()

zinbwave_method_definition <- function(...){
  zinbwave_parameters <- parameter_sets(
    param_reference(
      id = "counts",
      type = c("matrix"),
      default = NULL,
      force = TRUE,
      process = "estimation",
      description = "A counts matrix containing count data to estimate parameters from.",
      function_name = "zinbEstimate"
    ),
    param_others(
      id = "design.samples",
      type = "matrix",
      default = NULL,
      process = "estimation",
      description = "Design matrix of sample-level covariates.",
      function_name = "zinbEstimate"
    ),
    param_others(
      id = "design.genes",
      type = "matrix",
      default = NULL,
      process = "estimation",
      description = "Design matrix of gene-level covariates.",
      function_name = "zinbEstimate"
    ),
    param_Boolean(
      id = "common.disp",
      default = TRUE,
      description = "Logical. Whether or not a single dispersion for all features is estimated.",
      function_name = "zinbEstimate"
    ),
    param_integer(
      id = "iter.init",
      default = 2L,
      lower = 1L,
      description = "The number of iterations to use for initialization.",
      function_name = "zinbEstimate"
    ),
    param_integer(
      id = "iter.opt",
      default = 25L,
      lower = 1L,
      description = "The number of iterations to use for optimization.",
      function_name = "zinbEstimate"
    ),
    param_numeric(
      id = "stop.opt",
      default = 1e-04,
      lower = 0,
      border = FALSE,
      description = "Stopping criterion for optimization.",
      function_name = "zinbEstimate"
    ),
    param_others(
      id = "params",
      default = "newZINBParams()",
      type = "newZINBParams",
      description = "ZINBParams object to store estimated values in.",
      function_name = "zinbEstimate"
    ),
    param_Boolean(
      id = "verbose",
      default = TRUE,
      description = "Logical. Whether to print progress messages.",
      function_name = c("zinbEstimate", "zinbSimulate")
    ),
    param_others(
      id = "BPPARAM",
      default = "SerialParam()",
      type = "SerialParam",
      description = "	A BiocParallelParam instance giving the parallel back-end to be used. Default is SerialParam which uses a single core.",
      function_name = "zinbEstimate"
    ),
    param_others(
      id = "params",
      type = "ZINBParams",
      default = "ZINBParams()",
      description = "ZINBParams object containing simulation parameters.",
      function_name = "zinbSimulate"
    ),
    param_Boolean(
      id = "sparsify",
      default = TRUE,
      description = "logical. Whether to automatically convert assays to sparse matrices if there will be a size reduction."
    )
  )

  zinbwave_method <- method_definition(
    method = "zinbwave",
    programming = "R",
    url = "http://www.bioconductor.org/packages/release/bioc/html/zinbwave.html",
    authors = authors_definition(
      first = "Davide",
      last = "Risso",
      email = NULL,
      github = "https://github.com/drisso/zinbwave",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "A general and flexible method for signal extraction from single-cell RNA-seq data",
      doi = "10.1038/s41467-017-02554-5",
      journal = "Nature Communications",
      date = "2018",
      peer_review = TRUE
    ),
    description = "Implements a general and flexible zero-inflated negative binomial model that can be used to provide a low-dimensional representations of single-cell RNA-seq data. The model accounts for zero inflation (dropouts), over-dispersion, and the count nature of the data. The model also accounts for the difference in library sizes and optionally for batch effects and/or other covariates, avoiding the need for pre-normalize the data.",
    vignette = "http://47.254.148.113/software/Simsite/references/methods/14-zinbwave/")

  list(zinbwave_method = zinbwave_method,
       zinbwave_parameters = zinbwave_parameters)
}
