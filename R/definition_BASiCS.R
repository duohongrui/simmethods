#' Get Information of BASiCS
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' BASiCS_method_definition <- BASiCS_method_definition()
#'
BASiCS_method_definition <- function(...){

  BASiCS_parameters <- parameter_sets(
    param_reference(
      id = "counts",
      type = c("matrix", "SingleCellExperiment"),
      default = NULL,
      process = "estimation",
      force = TRUE,
      description = "Either a counts matrix or a SingleCellExperiment object containing count data to estimate parameters from.",
      function_name = "BASiCSEstimate"
    ),
    param_dataframe(
      id = "spike.info",
      process = "estimation",
      description = "data.frame describing spike-ins with two columns: 'Name' giving the names of the spike-in features (must match rownames(counts)) and 'Input' giving the number of input molecules.",
      function_name = "BASiCSEstimate"
    ),
    param_vector(
      id = "batch",
      process = "estimation",
      force = TRUE,
      description = "Vector giving the batch that each cell belongs to.",
      function_name = "BASiCSEstimate"
    ),
    param_numeric(
      id = "n",
      default = 20000,
      process = "estimation",
      description = "Total number of MCMC iterations. Must be >= max(4, thin) and a multiple of thin.",
      function_name = "BASiCSEstimate"
    ),
    param_numeric(
      id = "thin",
      default = 10,
      process = "estimation",
      description = "Thining period for the MCMC sampler. Must be >= 2.",
      function_name = "BASiCSEstimate"
    ),
    param_numeric(
      id = "burn",
      default = 5000,
      process = "estimation",
      description = "Burn-in period for the MCMC sampler. Must be in the range 1 <= burn < n and a multiple of thin.",
      function_name = "BASiCSEstimate"
    ),
    param_Boolean(
      id = "params",
      default = TRUE,
      process = "estimation",
      description = "Logical. Whether to use regression to identify over-dispersion. See BASiCS_MCMC for details.",
      function_name = "BASiCSEstimate"
    ),
    param_Boolean(
      id = "verbose",
      default = TRUE,
      process = "estimation",
      description = "Logical. Whether to print progress messages.",
      function_name = "BASiCSEstimate"
    ),
    param_Boolean(
      id = "progress",
      default = TRUE,
      process = "estimation",
      description = "Logical. Whether to print additional BASiCS progress messages.",
      function_name = "BASiCSEstimate"
    ),
    param_others(
      id = "params",
      type = "BASiCSParams",
      default = "newBASiCSParams()",
      description = "BASiCSParams object to store estimated values in.",
      function_name = c("BASiCSEstimate", "BASiCSSimulate")
    ),
    param_Boolean(
      id = "sparsify",
      default = TRUE,
      description = "Logical. Whether to automatically convert assays to sparse matrices if there will be a size reduction.",
      function_name = "BASiCSSimulate"
    )
  )

  BASiCS_method <- method_definition(
    method = "BASiCS",
    programming = "R",
    url = "https://bioconductor.org/packages/release/bioc/html/BASiCS.html",
    authors = authors_definition(
      first = "Catalina",
      last = "Vallejos",
      email = "catalina@mrc-bsu.cam.ac.uk",
      github = "https://github.com/catavallejos/BASiCS",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "BASiCS: Bayesian analysis of single-cell sequencing data",
      doi = "10.1371/journal.pcbi.1004333",
      journal = "PLoS Computational Biology",
      date = "2015",
      peer_review = TRUE
    ),
    description = "BASiCS (Bayesian Analysis of Single-Cell Sequencing data) is an integrated Bayesian hierarchical model to perform statistical analyses of single-cell RNA sequencing datasets in the context of supervised experiments (where the groups of cells of interest are known a priori.",
    vignette = "http://47.254.148.113/software/Simsite/references/methods/17-basics/")

  list(BASiCS_method = BASiCS_method,
       BASiCS_parameters = BASiCS_parameters)
}
