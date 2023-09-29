#' Get Information of POWSC
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' POWSC_method_definition <- POWSC_method_definition()
#'
POWSC_method_definition <- function(...){

  POWSC_parameters <- parameter_sets(
    param_reference(
      id = "sce",
      type = c("matrix", "SingleCellExperiment"),
      default = NULL,
      process = "estimation",
      description = "SingleCellExperiment object with assays(sce)[[1]] is the count matrix or input directly.",
      function_name = "Est2Phase"
    ),
    param_numeric(
      id = "low.prob",
      default = 0.99,
      lower = 0,
      upper = 1,
      process = "estimation",
      description = "lower bound probability for phase I.",
      function_name = "Est2Phase"
    ),
    param_integer(
      id = "n",
      default = 100L,
      lower = 0,
      border = FALSE,
      process = "simulation",
      description = "The number of total cells for two groups.",
      function_name = "Simulate2SCE"
    ),
    param_numeric(
      id = "perDE",
      default = 0.05,
      lower = 0,
      upper = 1,
      process = "simulation",
      description = "Percentage of DE genes.",
      function_name = "Simulate2SCE"
    ),
    param_others(
      id = "estParas1",
      type = "list",
      force = TRUE,
      process = "simulation",
      description = "The set of parameters corresponding to cell type I.",
      function_name = "Simulate2SCE"
    ),
    param_others(
      id = "estParas2",
      type = "list",
      force = TRUE,
      process = "simulation",
      description = "The set of parameters corresponding to cell type I.",
      function_name = "Simulate2SCE"
    )
  )

  POWSC_method <- method_definition(
    method = "POWSC",
    programming = "R",
    url = "http://www.bioconductor.org/packages/release/bioc/html/POWSC.html",
    authors = authors_definition(
      first = "Kenong",
      last = "Su",
      email = NULL,
      github = "https://github.com/suke18/POWSC",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "Simulation, power evaluation and sample size recommendation for single-cell RNA-seq",
      doi = "10.1093/bioinformatics/btaa607",
      journal = "Bioinformatics",
      date = "2020",
      peer_review = TRUE
    ),
    description = "POWSC, a simulation-based method, is designed to provide power evaluation and sample size recommendation for single-cell RNA sequencing DE analysis",
    vignette = "http://47.254.148.113/software/Simsite/references/methods/9-powsc/")

  list(POWSC_method = POWSC_method,
       POWSC_parameters = POWSC_parameters)
}
