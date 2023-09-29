#' Get Information of dropsim
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' dropsim_method_definition <- dropsim_method_definition()
#'
dropsim_method_definition <- function(...){

  dropsim_parameters <- parameter_sets(
    param_reference(
      id = "data",
      type = "matrix",
      default = NULL,
      process = "estimation",
      force = TRUE,
      description = "Matrix; gene by cell dge matrix of counts.",
      function_name = "fit_parameters "
    ),
    param_others(
      id = "parameters",
      type = "dropsim_parameters",
      force = TRUE,
      description = "Object of class dropsim containing parameters.",
      function_name = "simulateDGE"
    ),
    param_Boolean(
      id = "sparse",
      default = TRUE,
      description = "Logical; if true the counts are returned as a sparse matrix",
      function_name = "simulateDGE"
    ),
    param_others(
      id = "cell_prefix",
      type = "character",
      default = "cell",
      description = "Character; the default group name for cells",
      function_name = "simulateDGE"
    ),
    param_Boolean(
      id = "dge",
      default = TRUE,
      description = "Logical; if false counts a returned as a matrix cells by genes rather than genes by cells (If sparse=FALSE, setting dge=FALSE can save time on larger datasets.)",
      function_name = "simulateDGE"
    ),
    param_others(
      id = "seed",
      type = "integer",
      description = "Integer; seed for random number generation, set this for repoduciable simulations.",
      function_name = "simulateDGE"
    )
  )

  dropsim_method <- method_definition(
    method = "dropsim",
    programming = "R",
    url = "https://github.com/marchinilab/dropsim",
    authors = authors_definition(
      first = "Matt",
      last = "mkerin",
      email = "m4tt.kerin@gmail.com",
      github = "https://github.com/marchinilab/dropsim",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = NULL,
      doi = NULL,
      journal = NULL,
      date = NULL,
      peer_review = NULL
    ),
    description = "R Package for Single Cell RNAseq Synthetic Data Simulation.",
    vignette = "http://47.254.148.113/software/Simsite/references/methods/27-dropsim/")

  list(dropsim_method = dropsim_method,
       dropsim_parameters = dropsim_parameters)
}
