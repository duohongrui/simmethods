#' Get Information of dyntoy
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' dyntoy_method_definition <- dyntoy_method_definition()
dyntoy_method_definition <- function(...){

  dyntoy_parameters <- parameter_sets(
    param_reference(
      id = "ref_data",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      description = "Reference dataset.",
      function_name = "dynwrap::infer_trajectory"
    ),
    param_others(
      id = "group.condition",
      type = "vector",
      process = "estimation",
      description = "Which groups or clusters that each cell belongs to",
      function_name = "dynwrap::infer_trajectory"
    ),
    param_numeric(
      id = "de.group",
      default = 0.1,
      description = "Ge probability.",
      function_name = "dynwrap::infer_trajectory"
    ),
    param_others(
      id = "nGenes",
      type = "integer",
      default = "nrow(ref_data)",
      description = "Total number of genes.",
      function_name = "dynwrap::infer_trajectory"
    ),
    param_others(
      id = "nCells",
      type = "integer",
      default = "ncol(ref_data)",
      description = "Total number of cells.",
      function_name = "dynwrap::infer_trajectory"
    ),
    param_integer(
      id = "seed",
      default = 687680L,
      function_name = "PROSSTT_python"
    ))

  dyntoy_method <- method_definition(
    method = "dyntoy",
    programming = "R",
    url = "https://github.com/dynverse/dyntoy",
    authors = authors_definition(
      first = "Robrecht",
      last = "Cannoodt",
      email = NULL,
      github = "https://github.com/dynverse/dyntoy",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = NULL,
      doi = NULL,
      journal = NULL,
      date = NULL,
      peer_review = NULL
    ),
    description = NULL,
    vignette = "http://47.254.148.113/software/Simsite/references/methods/32-dyntoy/")

  list(dyntoy_method = dyntoy_method,
       dyntoy_parameters = dyntoy_parameters)
}
