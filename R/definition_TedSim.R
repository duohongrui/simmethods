#' Get Information of TedSim
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' TedSim_method_definition <- TedSim_method_definition()
TedSim_method_definition <- function(...){

  TedSim_parameters <- parameter_sets(
    param_reference(
      id = "input_refernece",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      force = T,
      description = "A matrix with cells in columns and genes in rows.",
      function_name = "simmethods::TedSim_est"
    ),
    param_others(
      id = "group.condition",
      type = "vector",
      process = "estimation",
      description = "Which groups or clusters that each cell belongs to",
      function_name = "simmethods::TedSim_est"
    ),
    param_numeric(
      id = "de.group",
      default = 0.1,
      description = "Ge probability.",
      function_name = "CIF2Truecounts"
    ),
    param_others(
      id = "cif_res",
      type = "list",
      description = "The CIFs simulated.",
      function_name = "CIF2Truecounts"
    ),
    param_others(
      id = "scale_s",
      type = "list",
      description = "Transcription rate scaler, or a vector to specify cell-type specific scale_s.",
      function_name = "CIF2Truecounts"
    ),
    param_numeric(
      id = "mean_hge",
      default = 5,
      lower = 0,
      border = FALSE,
      description = "The mean of hge, default is 5.",
      function_name = "CIF2Truecounts"
    ),
    param_numeric(
      id = "prob_hge",
      default = 0.015,
      lower = 0,
      border = FALSE,
      description = "The probability of hge.",
      function_name = "CIF2Truecounts"
    ),
    param_others(
      id = "nGenes",
      type = "integer",
      default = "nrow(ref_data)",
      description = "Total number of genes.",
      function_name = "CIF2Truecounts"
    ),
    param_integer(
      id = "seed",
      default = 687680L,
      function_name = "PROSSTT_python"
    ))

  TedSim_method <- method_definition(
    method = "TedSim",
    programming = "R",
    url = "https://github.com/Galaxeee/TedSim",
    authors = authors_definition(
      first = "Xinhai",
      last = "Pan",
      email = "xiuwei.zhang@gatech.edu",
      github = "https://github.com/Galaxeee/TedSim",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "TedSim: temporal dynamics simulation of single-cell RNA sequencing data and cell division history",
      doi = "10.1093/nar/gkac235",
      journal = "Nucleic Acids Research",
      date = "2022",
      peer_review = TRUE
    ),
    description = "Temporal dynamics simulation of single-cell RNA sequencing data and cell division history",
    vignette = "http://47.254.148.113/software/Simsite/references/methods/30-tedsim/")

  list(TedSim_method = TedSim_method,
       TedSim_parameters = TedSim_parameters)
}
