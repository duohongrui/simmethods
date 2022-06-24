#' Get Information of SPARSim
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' SPARSim_method_definition <- SPARSim_method_definition()
#'
SPARSim_method_definition <- function(...){

  SPARSim_parameters <- parameter_sets(
    param_reference(
      id = "raw_data",
      type = "matrix",
      default = NULL,
      process = "estimation",
      force = TRUE,
      description = "Count matrix (gene on rows, samples on columns) containing raw count data.",
      function_name = "SPARSim_estimate_parameter_from_data"
    ),
    param_reference(
      id = "norm_data",
      type = "matrix",
      default = NULL,
      process = "estimation",
      force = TRUE,
      description = "Count matrix (gene on rows, samples on columns) containing normalized count data.",
      function_name = "SPARSim_estimate_parameter_from_data"
    ),
    param_others(
      id = "dataset_parameter",
      type = "list",
      default = NULL,
      process = "simulation",
      force = TRUE,
      description = "List containing, the intensity, variability and lib sizes of each experimental condition. It is the return value of 'estimate_parameter_from_data' or could be created by the users",
      function_name = "SPARSim_simulation"
    ),
    param_others(
      id = "batch_parameter",
      type = "list",
      default = NULL,
      process = "simulation",
      description = "Bacth effect simulation parameter created by SPARSim_create_batch_parameter function",
      function_name = "SPARSim_simulation"
    ),
    param_others(
      id = "spikein_parameter",
      type = "list",
      default = NULL,
      process = "simulation",
      description = "Spike-in simulation parameter created by SPARSim_create_spikein_parameter function",
      function_name = "SPARSim_simulation"
    ),
    param_Boolean(
      id = "output_sim_param_matrices",
      default = FALSE,
      process = "simulation",
      description = "Boolean flag. If TRUE, the function will output two additional matrices, called abundance_matrix and variability_matrix, containing the gene intensities and gene variabilities used as simulation input. (Default: FALSE)",
      function_name = "SPARSim_simulation"
    ),
    param_Boolean(
      id = "output_batch_matrix",
      default = FALSE,
      process = "simulation",
      description = "Boolean flag. If TRUE, the function will output an additional matrix, called batch_factors_matrix, containing the multiplicative factors used in batch effect simulation. (Default: FALSE)",
      function_name = "SPARSim_simulation"
    )
  )

  SPARSim_method <- method_definition(
    method = "SPARSim",
    programming = "R",
    url = "https://gitlab.com/sysbiobig/sparsim",
    authors = authors_definition(
      first = "Giacomo",
      last = "Baruzzo",
      email = "barbara.dicamillo@unipd.it",
      github = NULL,
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "SPARSim single cell: a count data simulator for scRNA-seq data",
      doi = "10.1093/bioinformatics/btz752",
      journal = "Bioinformatics",
      date = "2019",
      peer_review = TRUE
    ),
    description = "SPARSim is an R tool for the simulation of single cell RNA-seq (scRNA-seq) count table.")

  list(SPARSim_method = SPARSim_method,
       SPARSim_parameters = SPARSim_parameters)
}
