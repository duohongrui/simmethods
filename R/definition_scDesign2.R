#' Get Information of scDesign2
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' scDesign2_method_definition <- scDesign2_method_definition()
#'
scDesign2_method_definition <- function(...){

  scDesign2_parameters <- parameter_sets(
    param_reference(
      id = "data_mat",
      type = "matrix",
      default = NULL,
      process = "estimation",
      force = TRUE,
      description = "A matrix of shape p by n that contains count values. Each of its column names should be the cell type names of that cell. Its column names should also match cell_type_sel.",
      function_name = "fit_model_scDesign2"
    ),
    param_vector(
      id = "cell_type_sel",
      process = "estimation",
      description = "A character vector that contains the selected cell types for which a model will be fitted.",
      function_name = "fit_model_scDesign2"
    ),
    param_character(
      id = "sim_method",
      alternatives = c("copula", "ind"),
      process = "estimation",
      description = "Specification of the type of model. Default value is 'copula', which selects the copula model. 'ind' will select the (w/o copula) model.",
      function_name = c("fit_model_scDesign2", "simulate_count_scDesign2")
    ),
    param_character(
      id = "marginal",
      alternatives = c("auto_choose", "zinb", "nb", "poisson"),
      force = TRUE,
      process = "estimation",
      description = "Specification of the types of marginal distribution. Default value is 'auto_choose' which chooses between ZINB, NB, ZIP and Poisson by a likelihood ratio test (lrt) and whether there is underdispersion. 'zinb' will fit the ZINB model. If there is underdispersion, it will choose between ZIP and Poisson by a lrt. Otherwise, it will try to fit the ZINB model. If in this case, there is no zero at all or an error occurs, it will fit an NB model instead. 'nb' fits the NB model that chooses between NB and Poisson depending on whether there is underdispersion. 'poisson' simply fits the Poisson model.",
      function_name = "fit_model_scDesign2"
    ),
    param_Boolean(
      id = "jitter",
      default = TRUE,
      process = "estimation",
      description = "Logical, whether a random projection should be performed in the distributional transform.",
      function_name = "fit_model_scDesign2"
    ),
    param_numeric(
      id = "zp_cutoff",
      process = "estimation",
      default = 0.8,
      lower = 0,
      upper = 1,
      description = "The maximum propotion of zero allowed for a gene to be included in the joint copula model.",
      function_name = "fit_model_scDesign2"
    ),
    param_integer(
      id = "ncores",
      process = "estimation",
      default = 1L,
      lower = 0L,
      description = "A numeric value that indicates the number of parallel cores for model fitting. One core for each cell type.",
      function_name = "fit_model_scDesign2"
    ),
    param_others(
      id = "model_params",
      force = TRUE,
      type = "list",
      description = "A list with the same length as cell_type_prop that contains the fitted model as each of its element (can be either the copula model or the (w/o copula) model).",
      function_name = "simulate_count_scDesign2"
    ),
    param_others(
      id = "n_cell_new",
      type = "integer",
      force = TRUE,
      process = "simulation",
      description = "The total number of cells in the simulated count matrix.",
      function_name = "simulate_count_scDesign2"
    ),
    param_numeric(
      id = "cell_type_prop",
      default = 1,
      lower = 0,
      upper = 1,
      process = "simulation",
      description = "The cell type proportion in the simulated count matrix.",
      function_name = "simulate_count_scDesign2"
    ),
    param_others(
      id = "total_count_new",
      type = "numeric",
      default = NULL,
      description = "The (expected) total number of reads or UMIs in the simulated count matrix.",
      process = "simulation",
      function_name = "simulate_count_scDesign2"
    ),
    param_others(
      id = "total_count_old",
      type = "numeric",
      default = NULL,
      description = "The total number of reads or UMIs in the original count matrix where model_params was fitted.",
      process = "simulation",
      function_name = "simulate_count_scDesign2"
    ),
    param_others(
      id = "n_cell_old",
      type = "numeric",
      default = NULL,
      description = "The The total number of cells in the original count matrix where model_params was fitted.",
      process = "simulation",
      function_name = "simulate_count_scDesign2"
    ),
    param_character(
      id = "reseq_method",
      alternatives = c("mean_scale", "multinomial"),
      description = "Specification of how the new count matrix should be derived under the new sequencing depth. Default is 'mean_scale', which scales the original parameters and then simulate new data. 'multinomial' will do a resampling. It ensures that the simulated count matrix has the exact total number of reads as specified in total_count_new.",
      function_name = "simulate_count_scDesign2"
    ),
    param_Boolean(
      id = "cell_sample",
      default = FALSE,
      description = "Logical, whether cells for each cell type should be sampled from a multinomial distribution or follows the exact same proportion as specified in cell_type_prop.",
      function_name = "simulate_count_scDesign2"
    )
  )

  scDesign2_method <- method_definition(
    method = "scDesign2",
    programming = "R",
    url = "https://github.com/JSB-UCLA/scDesign2",
    authors = authors_definition(
      first = "Tianyi",
      last = "Sun",
      email = NULL,
      github = "https://github.com/JSB-UCLA/scDesign2",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "scDesign2: a transparent simulator that generates high-fidelity single-cell gene expression count data with gene correlations captured",
      doi = "10.1186/s13059-021-02367-2",
      journal = "Genome Biology",
      date = "2021",
      peer_review = TRUE
    ),
    description = "An interpretable simulator that generates realistic single-cell gene expression count data with gene correlations recapitulated.")

  list(scDesign2_method = scDesign2_method,
       scDesign2_parameters = scDesign2_parameters)
}
