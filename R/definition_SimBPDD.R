#' Get Information of SimBPDD
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' muscat_method_definition <- muscat_method_definition()
#'
SimBPDD_method_definition <- function(...){

  SimBPDD_parameters <- parameter_sets(
    param_others(
      id = "DATA",
      type = "matrix",
      process = "simulation",
      description = "Matrix of single-cell RNA-sequencing expression data with genes in rows and samples (cells) in columns.",
      function_name = "bp.sim.DD"
    ),
    param_vector(
      id = "N1",
      type = "integer",
      process = "simulation",
      description = "size of the samples drawn from the control BP models.",
      function_name = "bp.sim.DD"
    ),
    param_vector(
      id = "N2",
      type = "integer",
      process = "simulation",
      description = "size of the samples drawn from the manipulated BP models.",
      function_name = "bp.sim.DD"
    ),
    param_character(
      id = "case",
      default = "DBeta",
      alternatives = c("DBeta", "DLambda", "DAlpha", "DAlphaBeta", "DPZ"),
      process = "simulation",
      description = "A cell by feature matrix of the multivariate quantile.",
      function_name = "bp.sim.DD"
    ),
    param_numeric(
      id = "degree",
      default = 1/2,
      lower = 0,
      border = FALSE,
      process = "simulation",
      description = "Parameter to set the degree of the created difference (low to strong), see Schefzik (2021) for details and the choice of a range of possible values.",
      function_name = "bp.sim.DD"
    ),
    param_vector(
      id = "seedex",
      description = "Seed used for sampling from the fitted BP models to ensure reproducibility.",
      process = "simulation",
      function_name = "bp.sim.DD"
    )
  )

  SimBPDD_method <- method_definition(
    method = "SimBPDD",
    programming = "R",
    url = "https://github.com/RomanSchefzik/SimBPDD",
    authors = authors_definition(
      first = "Roman",
      last = "Schefzik",
      email = "Roman.Schefzik@medma.uni-heidelberg.de",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "Simulating differential distributions in Beta-Poisson models, in particular for single-cell RNA sequencing data",
      doi = "10.33039/ami.2021.03.003",
      journal = "Annales Mathematicae et Informaticae",
      date = "2021",
      peer_review = TRUE
    ),
    description = "Simulating differential distributions for Beta-Poisson models, in particular for single-cell RNA-sequencing (scRNA-seq) data.")

  list(SimBPDD_method = SimBPDD_method,
       SimBPDD_parameters = SimBPDD_parameters)
}
