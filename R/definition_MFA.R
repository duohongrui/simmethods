#' Get Information of MFA
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @importFrom splatter newMFAParams
#' @export
#'
#' @examples
#' MFA_method_definition <- MFA_method_definition()

MFA_method_definition <- function(...){
  MFA_parameters <- parameter_sets(
    param_reference(
      id = "input_refernece",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      description = "The reference data to go through estimation and evaluation process. Usually, no reference data is also OK.",
      function_name = "mfaEstimate"
    ),
    param_others(
      id = "MFAParams",
      type = "MFAParams",
      default = "splatter::newMFAParams()",
      process = "estimation",
      description = "Usually it is default by splatter::newMFAParams function. Users can change the parameters by splatter::setParam function.",
      function_name = "mfaEstimate"
    ),
    param_integer(
      id = "nGenes",
      default = 10000L,
      lower = 1L,
      description = "The number of genes to simulate.",
      function_name = "mfaSimulate"
    ),
    param_integer(
      id = "nCells",
      default = 100L,
      lower = 1L,
      description = "The number of cells to simulate.",
      function_name = "mfaSimulate"
    ),
    param_integer(
      id = "seed",
      force = TRUE,
      description = "Seed to use for generating random numbers.",
      function_name = "mfaSimulate"
    ),
    param_numeric(
      id = "trans.prop",
      default = 0,
      lower = 0,
      upper = 1,
      description = "Proportion of genes that show transient expression. These genes are briefly up or down-regulated before returning to their initial state.",
      function_name = "mfaSimulate"
    ),
    param_Boolean(
      id = "zero.neg",
      description = "Logical. Whether to set negative expression values to zero. This will zero-inflate the data.",
      default = TRUE,
      function_name = "mfaSimulate"
    ),
    param_Boolean(
      id = "dropout.present",
      description = "Logical. Whether to simulate dropout.",
      default = FALSE,
      function_name = "mfaSimulate"
    ),
    param_numeric(
      id = "dropout.lambda",
      default = 1,
      description = "Lambda parameter for the exponential dropout function.",
      function_name = "mfaSimulate"
    ),
    param_Boolean(
      id = "verbose",
      description = "Logical. Whether to print progress messages.",
      default = FALSE,
      function_name = "mfaSimulate"
    )
  )

  MFA_method <- method_definition(
    method = "MFA",
    programming = "R",
    url = "https://github.com/kieranrcampbell/mfa",
    authors = authors_definition(
      first = "Kieran R",
      last = "Campbell",
      email = "cyau@well.ox.ac.uk",
      github = "https://github.com/kieranrcampbell/mfa",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "Probabilistic modeling of bifurcations in single-cell gene expression data using a Bayesian mixture of factor analyzers",
      doi = "10.12688/wellcomeopenres.11087.1",
      journal = "Wellcome Open Research",
      date = "2017",
      peer_review = TRUE
    ),
    description = "Probabilistic inference of single-cell bifurcations.",
    vignette = "http://47.254.148.113/software/Simsite/references/methods/36-mfa/")

  list(MFA_method = MFA_method,
       MFA_parameters = MFA_parameters)
}
