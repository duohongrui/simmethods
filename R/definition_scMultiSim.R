#' Get Information of scMultiSim
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' scMultiSim_method_definition <- scMultiSim_method_definition()
#'
scMultiSim_method_definition <- function(...){

  scMultiSim_parameters <- parameter_sets(
    param_others(
      id = "simsrt",
      type = "SRTsim",
      default = NULL,
      process = "estimation",
      force = TRUE,
      description = "A SRTsim object",
      function_name = "srtsim_fit"
    ),
    param_character(
      id = "marginal",
      default = "auto_choose",
      alternatives = c("auto_choose", "zinb", "nb", "poisson", "zip"),
      process = "estimation",
      description = "Specification of the types of marginal distribution.Default value is 'auto_choose' which chooses between ZINB, NB, ZIP, and Poisson by a likelihood ratio test (lrt),AIC and whether there is underdispersion.'zinb' will fit the ZINB model. If there is underdispersion, it will fit the Poisson model. If there is no zero at all or an error occurs, it will fit an NB model instead.'nb' fits the NB model and chooses between NB and Poisson depending on whether there is underdispersion. 'poisson' simply fits the Poisson model.'zip' fits the ZIP model and chooses between ZIP and Poisson by a likelihood ratio test.",
      function_name = "srtsim_fit"
    ),
    param_character(
      id = "sim_scheme",
      default = "domain",
      alternatives = c("domain", "tissue"),
      process = "estimation",
      description = "A character string specifying simulation scheme. 'tissue' stands for tissue-based simulation; 'domain' stands for domain-specific simulation. Default is domain.",
      function_name = "srtsim_fit"
    ),
    param_numeric(
      id = "min_nonzero_num",
      default = 2,
      process = "estimation",
      description = "The minimum number of non-zero values required for a gene to be fitted. Default is 2.",
      function_name = "srtsim_fit"
    ),
    param_numeric(
      id = "maxiter",
      default = 500,
      process = "estimation",
      description = "The number of iterations for the model-fitting. Default is 500.",
      function_name = "srtsim_fit"
    ),
    param_Boolean(
      id = "write.noise.model",
      process = "estimation",
      default = TRUE,
      description = "When write.noise.model=TRUE outputs two tab-delimited files containing the dropout effects and noise model parameters; this allows users to apply the noise generation on a seperate high compute node. The root file name is set by file option.",
      function_name = "srtsim_count"
    ),
    param_others(
      id = "simsrt",
      type = "SRTsim",
      default = NULL,
      process = "simulation",
      force = TRUE,
      description = "A object with estimated parameters from fitting step.",
      function_name = "srtsim_count"
    ),
    param_others(
      id = "total_count_new",
      type = c("numeric", "NULL"),
      default = NULL,
      process = "simulation",
      description = "The (expected) total number of reads or UMIs in the simulated count matrix.",
      function_name = "srtsim_count"
    ),
    param_others(
      id = "total_count_old",
      type = c("numeric", "NULL"),
      default = NULL,
      process = "simulation",
      description = "The total number of reads or UMIs in the original count matrix.",
      function_name = "srtsim_count"
    ),
    param_others(
      id = "rrr",
      type = c("numeric", "NULL"),
      default = NULL,
      process = "simulation",
      description = "The ratio applies to the gene-specific mean estimate, used for the fixing average sequencing depth simulation. Default is null. Its specification will override the specification of total_count_new and total_count_old.",
      function_name = "srtsim_count"
    ),
    param_numeric(
      id = "nn_num",
      default = 5,
      process = "simulation",
      description = "A integer of nearest neighbors, default is 5.",
      function_name = "srtsim_count"
    ),
    param_integer(
      id = "numCores",
      default = 1L,
      description = "The number of cores to use.",
      process = "simulation",
      function_name = "srtsim_count"
    ),
    param_character(
      id = "nn_func",
      default = "mean",
      alternatives = c("mean", "median", "ransam"),
      process = "simulation",
      description = "A character string specifying how the psedo-count to be generated. options include 'mean','median' and 'ransam'.",
      function_name = "srtsim_count"
    ),
    param_Boolean(
      id = "verbose",
      default = FALSE,
      process = "simulation",
      description = "Whether to show running information for srtsim_count.",
      function_name = "srtsim_count"
    )
  )

  scMultiSim_method <- method_definition(
    method = "scMultiSim",
    programming = "R",
    url = "https://github.com/ZhangLabGT/scMultiSim",
    authors = authors_definition(
      first = "Hechen",
      last = "Li",
      email = NULL,
      github = "https://github.com/ZhangLabGT/scMultiSim",
      orcid = "0000-0002-1455-0041"
    ),
    manuscript = manuscript_definition(
      title = "SRTsim: spatial pattern preserving simulations for spatially resolved transcriptomics",
      doi = "10.1186/s13059-023-02879-z",
      journal = "Genome Bilogy",
      date = "2023",
      peer_review = TRUE
    ),
    description = "An independent, reproducible, and flexible Spatially Resolved Transcriptomics (SRT) simulation framework that can be used to facilitate the development of SRT analytical methods for a wide variety of SRT-specific analyses.")

  list(scMultiSim_method = scMultiSim_method,
       scMultiSim_parameters = scMultiSim_parameters)
}
