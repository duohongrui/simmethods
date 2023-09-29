#' Get Information of dyngen
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' dyngen_method_definition <- dyngen_method_definition()
#'
dyngen_method_definition <- function(...){

  dyngen_parameters <- parameter_sets(
    param_reference(
      id = "input_refernece",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      description = "The reference data to go through estimation and evaluation process. Usually, no reference data is also OK.",
      function_name = "dynwrap::infer_trajectory"
    ),
    param_others(
      id = "backbone",
      type = "list",
      force = TRUE,
      description = "The gene module configuration that determines the type of dynamic process being simulated. See list_backbones() for a full list of different backbones available in this package.",
      function_name = "initialise_model"
    ),
    param_integer(
      id = "num_cells",
      default = 1000L,
      lower = 0L,
      description = "The number of cells to sample.",
      function_name = "initialise_model"
    ),
    param_others(
      id = "num_tfs",
      type = "integer",
      default = "nrow(backbone$module_info)",
      description = "The number of transcription factors (TFs) to generate. TFs are the main drivers of the changes that occur in a cell. TFs are regulated only by other TFs.",
      function_name = "initialise_model"
    ),
    param_integer(
      id = "num_targets",
      default = 100L,
      lower = 0,
      description = "The number of target genes to generate. Target genes are regulated by TFs and sometimes by other target genes.",
      function_name = "initialise_model"
    ),
    param_integer(
      id = "num_hks",
      default = 50L,
      lower = 0,
      description = "The number of housekeeping genes (HKs) to generate. HKs are typically highly expressed, and are not regulated by the TFs or targets.",
      function_name = "initialise_model"
    ),
    param_character(
      id = "distance_metric",
      default = "euclidean",
      alternatives = c("pearson", "spearman", "cosine", "euclidean", "manhattan"),
      description = "The distance metric to be used to calculate the distance between cells. See dynutils::calculate_distance() for a list of possible distance metrics.",
      function_name = "initialise_model"
    ),
    param_others(
      id = "tf_network_params",
      default = "tf_network_default()",
      type = "list",
      description = "Settings for generating the TF network with generate_tf_network(), see tf_network_default().",
      function_name = "initialise_model"
    ),
    param_others(
      id = "feature_network_params",
      type = "list",
      default = "feature_network_default()",
      description = "Settings for generating the feature network with generate_feature_network(), see feature_network_default().",
      function_name = "initialise_model"
    ),
    param_others(
      id = "kinetics_params",
      type = "list",
      default = "kinetics_default()",
      description = "Settings for determining the kinetics of the feature network with generate_kinetics(), see kinetics_default().",
      function_name = "initialise_model"
    ),
    param_others(
      id = "gold_standard_params",
      type = "list",
      default = "gold_standard_default()",
      description = "Settings pertaining simulating the gold standard with generate_gold_standard(), see gold_standard_default().",
      function_name = "initialise_model"
    ),
    param_others(
      id = "simulation_params",
      type = "list",
      default = "simulation_default()",
      description = "Settings pertaining the simulation itself with generate_cells(), see simulation_default().",
      function_name = "initialise_model"
    ),
    param_others(
      id = "experiment_params",
      type = "list",
      default = "experiment_snapshot()",
      description = "Settings related to how the experiment is simulated with generate_experiment(), see experiment_snapshot() or experiment_synchronised().",
      function_name = "initialise_model"
    ),
    param_others(
      id = "download_cache_dir",
      type = "character",
      default = "tools::R_user_dir('dyngen', 'data')",
      description = "If not NULL, temporary downloaded files will be cached in this directory.",
      function_name = "initialise_model"
    ),
    param_integer(
      id = "num_cores",
      default = 1L,
      lower = 1L,
      description = "Parallellisation parameter for various steps in the pipeline.",
      function_name = "initialise_model"
    ),
    param_others(
      id = "id",
      type = "character",
      description = "An identifier for the model.",
      function_name = "initialise_model"
    )
  )

  dyngen_method <- method_definition(
    method = "dyngen",
    programming = "R",
    url = "https://cran.r-project.org/web/packages/dyngen/index.html",
    authors = authors_definition(
      first = "Robrecht",
      last = "Cannoodt",
      email = "robrecht@cannoodt.dev",
      github = "https://github.com/dynverse/dyngen",
      orcid = "0000-0003-3641-729X"
    ),
    manuscript = manuscript_definition(
      title = "Spearheading future omics analyses using dyngen, a multi-modal simulator of single cells",
      doi = "10.1038/s41467-021-24152-2",
      journal = "Nature Communications",
      date = "2021",
      peer_review = TRUE
    ),
    description = "Simulating single-cell data using gene regulatory networks.",
    vignette = "http://47.254.148.113/software/Simsite/references/methods/33-dyngen/")

  list(dyngen_method = dyngen_method,
       dyngen_parameters = dyngen_parameters)
}
