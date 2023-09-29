#' Get Information of VeloSim
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' VeloSim_method_definition <- VeloSim_method_definition()
VeloSim_method_definition <- function(...){

  VeloSim_parameters <- parameter_sets(
    param_reference(
      id = "ref_data",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      description = "Reference dataset.",
      function_name = "simutils::make_trees"
    ),
    param_others(
      id = "group.condition",
      type = "vector",
      process = "estimation",
      description = "Which groups or clusters that each cell belongs to",
      function_name = "simutils::make_trees"
    ),
    param_others(
      id = "ngenes",
      type = "integer",
      default = "nrow(ref_data)",
      description = "Total number of genes.",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_others(
      id = "ncells_total",
      type = "integer",
      default = "ncol(ref_data)",
      description = "Total number of cells.",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_others(
      id = "start_s",
      type = "matrix",
      description = "Initial spliced count, NULL if not specified.",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_others(
      id = "start_u",
      type = "matrix",
      description = "Initial unspliced count, NULL if not specified.",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_numeric(
      id = "evf_center",
      default = 1,
      description = "The value which evf mean is generated from",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_integer(
      id = "nevf",
      default = 20L,
      description = "Number of evfs.",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_others(
      id = "phyla",
      type = "tree",
      description = "The cell developmental tree if chosing 'discrete' or 'continuous' evf type. Can either be generated randomly (using pbtree(nclusters) function from phytools package) or read from newick format file using the ape package.",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_integer(
      id = "randseed",
      default = 687680L,
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_integer(
      id = "n_de_evf",
      default = 12L,
      description = "Number of differential evfs between populations.",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_character(
      id = "vary",
      description = "Which kinetic parameters should the differential evfs affect. Default is 's'. Can be 'kon', 'koff', 's', 'all', 'except_kon', 'except_koff', 'except_s'. Suggestions are 'all' or 's'.",
      alternatives = c("kon", "kon", "s", "all", "except_kon", "except_koff", "except_s"),
      default = "s",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_numeric(
      id = "Sigma",
      default = 0.1,
      description = "Parameter of the std of evf values within the same population.",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_numeric(
      id = "geffect_mean",
      default = 0,
      description = "The mean of the normal distribution where the non-zero gene effect sizes are sampled from.",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_numeric(
      id = "gene_effect_prob",
      default = 0.3,
      description = "The probability that the effect size is not 0.",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_numeric(
      id = "bimod",
      default = 0,
      description = "The amount of increased bimodality in the transcript distribution, 0 being not changed from the results calculated using evf and gene effects, and 1 being all genes are bimodal.",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_others(
      id = "param_realdata",
      type = "character",
      default = "zeisel.imputed",
      description = "Pick from zeisel.imputed or NULL; zeisel.imputed means using the distribution of kinetic parameters learned from the Zeisel 2015 dataset. This option is recommended.",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_numeric(
      id = "scale_s",
      default = 1,
      description = "A factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total).",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_numeric(
      id = "prop_hge",
      default = 0.015,
      description = "The proportion of very highly expressed genes.",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_numeric(
      id = "mean_hge",
      default = 5,
      description = "The parameter to amplify the gene-expression levels of the very highly expressed genes.",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_integer(
      id = "n_unstable",
      default = 0L,
      description = "The number of cells to be removed at the beginning of the simulation, sutract from ncells_total.",
      function_name = "VeloSim::SimulateVeloTree"
    ),
    param_Boolean(
      id = "plot",
      default = FALSE,
      description = "Plot the kinetic or not",
      function_name = "VeloSim::SimulateVeloTree"
    ))

  VeloSim_method <- method_definition(
    method = "VeloSim",
    programming = "R",
    url = "https://github.com/PeterZZQ/VeloSim",
    authors = authors_definition(
      first = "Ziqi",
      last = "Zhang",
      email = "xiuwei.zhang@gatech.edu",
      github = "https://github.com/PeterZZQ/VeloSim",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "VeloSim: Simulating single cell gene-expression and RNA velocity",
      doi = "10.1101/2021.01.11.426277",
      journal = "bioRxiv",
      date = 2021,
      peer_review = FALSE
    ),
    description = "This package simulates single cell RNA sequencing data during cell developmental process, the simulation includes unspliced RNA count, spliced RNA count, RNA velocity and cell developmental pseudo-time. The simulation can be used to benchmark computational methods including trajectory inference and RNA velocity inference.",
    vignette = "http://47.254.148.113/software/Simsite/references/methods/35-velosim/")

  list(VeloSim_method = VeloSim_method,
       VeloSim_parameters = VeloSim_parameters)
}
