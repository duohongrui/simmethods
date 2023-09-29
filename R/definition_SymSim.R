#' Get Information of SymSim
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' SymSim_method_definition <- SymSim_method_definition()
SymSim_method_definition <- function(...){

  SymSim_parameters <- parameter_sets(
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
    param_integer(
      id = "min_popsize",
      default = 1L,
      description = "The number of cells in the smallest population.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_others(
      id = "ngenes",
      type = "integer",
      default = "nrow(ref_data)",
      description = "Total number of genes.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_others(
      id = "ncells_total",
      type = "integer",
      default = "ncol(ref_data)",
      description = "Total number of cells.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_integer(
      id = "i_minpop",
      default = 1L,
      description = "Specifies which population has the smallest size.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_numeric(
      id = "evf_center",
      default = 1,
      description = "The value which evf mean is generated from",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_character(
      id = "evf_type",
      description = "string that is one of the following: 'one.population','discrete','continuous'",
      alternatives = c("one.population", "discrete", "continuous"),
      default = "continuous",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_integer(
      id = "nevf",
      default = 10L,
      description = "Number of evfs.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_others(
      id = "phyla",
      type = "tree",
      description = "The cell developmental tree if chosing 'discrete' or 'continuous' evf type. Can either be generated randomly (using pbtree(nclusters) function from phytools package) or read from newick format file using the ape package.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_integer(
      id = "randseed",
      default = 687680L,
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_integer(
      id = "n_de_evf",
      default = 0L,
      description = "Number of differential evfs between populations.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_character(
      id = "vary",
      description = "Which kinetic parameters should the differential evfs affect. Default is 's'. Can be 'kon', 'koff', 's', 'all', 'except_kon', 'except_koff', 'except_s'. Suggestions are 'all' or 's'.",
      alternatives = c("kon", "kon", "s", "all", "except_kon", "except_koff", "except_s"),
      default = "s",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_numeric(
      id = "Sigma",
      default = 0.4,
      description = "Parameter of the std of evf values within the same population.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_numeric(
      id = "geffect_mean",
      default = 0,
      description = "The mean of the normal distribution where the non-zero gene effect sizes are sampled from.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_numeric(
      id = "gene_effect_prob",
      default = 0.3,
      description = "The probability that the effect size is not 0.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_numeric(
      id = "bimod",
      default = 0,
      description = "The amount of increased bimodality in the transcript distribution, 0 being not changed from the results calculated using evf and gene effects, and 1 being all genes are bimodal.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_others(
      id = "param_realdata",
      type = "character",
      default = "zeisel.imputed",
      description = "Pick from zeisel.imputed or NULL; zeisel.imputed means using the distribution of kinetic parameters learned from the Zeisel 2015 dataset. This option is recommended.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_numeric(
      id = "scale_s",
      default = 1,
      description = "A factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total).",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_Boolean(
      id = "impulse",
      default = FALSE,
      description = "Use the impulse function when generating continuous population or not. Default is F.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_numeric(
      id = "gene_module_prop",
      default = 0,
      description = "Proportion of genes which are in co-expressed gene module.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_numeric(
      id = "prop_hge",
      default = 0.015,
      description = "The proportion of very highly expressed genes.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_numeric(
      id = "mean_hge",
      default = 5,
      description = "The parameter to amplify the gene-expression levels of the very highly expressed genes.",
      function_name = "SymSim::SimulateTrueCounts"
    ),
    param_numeric(
      id = "gene_effects_sd",
      default = 1,
      description = "The standard deviation of the normal distribution where the non-zero gene effect sizes are sampled from.",
      function_name = "SymSim::SimulateTrueCounts"
    ))

  SymSim_method <- method_definition(
    method = "SymSim",
    programming = "R",
    url = "https://github.com/YosefLab/SymSim",
    authors = authors_definition(
      first = "Xiuwei",
      last = "Zhang",
      email = "niryosef@berkeley.edu",
      github = "https://github.com/YosefLab/SymSim",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "Simulating multiple faceted variability in single cell RNA sequencing",
      doi = "10.1038/s41467-019-10500-w",
      journal = "Nature Communications",
      date = 2019,
      peer_review = TRUE
    ),
    description = NULL,
    vignette = "http://47.254.148.113/software/Simsite/references/methods/34-symsim/")

  list(SymSim_method = SymSim_method,
       SymSim_parameters = SymSim_parameters)
}
