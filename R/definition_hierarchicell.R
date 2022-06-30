#' Get Information of hierarchicell
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' hierarchicell_method_definition <- hierarchicell_method_definition()
#'
hierarchicell_method_definition <- function(...){

  muscat_parameters <- parameter_sets(
    param_reference(
      id = "x",
      type = "SingleCellExperiment",
      default = NULL,
      process = "estimation",
      force = TRUE,
      description = "A SingleCellExperiment.",
      function_name = c("prepSim", "simData")
    ),
    param_integer(
      id = "min_count",
      default = 1L,
      process = "estimation",
      description = "Used for filtering of genes; only genes with a count > min_count in >= min_cells will be retained.",
      function_name = "prepSim"
    ),
    param_integer(
      id = "min_cell",
      default = 10L,
      process = "estimation",
      description = "Used for filtering of genes; only genes with a count > min_count in >= min_cells will be retained.",
      function_name = "prepSim"
    ),
    param_integer(
      id = "min_genes",
      default = 100L,
      process = "estimation",
      description = "Used for filtering cells; only cells with a count > 0 in >= min_genes will be retained.",
      function_name = "prepSim"
    ),
    param_integer(
      id = "min_size",
      default = 100L,
      process = "estimation",
      description = "Used for filtering subpopulation-sample combinations; only instances with >= min_size cells will be retained. Specifying min_size = NULL skips this step.",
      function_name = "prepSim"
    ),
    param_others(
      id = "group_keep",
      type = "character",
      process = "estimation",
      description = "Character string; if nlevels(x$group_id) > 1, specifies which group of samples to keep (see details). The default NULL retains samples from levels(x$group_id)[1]; otherwise, if 'colData(x)$group_id' is not specified, all samples will be kept.",
      function_name = "prepSim"
    ),
    param_Boolean(
      id = "verbose",
      default = TRUE,
      process = "estimation",
      description = "Logical; should information on progress be reported?",
      function_name = "prepSim"
    ),
    param_others(
      id = "ng",
      type = "integer",
      process = "simulation",
      default = "nrow(x)",
      description = "Number of genes to simulate. Importantly, for the library sizes computed by prepSim (= exp(x$offset)) to make sense, the number of simulated genes should match with the number of genes in the reference. To simulate a reduced number of genes, e.g. for testing and development purposes, please set force = TRUE.",
      function_name = "simData"
    ),
    param_others(
      id = "nc",
      type = "integer",
      process = "simulation",
      default = "ncol(x)",
      description = "Number of cells to simulate.",
      function_name = "simData"
    ),
    param_others(
      id = "ns",
      type = "integer",
      process = "simulation",
      description = "Number of samples to simulate; defaults to as many as available in the reference to avoid duplicated reference samples. Specifically, the number of samples will be set to n = nlevels(x$sample_id) when dd = FALSE, n per group when dd, paired = TRUE, and floor(n/2) per group when dd = TRUE, paired = FALSE. When a larger number samples should be simulated, set force = TRUE.",
      function_name = "simData"
    ),
    param_others(
      id = "nk",
      type = "integer",
      process = "simulation",
      description = "Number of clusters to simulate; defaults to the number of available reference clusters.",
      function_name = "simData"
    ),
    param_others(
      id = "probs",
      type = "list",
      process = "simulation",
      description = "A list of length 3 containing probabilities of a cell belonging to each cluster, sample, and group, respectively. List elements must be NULL (equal probabilities) or numeric values in [0, 1] that sum to 1.",
      function_name = "simData"
    ),
    param_Boolean(
      id = "dd",
      default = TRUE,
      process = "simulation",
      description = "Whether or not to simulate differential distributions; if TRUE, two groups are simulated and ns corresponds to the number of samples per group, else one group with ns samples is simulated.",
      function_name = "simData"
    ),
    param_vector(
      id = "p_dd",
      process = "simulation",
      description = "Numeric vector of length 6. Specifies the probability of a gene being EE, EP, DE, DP, DM, or DB, respectively.",
      function_name = "simData"
    ),
    param_Boolean(
      id = "paired",
      default = FALSE,
      process = "simulation",
      description = "Logical specifying whether a paired design should be simulated (both groups use the same set of reference samples) or not (reference samples are drawn at random).",
      function_name = "simData"
    ),
    param_numeric(
      id = "p_ep",
      default = 0.5,
      lower = 0,
      upper = 1,
      description = "Numeric specifying the proportion of cells to be shifted to a different expression state in one group.",
      process = "simulation",
      function_name = "simData"
    ),
    param_numeric(
      id = "p_dp",
      default = 0.3,
      lower = 0,
      upper = 1,
      description = "Numeric specifying the proportion of cells to be shifted to a different expression state in one group.",
      process = "simulation",
      function_name = "simData"
    ),
    param_numeric(
      id = "p_dp",
      default = 0.5,
      lower = 0,
      upper = 1,
      description = "Numeric specifying the proportion of cells to be shifted to a different expression state in one group.",
      process = "simulation",
      function_name = "simData"
    ),
    param_numeric(
      id = "p_type",
      default = 0,
      lower = 0,
      upper = 1,
      description = "Numeric. Probability of EE/EP gene being a type-gene. If a gene is of class 'type' in a given cluster, a unique mean will be used for that gene in the respective cluster.",
      process = "simulation",
      function_name = "simData"
    ),
    param_numeric(
      id = "lfc",
      default = 2,
      description = "Numeric value to use as mean logFC (logarithm base 2) for DE, DP, DM, and DB type of genes.",
      process = "simulation",
      function_name = "simData"
    ),
    param_others(
      id = "rel_lfc",
      type = "numeric",
      default = NULL,
      description = "Numeric vector of relative logFCs for each cluster. Should be of length nlevels(x$cluster_id) with levels(x$cluster_id) as names. Defaults to factor of 1 for all clusters.",
      process = "simulation",
      function_name = "simData"
    ),
    param_others(
      id = "phylo_tree",
      type = "newick tree",
      default = "c(ifelse(is.null(phylo_tree), 0, 0.1), 3)",
      description = "Newick tree text representing cluster relations and their relative distance. An explanation of the syntax can be found here. The distance between the nodes, except for the original branch, will be translated in the number of shared genes between the clusters belonging to these nodes (this relation is controlled with phylo_pars). The distance between two clusters is defined as the sum of the branches lengths separating them.",
      function_name = "simData"
    ),
    param_Boolean(
      id = "force",
      default = FALSE,
      process = "simulation",
      description = "Logical specifying whether to force simulation when ng and/or ns don't match the number of available reference genes and samples, respectively.",
      function_name = "simData"
    )
  )

  hierarchicell_method <- method_definition(
    method = "hierarchicell",
    programming = "R",
    url = "https://github.com/HelenaLC/muscat",
    authors = authors_definition(
      first = "Helena L.",
      last = "Crowell",
      email = "helena.crowell@uzh.ch",
      github = "https://github.com/HelenaLC/muscat",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data",
      doi = "10.1038/s41467-020-19894-4",
      journal = "Nature Communications",
      date = "2020",
      peer_review = TRUE
    ),
    description = "Multi-sample multi-group scRNA-seq analysis tools.")

  list(hierarchicell_method = hierarchicell_method,
       hierarchicell_parameters = hierarchicell_parameters)
}
