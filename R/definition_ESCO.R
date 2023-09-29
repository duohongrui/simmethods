#' Get Information of ESCO
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' ESCO_method_definition <- ESCO_method_definition()
#'
ESCO_method_definition <- function(...){

  ESCO_parameters <- parameter_sets(
    param_reference(
      id = "counts",
      type = "matrix",
      default = NULL,
      process = "estimation",
      force = TRUE,
      description = "Either a counts matrix or a SingleCellExperiment object containing count data to estimate parameters from.",
      function_name = "escoEstimate"
    ),
    param_others(
      id = "dirname",
      type = "character",
      default = NULL,
      process = "estimation",
      description = "A tring of directory name to indicate where to save the results.",
      function_name = "escoEstimate"
    ),
    param_Boolean(
      id = "group",
      default = FALSE,
      process = "estimation",
      description = "Whether the data is believed to be of discrete cell groups or not, if yes, the corresponding cellinfo indicating the cell group labels need to be input as well.",
      function_name = "escoEstimate"
    ),
    param_vector(
      id = "cellinfo",
      process = "estimation",
      description = "A vector of length n, where n is the number of cells. Each entries is the group identity of a cell.",
      function_name = "escoEstimate"
    ),
    param_others(
      id = "params",
      type = "newescoParams",
      default = "newescoParams()",
      process = c("estimation", "simulation"),
      description = "escoParams object to store estimated values in.",
      function_name = c("escoEstimate", "escoSimulate")
    ),
    param_character(
      id = "type",
      default = "single",
      alternatives = c("single", "group", "tree", "traj"),
      process = "simulation",
      description = "Which type of heterogenounity to use. Options are : 'single' which produces a single population; 'group' which produces distinct groups; 'tree' which produces distinct groups but admits a tree structure; 'traj' which produces distinct groups but admits a smooth trajectory structure.",
      function_name = "escoSimulate"
    ),
    param_Boolean(
      id = "paths",
      default = FALSE,
      description = "Whether to simulation trajectory datasets",
      function_name = "escoSimulate"
    ),
    param_Boolean(
      id = "tree",
      default = FALSE,
      description = "Whether to simulation datasets from tree format files",
      function_name = "escoSimulate"
    ),
    param_Boolean(
      id = "verbose",
      default = TRUE,
      process = "simulation",
      description = "Logical. Whether to print progress messages.",
      function_name = "escoSimulate"
    ),
    param_integer(
      id = "numCores",
      default = 2L,
      lower = 1L,
      description = "The number of cores used for parallelization, default is 2. If set as NULL, then all the detected cores are used.",
      function_name = "escoSimulate"
    ),
    param_integer(
      id = "nGenes",
      default = 10000L,
      lower = 1L,
      description = "The number of genes to simulate.",
      function_name = "escoSimulate"
    ),
    param_integer(
      id = "nCells",
      default = 100L,
      lower = 1L,
      description = "The number of cells to simulate.",
      function_name = "escoSimulate"
    ),
    param_integer(
      id = "seed",
      force = TRUE,
      description = "Seed to use for generating random numbers.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "mean.shape",
      default = 0.6,
      lower = 0,
      border = FALSE,
      description = "Shape parameter for the mean gamma distribution.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "mean.rate",
      default = 0.3,
      lower = 0,
      border = FALSE,
      description = "Rate parameter for the mean gamma distribution.",
      function_name = "escoSimulate"
    ),
    param_others(
      id = "mean.dens",
      type = "list",
      description = "A density object.",
      function_name = "escoSimulate"
    ),
    param_integer(
      id = "lib.loc",
      default = 11L,
      lower = 1L,
      description = "Location (meanlog) parameter for the library size log-normal distribution, or mean parameter if a normal distribution is used.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "lib.scale",
      default = 0.2,
      lower = 0,
      border = FALSE,
      description = "Scale (sdlog) parameter for the library size log-normal distribution, or sd parameter if a normal distribution is used.",
      function_name = "escoSimulate"
    ),
    param_Boolean(
      id = "lib.norm",
      default = FALSE,
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "out.prob",
      default = 0.05,
      lower = 0,
      upper = 1,
      border = FALSE,
      description = "Probability that a gene is an expression outlier.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "out.facLoc",
      default = 4L,
      lower = 0,
      border = FALSE,
      description = "Location (meanlog) parameter for the expression outlier factor log-normal distribution.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "out.facScale",
      default = 0.5,
      lower = 0,
      border = FALSE,
      upper = 1,
      description = "Scale (sdlog) parameter for the expression outlier factor log-normal distribution.",
      function_name = "escoSimulate"
    ),
    param_integer(
      id = "nGroups",
      default = 1L,
      lower = 1L,
      description = "The number of groups or paths to simulate.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "group.prob",
      default = 1,
      lower = 0,
      border = TRUE,
      upper = 1,
      description = "Probability that a cell comes from a group.",
      function_name = "escoSimulate"
    ),
    param_others(
      id = "tree",
      description = "The tree structure to simulate.",
      type = "list",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "de.prob",
      default = 0.1,
      lower = 0,
      border = TRUE,
      upper = 1,
      description = "Probability that a gene is differentially expressed in a group. Can be a vector.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "de.facLoc",
      default = 0.1,
      lower = 0,
      border = FALSE,
      description = "Location (meanlog) parameter for the differential expression factor log-normal distribution. Can be a vector.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "de.facScale",
      default = 0.4,
      lower = 0,
      border = FALSE,
      description = "Scale (sdlog) parameter for the differential expression factor log-normal distribution. Can be a vector.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "de.center",
      default = 0,
      lower = 0,
      description = "The mean of the tree DE factors.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "de.loProb",
      default = 0,
      lower = 0,
      description = "The mean of the tree DE factors.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "de.downProb",
      default = 0.5,
      lower = 0,
      border = TRUE,
      upper = 1,
      description = "Probability that a differentially expressed gene is down-regulated. Can be a vector.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "bcv.common",
      default = 0.1,
      lower = 0,
      border = FALSE,
      description = "Underlying common dispersion across all genes.",
      function_name = "escoSimulate"
    ),
    param_integer(
      id = "bcv.df",
      default = 60L,
      lower = 1L,
      description = "Degrees of Freedom for the BCV inverse chi-squared distribution.",
      function_name = "escoSimulate"
    ),
    param_character(
      id = "dropout.type",
      default = "none",
      alternatives = c("none", "experiment", "batch", "group", "cell"),
      description = "The type of dropout to simulate. 'none' indicates no dropout, 'experiment' is global dropout using the same parameters for every cell, 'batch' uses the same parameters for every cell in each batch, 'group' uses the same parameters for every cell in each groups and 'cell' uses a different set of parameters for each cell.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "dropout.mid",
      default = 0,
      description = "Midpoint parameter for the dropout logistic function.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "dropout.shape",
      default = -1,
      description = "Shape parameter for the dropout logistic function.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "alpha_mean",
      default = 0.3,
      description = "Mean parameter for the dwonsampling gamma function.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "alpha_sd",
      default = 0.002,
      description = "Mean parameter for the dwonsampling gamma function.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "lenslope",
      default = 0.02,
      description = "Shape parameter for the dropout logistic function.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "nbins",
      default = 20,
      description = "Shape parameter for the dropout logistic function.",
      function_name = "escoSimulate"
    ),
    param_vector(
      id = "amp_bias_limt",
      default = c(-0.2, 0.2),
      description = "Shape parameter for the dropout logistic function.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "rate_2PCR",
      default = 0.8,
      description = "PCR efficiency, usually very high.",
      function_name = "escoSimulate"
    ),
    param_Boolean(
      id = "LinearAmp",
      default = FALSE,
      description = "If linear amplification is used for pre-amplification step, default is FALSE.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "LinearAmp_coef",
      default = 2000,
      description = "The coeficient of linear amplification, that is, how many times each molecule is amplified by.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "depth_mean",
      default = 1e+05,
      description = "Mean parameter of the sequencing depths.",
      function_name = "escoSimulate"
    ),
    param_numeric(
      id = "depth_sd",
      default = 1500,
      description = "Standard variance parameter of sequencing depths.",
      function_name = "escoSimulate"
    )
  )

  ESCO_method <- method_definition(
    method = "ESCO",
    programming = "R",
    url = "https://github.com/JINJINT/ESCO",
    authors = authors_definition(
      first = "Jinjin",
      last = "Tian",
      email = "roeder@andrew.cmu.edu",
      github = "https://github.com/JINJINT/ESCO",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "ESCO: single cell expression simulation incorporating gene co-expression",
      doi = "10.1093/bioinformatics/btab116",
      journal = "Bioinformatics",
      date = "2021",
      peer_review = TRUE
    ),
    description = "Single cell simulator with gene co-expression.",
    vignette = "http://47.254.148.113/software/Simsite/references/methods/21-esco/")

  list(ESCO_method = ESCO_method,
       ESCO_parameters = ESCO_parameters)
}
