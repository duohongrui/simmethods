#' Get Information of SCRIP
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @importFrom splatter newSplatParams
#' @export
#'
#' @examples
#' Splat_method_definition <- Splat_method_definition()
#'
SCRIP_method_definition <- function(...){

  SCRIP_parameters <- parameter_sets(
    param_reference(
      id = "input_refernece",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      description = "The reference data to go through estimation and evaluation process. Usually, no reference data is also OK.",
      function_name = "splatEstimate"
    ),
    param_others(
      id = "params",
      type = "SplatParams",
      default = "splatter::newSplatParams()",
      process = "simulation",
      description = "Usually it is default by splatter::newSplatParams function. Users can change the parameters by splatter::setParam function.",
      function_name = "splatEstimate"
    ),
    param_character(
      id = "method",
      default = "single",
      alternatives = c("single", "groups", "paths"),
      description = "Which simulation method to use. Options are 'single' which produces a single population, 'groups' which produces distinct groups (eg. cell types), 'paths' which selects cells from continuous trajectories (eg. differentiation processes).",
      process = "simulation",
      function_name = "SCRIPsimu"
    ),
    param_integer(
      id = "nBatches",
      default = 1L,
      lower = 1L,
      description = "The number of batches to simulate.",
      function_name = "SCRIPsimu"
    ),
    param_character(
      id = "mode",
      alternatives = c("GP-commonBCV", "BP-commonBCV", "BP", "BGP-commonBCV", "BGP-trendedBCV"),
      force = TRUE,
      description = "mode",
      default = "GP-commonBCV",
      process = "simulation",
      function_name = "SCRIPsimu"
    ),
    param_vector(
      id = "base_allcellmeans_SC",
      default = NULL,
      description = "base mean vector provided to help setting DE analysis",
      process = "simulation",
      function_name = "SCRIPsimu"
    ),
    param_dataframe(
      id = "pre.bcv.df",
      description = "BCV.df enables us to change the variation of BCV values",
      process = "simulation",
      function_name = "SCRIPsimu"
    ),
    param_vector(
      id = "libsize",
      default = NULL,
      description = "library size can be provided directly",
      process = "simulation",
      function_name = "SCRIPsimu"
    ),
    param_integer(
      id = "batchCells",
      default = 100L,
      lower = 1L,
      description = "Vector giving the number of cells in each batch.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "bcv.shrink",
      default = 1,
      description = "factor to control the BCV levels",
      function_name = "SCRIPsimu"
    ),
    param_others(
      id = "Dropout_rate",
      default = NULL,
      type = "nnumeric",
      description = "factor to control the dropout rate directly",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "batch.facLoc",
      default = 0.1,
      lower = 0,
      border = FALSE,
      upper = 1,
      description = "Location (meanlog) parameter for the batch effect factor log-normal distribution. Can be a vector.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "batch.facScale",
      default = 0.1,
      lower = 0,
      border = FALSE,
      upper = 1,
      description = "Scale (sdlog) parameter for the batch effect factor log-normal distribution. Can be a vector.",
      function_name = "SCRIPsimu"
    ),
    param_Boolean(
      id = "batch.rmEffect",
      default = FALSE,
      description = "Logical, removes the batch effect and continues with the simulation when TRUE. This allows the user to test batch removal algorithms without having to calculate the new expected cell means with batch removed.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "mean.shape",
      default = 0.6,
      lower = 0,
      border = FALSE,
      description = "Shape parameter for the mean gamma distribution.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "mean.rate",
      default = 0.3,
      lower = 0,
      border = FALSE,
      description = "Rate parameter for the mean gamma distribution.",
      function_name = "SCRIPsimu"
    ),
    param_integer(
      id = "lib.loc",
      default = 11L,
      lower = 1L,
      description = "Location (meanlog) parameter for the library size log-normal distribution, or mean parameter if a normal distribution is used.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "lib.scale",
      default = 0.2,
      lower = 0,
      border = FALSE,
      description = "Scale (sdlog) parameter for the library size log-normal distribution, or sd parameter if a normal distribution is used.",
      function_name = "SCRIPsimu"
    ),
    param_Boolean(
      id = "lib.norm",
      default = FALSE,
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "out.prob",
      default = 0.05,
      lower = 0,
      upper = 1,
      border = FALSE,
      description = "Probability that a gene is an expression outlier.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "out.facLoc",
      default = 4L,
      lower = 0,
      border = FALSE,
      description = "Location (meanlog) parameter for the expression outlier factor log-normal distribution.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "out.facScale",
      default = 0.5,
      lower = 0,
      border = FALSE,
      upper = 1,
      description = "Scale (sdlog) parameter for the expression outlier factor log-normal distribution.",
      function_name = "SCRIPsimu"
    ),
    param_integer(
      id = "nGroups",
      default = 1L,
      lower = 1L,
      description = "The number of groups or paths to simulate.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "group.prob",
      default = 1,
      lower = 0,
      border = TRUE,
      upper = 1,
      description = "Probability that a cell comes from a group.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "de.prob",
      default = 0.1,
      lower = 0,
      border = TRUE,
      upper = 1,
      description = "Probability that a gene is differentially expressed in a group. Can be a vector.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "de.downProb",
      default = 0.5,
      lower = 0,
      border = TRUE,
      upper = 1,
      description = "Probability that a differentially expressed gene is down-regulated. Can be a vector.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "de.facLoc",
      default = 0.1,
      lower = 0,
      border = FALSE,
      description = "Location (meanlog) parameter for the differential expression factor log-normal distribution. Can be a vector.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "de.facScale",
      default = 0.4,
      lower = 0,
      border = FALSE,
      description = "Scale (sdlog) parameter for the differential expression factor log-normal distribution. Can be a vector.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "bcv.common",
      default = 0.1,
      lower = 0,
      border = FALSE,
      description = "Underlying common dispersion across all genes.",
      function_name = "SCRIPsimu"
    ),
    param_integer(
      id = "bcv.df",
      default = 60L,
      lower = 1L,
      description = "Degrees of Freedom for the BCV inverse chi-squared distribution.",
      function_name = "SCRIPsimu"
    ),
    param_character(
      id = "dropout.type",
      default = "none",
      alternatives = c("none", "experiment", "batch", "group", "cell"),
      description = "The type of dropout to simulate. 'none' indicates no dropout, 'experiment' is global dropout using the same parameters for every cell, 'batch' uses the same parameters for every cell in each batch, 'group' uses the same parameters for every cell in each groups and 'cell' uses a different set of parameters for each cell.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "dropout.mid",
      default = 0,
      description = "Midpoint parameter for the dropout logistic function.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "dropout.shape",
      default = -1,
      description = "Shape parameter for the dropout logistic function.",
      function_name = "SCRIPsimu"
    ),
    param_vector(
      id = "path.from",
      default = 0,
      description = "Vector giving the originating point of each path. This allows path structure such as a cell type which differentiates into an intermediate cell type that then differentiates into two mature cell types. A path structure of this form would have a 'from' parameter of c(0, 1, 1) (where 0 is the origin). If no vector is given all paths will start at the origin.",
      function_name = "SCRIPsimu"
    ),
    param_integer(
      id = "path.nSteps",
      default = 100L,
      lower = 1L,
      description = "Vector giving the number of steps to simulate along each path. If a single value is given it will be applied to all paths. This parameter was previously called path.length.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "path.skew",
      default = 0.5,
      lower = 0,
      upper = 1,
      description = "Vector giving the skew of each path. Values closer to 1 will give more cells towards the starting population, values closer to 0 will give more cells towards the final population. If a single value is given it will be applied to all paths.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "path.nonlinearProb",
      default = 0.1,
      lower = 0,
      upper = 1,
      description = "Probability that a gene follows a non-linear path along the differentiation path. This allows more complex gene patterns such as a gene being equally expressed at the beginning an end of a path but lowly expressed in the middle.",
      function_name = "SCRIPsimu"
    ),
    param_numeric(
      id = "path.sigmaFac",
      default = 0.8,
      lower = 0,
      upper = 1,
      description = "Sigma factor for non-linear gene paths. A higher value will result in more extreme non-linear variations along a path.",
      function_name = "SCRIPsimu"
    ),
    param_integer(
      id = "nGenes",
      default = 10000L,
      lower = 1L,
      description = "The number of genes to simulate.",
      function_name = "SCRIPsimu"
    ),
    param_integer(
      id = "nCells",
      default = 100L,
      lower = 1L,
      description = "The number of cells to simulate.",
      function_name = "SCRIPsimu"
    ),
    param_integer(
      id = "seed",
      default = 687680L,
      description = "Seed to use for generating random numbers.",
      function_name = "SCRIPsimu"
    )
  )

  SCRIP_method <- method_definition(
    method = "SCRIP",
    programming = "R",
    url = "https://cran.r-project.org/web/packages/SCRIP/index.html",
    authors = authors_definition(
      first = "Fei",
      last = "Qin",
      email = "fqin@email.sc.edu",
      github = "https://github.com/thecailab/SCRIP",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "SCRIP: an accurate simulator for single-cell RNA sequencing data",
      doi = "10.1093/bioinformatics/btab824",
      journal = "Bioinformatics",
      date = "2022",
      peer_review = TRUE
    ),
    description = "We provide a comprehensive scheme that is capable of simulating Single Cell RNA Sequencing data for various parameters of Biological Coefficient of Variation, busting kinetics, differential expression (DE), cell or sample groups, cell trajectory, batch effect and other experimental designs. 'SCRIP' proposed and compared two frameworks with Gamma-Poisson and Beta-Gamma-Poisson models for simulating Single Cell RNA Sequencing data.")

  list(SCRIP_method = SCRIP_method,
       SCRIP_parameters = SCRIP_parameters)
}
