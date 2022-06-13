#' Get Information of Splat
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
Splat_method_definition <- function(...){

  Splat_parameters <- parameter_sets(
    param_reference(
      id = "input_refernece",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      description = "The reference data to go through estimation and evaluation process. Usually, no reference data is also OK.",
      function_name = "splatEstimate"
    ),
    param_others(
      id = "SplatParams",
      type = "SplatParams",
      default = "splatter::newSplatParams()",
      process = "estimation",
      description = "Usually it is default by splatter::newSplatParams function. Users can change the parameters by splatter::setParam function.",
      function_name = "splatEstimate"
    ),
    param_character(
      id = "method",
      default = "single",
      alternatives = c("single", "groups", "paths"),
      description = "Which simulation method to use. Options are 'single' which produces a single population, 'groups' which produces distinct groups (eg. cell types), 'paths' which selects cells from continuous trajectories (eg. differentiation processes).",
      function_name = "splatSimulate"
    ),
    param_integer(
      id = "nBatches",
      default = 1L,
      lower = 1L,
      description = "The number of batches to simulate.",
      function_name = "splatSimulate"
    ),
    param_integer(
      id = "batchCells",
      default = 100L,
      lower = 1L,
      description = "Vector giving the number of cells in each batch.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "batch.facLoc",
      default = 0.1,
      lower = 0,
      border = FALSE,
      upper = 1,
      description = "Location (meanlog) parameter for the batch effect factor log-normal distribution. Can be a vector.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "batch.facScale",
      default = 0.1,
      lower = 0,
      border = FALSE,
      upper = 1,
      description = "Scale (sdlog) parameter for the batch effect factor log-normal distribution. Can be a vector.",
      function_name = "splatSimulate"
    ),
    param_Boolean(
      id = "batch.rmEffect",
      default = FALSE,
      description = "Logical, removes the batch effect and continues with the simulation when TRUE. This allows the user to test batch removal algorithms without having to calculate the new expected cell means with batch removed.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "mean.shape",
      default = 0.6,
      lower = 0,
      border = FALSE,
      description = "Shape parameter for the mean gamma distribution.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "mean.rate",
      default = 0.3,
      lower = 0,
      border = FALSE,
      description = "Rate parameter for the mean gamma distribution.",
      function_name = "splatSimulate"
    ),
    param_integer(
      id = "lib.loc",
      default = 11L,
      lower = 1L,
      description = "Location (meanlog) parameter for the library size log-normal distribution, or mean parameter if a normal distribution is used.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "lib.scale",
      default = 0.2,
      lower = 0,
      border = FALSE,
      description = "Scale (sdlog) parameter for the library size log-normal distribution, or sd parameter if a normal distribution is used.",
      function_name = "splatSimulate"
    ),
    param_Boolean(
      id = "lib.norm",
      default = FALSE,
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "out.prob",
      default = 0.05,
      lower = 0,
      upper = 1,
      border = FALSE,
      description = "Probability that a gene is an expression outlier.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "out.facLoc",
      default = 4L,
      lower = 0,
      border = FALSE,
      description = "Location (meanlog) parameter for the expression outlier factor log-normal distribution.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "out.facScale",
      default = 0.5,
      lower = 0,
      border = FALSE,
      upper = 1,
      description = "Scale (sdlog) parameter for the expression outlier factor log-normal distribution.",
      function_name = "splatSimulate"
    ),
    param_integer(
      id = "nGroups",
      default = 1L,
      lower = 1L,
      description = "The number of groups or paths to simulate.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "group.prob",
      default = 1,
      lower = 0,
      border = TRUE,
      upper = 1,
      description = "Probability that a cell comes from a group.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "de.prob",
      default = 0.1,
      lower = 0,
      border = TRUE,
      upper = 1,
      description = "Probability that a gene is differentially expressed in a group. Can be a vector.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "de.downProb",
      default = 0.5,
      lower = 0,
      border = TRUE,
      upper = 1,
      description = "Probability that a differentially expressed gene is down-regulated. Can be a vector.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "de.facLoc",
      default = 0.1,
      lower = 0,
      border = FALSE,
      description = "Location (meanlog) parameter for the differential expression factor log-normal distribution. Can be a vector.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "de.facScale",
      default = 0.4,
      lower = 0,
      border = FALSE,
      description = "Scale (sdlog) parameter for the differential expression factor log-normal distribution. Can be a vector.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "bcv.common",
      default = 0.1,
      lower = 0,
      border = FALSE,
      description = "Underlying common dispersion across all genes.",
      function_name = "splatSimulate"
    ),
    param_integer(
      id = "bcv.df",
      default = 60L,
      lower = 1L,
      description = "Degrees of Freedom for the BCV inverse chi-squared distribution.",
      function_name = "splatSimulate"
    ),
    param_character(
      id = "dropout.type",
      default = "none",
      alternatives = c("none", "experiment", "batch", "group", "cell"),
      description = "The type of dropout to simulate. 'none' indicates no dropout, 'experiment' is global dropout using the same parameters for every cell, 'batch' uses the same parameters for every cell in each batch, 'group' uses the same parameters for every cell in each groups and 'cell' uses a different set of parameters for each cell.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "dropout.mid",
      default = 0,
      description = "Midpoint parameter for the dropout logistic function.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "dropout.shape",
      default = -1,
      description = "Shape parameter for the dropout logistic function.",
      function_name = "splatSimulate"
    ),
    param_vector(
      id = "path.from",
      default = 0,
      description = "Vector giving the originating point of each path. This allows path structure such as a cell type which differentiates into an intermediate cell type that then differentiates into two mature cell types. A path structure of this form would have a 'from' parameter of c(0, 1, 1) (where 0 is the origin). If no vector is given all paths will start at the origin.",
      function_name = "splatSimulate"
    ),
    param_integer(
      id = "path.nSteps",
      default = 100L,
      lower = 1L,
      description = "Vector giving the number of steps to simulate along each path. If a single value is given it will be applied to all paths. This parameter was previously called path.length.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "path.skew",
      default = 0.5,
      lower = 0,
      upper = 1,
      description = "Vector giving the skew of each path. Values closer to 1 will give more cells towards the starting population, values closer to 0 will give more cells towards the final population. If a single value is given it will be applied to all paths.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "path.nonlinearProb",
      default = 0.1,
      lower = 0,
      upper = 1,
      description = "Probability that a gene follows a non-linear path along the differentiation path. This allows more complex gene patterns such as a gene being equally expressed at the beginning an end of a path but lowly expressed in the middle.",
      function_name = "splatSimulate"
    ),
    param_numeric(
      id = "path.sigmaFac",
      default = 0.8,
      lower = 0,
      upper = 1,
      description = "Sigma factor for non-linear gene paths. A higher value will result in more extreme non-linear variations along a path.",
      function_name = "splatSimulate"
    ),
    param_integer(
      id = "nGenes",
      default = 10000L,
      lower = 1L,
      description = "The number of genes to simulate.",
      function_name = "splatSimulate"
    ),
    param_integer(
      id = "nCells",
      default = 100L,
      lower = 1L,
      description = "The number of cells to simulate.",
      function_name = "splatSimulate"
    ),
    param_integer(
      id = "seed",
      default = 687680L,
      description = "Seed to use for generating random numbers.",
      function_name = "splatSimulate"
    )
  )

  Splat_method <- method_definition(
    method = "Splat",
    programming = "R",
    url = "https://bioconductor.org/packages/release/bioc/html/splatter.html",
    authors = authors_definition(
      first = "Luke",
      last = "Zappia",
      email = NULL,
      github = "https://github.com/Oshlack/splatter",
      orcid = "0000-0001-7744-8565"
    ),
    manuscript = manuscript_definition(
      title = "Splatter: simulation of single-cell RNA sequencing data",
      doi = "10.1186/s13059-017-1305-0",
      journal = "Genome Biology",
      date = "2017",
      peer_review = TRUE
    ),
    description = "Splatter is a package for the simulation of single-cell RNA sequencing count data")

  list(Splat_method = Splat_method,
       Splat_parameters = Splat_parameters)
}
