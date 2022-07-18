#' Get Information of Lun
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @importFrom splatter newLunParams
#' @export
#'
#' @examples
#' Lun_method_definition <- Lun_method_definition()

Lun_method_definition <- function(...){
  Lun_parameters <- parameter_sets(
    param_reference(
      id = "input_refernece",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      description = "The reference data to go through estimation and evaluation process. Usually, no reference data is also OK.",
      function_name = "lunEstimate"
    ),
    param_others(
      id = "LunParams",
      type = "LunParams",
      default = "splatter::newLunParams()",
      process = "estimation",
      description = "Usually it is default by splatter::newLunParams function. Users can change the parameters by splatter::setParam function.",
      function_name = "lunEstimate"
    ),
    param_integer(
      id = "nGroups",
      default = 1L,
      lower = 1L,
      description = "The number of groups to simulate.",
      function_name = "lunSimulate"
    ),
    param_integer(
      id = "groupCells",
      default = 100L,
      lower = 1L,
      description = "Vector giving the number of cells in each simulation group/path.",
      function_name = "lunSimulate"
    ),
    param_integer(
      id = "de.nGenes",
      default = 1000L,
      lower = 1L,
      description = "The number of genes that are differentially expressed in each group",
      function_name = "lunSimulate"
    ),
    param_numeric(
      id = "de.upProp",
      default = 0.5,
      lower = 0,
      upper = 1,
      description = "The proportion of differentially expressed genes that are up-regulated in each group.",
      function_name = "lunSimulate"
    ),
    param_numeric(
      id = "de.upFC",
      default = 5,
      lower = 0,
      description = "The fold change for up-regulated genes.",
      function_name = "lunSimulate"
    ),
    param_numeric(
      id = "de.downFC",
      default = 0,
      upper = 0,
      description = "The fold change for down-regulated genes.",
      function_name = "lunSimulate"
    ),
    param_numeric(
      id = "mean.shape",
      default = 2,
      lower = 0,
      border = FALSE,
      description = "Shape parameter for the mean gamma distribution.",
      function_name = "lunSimulate"
    ),
    param_numeric(
      id = "mean.rate",
      default = 2,
      lower = 0,
      border = FALSE,
      description = "Rate parameter for the mean gamma distribution.",
      function_name = "lunSimulate"
    ),
    param_numeric(
      id = "count.disp",
      default = 0.1,
      lower = 0,
      border = FALSE,
      description = "The dispersion parameter for the counts negative binomial distribution.",
      function_name = "lunSimulate"
    ),
    param_integer(
      id = "nGenes",
      default = 10000L,
      lower = 1L,
      description = "The number of genes to simulate.",
      function_name = "lunSimulate"
    ),
    param_integer(
      id = "nCells",
      default = 100L,
      lower = 1L,
      description = "The number of cells to simulate.",
      function_name = "lunSimulate"
    ),
    param_integer(
      id = "seed",
      force = TRUE,
      description = "Seed to use for generating random numbers.",
      function_name = "lunSimulate"
    )
  )

  Lun_method <- method_definition(
    method = "Lun",
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

  list(Lun_method = Lun_method,
       Lun_parameters = Lun_parameters)
}
