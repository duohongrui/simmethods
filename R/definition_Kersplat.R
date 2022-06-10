#' Get Information of Kersplat
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @importFrom splatter newKersplatParams
#' @export
#'
#' @examples
#' Kersplat_method_definition <- Kersplat_method_definition()

Kersplat_method_definition <- function(...){
  Kersplat_parameters <- parameter_sets(
    param_reference(
      id = "input_refernece",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      description = "The reference data to go through estimation and evaluation process. Usually, no reference data is also OK."
    ),
    param_others(
      id = "KersplatParams",
      type = "KersplatParams",
      default = "splatter::newKersplatParams()",
      process = "estimation",
      description = "Usually it is default by splatter::newKersplatParams function. Users can change the parameters by splatter::setParam function."
    ),
    param_numeric(
      id = "mean.shape",
      default = 0.6,
      lower = 0,
      border = FALSE,
      description = "Shape parameter for the mean gamma distribution."
    ),
    param_numeric(
      id = "mean.rate",
      default = 0.3,
      lower = 0,
      border = FALSE,
      description = "Rate parameter for the mean gamma distribution."
    ),
    param_numeric(
      id = "mean.outProb",
      default = 0.05,
      lower = 0,
      upper = 1,
      border = FALSE,
      description = "Probability that a gene is an expression outlier."
    ),
    param_numeric(
      id = "mean.outLoc",
      default = 4L,
      lower = 0,
      border = FALSE,
      description = "Location (meanlog) parameter for the expression outlier factor log-normal distribution."
    ),
    param_numeric(
      id = "mean.outScale",
      default = 0.5,
      lower = 0,
      border = FALSE,
      upper = 1,
      description = "Scale (sdlog) parameter for the expression outlier factor log-normal distribution."
    ),
    param_character(
      id = "mean.method",
      alternatives = c("fit", "density"),
      default = "fit",
      description = "Method to use for simulating gene means. Either 'fit' to sample from a gamma distribution (with expression outliers) or 'density' to sample from the provided density object."
    ),
    param_vector(
      id = "mean.values",
      default = NULL,
      description = "Vector of means for each gene."
    ),
    param_numeric(
      id = "bcv.common",
      default = 0.1,
      lower = 0,
      border = FALSE,
      description = "Underlying common dispersion across all genes."
    ),
    param_integer(
      id = "bcv.df",
      default = 60L,
      lower = 1L,
      description = "Degrees of Freedom for the BCV inverse chi-squared distribution."
    ),
    param_others(
      id = "network.graph",
      default = NULL,
      type = "list",
      description = "Graph containing the gene network."
    ),
    param_numeric(
      id = "network.nRegs",
      default = 100,
      description = "Number of regulators in the network."
    ),
    param_Boolean(
      id = "network.regsSet",
      default = FALSE
    ),
    param_numeric(
      id = "paths.nPrograms",
      default = 10,
      description = "Number of expression programs."
    ),
    param_dataframe(
      id = "paths.design",
      type = "data.frame",
      description = "data.frame describing path structure. See kersplatSimPaths for details."
    ),
    param_others(
      id = "paths.means",
      type = "list",
      default = "list()"
    ),
    param_integer(
      id = "lib.loc",
      default = 11L,
      lower = 1L,
      description = "Location (meanlog) parameter for the library size log-normal distribution, or mean parameter if a normal distribution is used."
    ),
    param_numeric(
      id = "lib.scale",
      default = 0.2,
      lower = 0,
      border = FALSE,
      description = "Scale (sdlog) parameter for the library size log-normal distribution, or sd parameter if a normal distribution is used."
    ),
    param_character(
      id = "lib.method",
      default = "fit",
      description = "Method to use for simulating library sizes. Either 'fit' to sample from a log-normal distribution or 'density' to sample from the provided density object.",
      alternatives = c("fit", "density")
    ),
    param_dataframe(
      id = "cells.design",
      type = "data.frame",
      description = "data.frame describing cell structure. See kersplatSimCellMeans for details."
    ),
    param_numeric(
      id = "doublet.prop",
      default = 0,
      description = "Proportion of cells that are doublets."
    ),
    param_numeric(
      id = "ambient.scale",
      default = 0.05,
      description = "Scaling factor for the library size log-normal distribution when generating ambient library sizes."
    ),
    param_numeric(
      id = "ambient.nEmpty",
      default = 0,
      description = "Number of empty cells to simulate."
    ),
    param_integer(
      id = "nGenes",
      default = 10000L,
      lower = 1L
    ),
    param_integer(
      id = "nCells",
      default = 100L,
      lower = 1L
    ),
    param_integer(
      id = "seed",
      default = 633483L
    )
  )

  Kersplat_method <- method_definition(
    method = "Kersplat",
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

  list(Kersplat_method = Kersplat_method,
       Kersplat_parameters = Kersplat_parameters)
}
