#' Get Information of Lun2
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @importFrom splatter newLun2Params
#' @export
#'
#' @examples
#' Lun2_method_definition <- Lun2_method_definition()

Lun2_method_definition <- function(...){
  Lun2_parameters <- parameter_sets(
    param_reference(
      id = "input_refernece",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      description = "The reference data to go through estimation and evaluation process. Usually, no reference data is also OK."
    ),
    param_others(
      id = "Lun2Params",
      type = "Lun2Params",
      default = "splatter::newLun2Params()",
      process = "estimation",
      description = "Usually it is default by splatter::newLun2Params function. Users can change the parameters by splatter::setParam function."
    ),
    param_vector(
      id = "plates",
      description = "Integer vector giving the plate that each cell originated from.",
      process = "estimation"
    ),
    param_integer(
      id = "min.size",
      default = 200L,
      lower = 1L,
      description = "Minimum size of clusters when identifying group of cells in the data.",
      process = "estimation"
    ),
    param_others(
      id = "BPPARAM",
      default = "SerialParam()",
      type = "SerialParam()",
      process = "estimation",
      description = "A BiocParallelParam instance giving the parallel back-end to be used. Default is SerialParam which uses a single core."
    ),
    param_Boolean(
      id = "zinb",
      default = FALSE,
      description = "logical. Whether to use a zero-inflated model."
    ),
    param_integer(
      id = "nPlates",
      default = 1L,
      lower = 1L,
      description = "The number of plates to simulate."
    ),
    param_vector(
      id = "plate.ingroup",
      default = "1",
      description = "Character vector giving the plates considered to be part of the 'ingroup'."
    ),
    param_numeric(
      id = "plate.mod",
      default = 1,
      description = "Plate effect modifier factor. The plate effect variance is divided by this value."
    ),
    param_numeric(
      id = "plate.var",
      default = 14,
      description = "Plate effect variance."
    ),
    param_dataframe(
      id = "gene.params",
      description = "A data.frame containing gene parameters with two columns: Mean (mean expression for each gene) and Disp (dispersion for each gene)."
    ),
    param_dataframe(
      id = "zi.params",
      description = "A data.frame containing zero-inflated gene parameters with three columns: Mean (mean expression for each gene), Disp (dispersion for each, gene), and Prop (zero proportion for each gene)."
    ),
    param_others(
      id = "cell.plates",
      type = "factor",
      default = "factor(rep(1,100))",
      description = "Factor giving the plate that each cell comes from."
    ),
    param_vector(
      id = "cell.libSizes",
      default = "rep(70000, 100)",
      description = "Library size for each cell."
    ),
    param_numeric(
      id = "cell.libMod",
      default = 1,
      upper = 1,
      lower = 0,
      description = "Modifier factor for library sizes. The library sizes are multiplied by this value."
    ),
    param_integer(
      id = "de.nGenes",
      default = 0L,
      lower = 0L,
      description = "Number of differentially expressed genes."
    ),
    param_numeric(
      id = "de.fc",
      default = 3,
      description = "Fold change for differentially expressed genes."
    ),
    param_integer(
      id = "nGenes",
      default = 10000L,
      lower = 1L,
      description = "The number of genes to simulate."
    ),
    param_integer(
      id = "nCells",
      default = 100L,
      lower = 1L,
      description = "The number of cells to simulate."
    ),
    param_integer(
      id = "seed",
      default = 41234L,
      description = "Seed to use for generating random numbers."
    )
  )

  Lun2_method <- method_definition(
    method = "Lun2",
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

  list(Lun2_method = Lun2_method,
       Lun2_parameters = Lun2_parameters)
}
