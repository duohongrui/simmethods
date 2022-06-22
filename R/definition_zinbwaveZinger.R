#' Get Information of zinbwaveZinger
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' zinbwaveZinger_method_definition <- zinbwaveZinger_method_definition()

zinbwaveZinger_method_definition <- function(...){
  zinbwaveZinger_parameters <- parameter_sets(
    param_reference(
      id = "counts",
      type = c("matrix"),
      default = NULL,
      force = TRUE,
      process = "estimation",
      description = "A numeric matrix containing gene expression counts. Note that every gene in this matrix must have at least $p+1$ positive counts, with $p$ the number of columns in the design matrix.",
      function_name = "getDatasetMoMPositive"
    ),
    param_others(
      id = "drop.extreme.dispersion",
      type = c("numeric", "Boolean"),
      default = FALSE,
      process = "estimation",
      description = "Either a numeric value between $0$ and $1$, stating the proportion of genes with extreme (high) dispersions to remove for simulation, or FALSE (default), if no dispersions should be removed for the analysis.",
      function_name = c("getDatasetMoMPositive", "NBsimSingleCell_zinbwaveZinger")
    ),
    param_character(
      id = "cpm",
      default = "AveLogCPM",
      alternatives = c("AveLogCPM"),
      description = "Normalization method",
      function_name = c("getDatasetMoMPositive", "NBsimSingleCell_zinbwaveZinger")
    ),
    param_numeric(
      id = "MoMIter",
      default = 10,
      function_name = "getDatasetMoMPositive"
    ),
    param_reference(
      id = "dataset",
      type = "matrix",
      force = TRUE,
      description = "An expression matrix representing the dataset on which the simulation is based.",
      function_name = "NBsimSingleCell_zinbwaveZinger"
    ),
    param_vector(
      id = "group",
      force = TRUE,
      description = "Group indicator specifying the attribution of the samples to the different conditions of interest that are being simulated.",
      function_name = "NBsimSingleCell_zinbwaveZinger"
    ),
    param_integer(
      id = "nTags",
      default = 10000L,
      lower = 1,
      description = "The number of features (genes) to simulate. $1000$ by default",
      function_name = "NBsimSingleCell_zinbwaveZinger"
    ),
    param_others(
      id = "nlibs",
      default = "length(group)",
      type = "integer",
      description = "The number of samples to simulate. Defaults to length(group).",
      function_name = "NBsimSingleCell_zinbwaveZinger"
    ),
    param_others(
      id = "lib.size",
      type = "numeric",
      description = "The library sizes for the simulated samples. If NULL (default), library sizes are resampled from the original datset.",
      function_name = "NBsimSingleCell_zinbwaveZinger"
    ),
    param_numeric(
      id = "pUp",
      default = 0.5,
      lower = 0,
      upper = 1,
      description = "Numeric value between $0$ and $1$ ($0.5$ by default) specifying the proportion of differentially expressed genes that show an upregulation in the second condition.",
      function_name = "NBsimSingleCell_zinbwaveZinger"
    ),
    param_numeric(
      id = "foldDiff",
      default = 3,
      description = "The fold changes used in simulating the differentially expressed genes. Either one numeric value for specifying the same fold change for all DE genes, or a vector of the same length as ind to specify fold changes for all differentially expressed genes. Note that fold changes above $1$ should be used as input of which a fraction will be inversed (i.e. simulation downregulation) according to 'pUp'. Defaults to $3$.",
      function_name = "NBsimSingleCell_zinbwaveZinger"
    ),
    param_Boolean(
      id = "verbose",
      default = TRUE,
      description = "Logical, stating whether progress be printed.",
      function_name = "NBsimSingleCell_zinbwaveZinger"
    ),
    param_vector(
      id = "ind",
      description = "Integer vector specifying the rows of the count matrix that represent differential features.",
      function_name = "NBsimSingleCell_zinbwaveZinger"
    ),
    param_others(
      id = "params",
      type = "created by getDatasetZTNB",
      description = "An object containing feature-wise parameters used for simulation as created by getDatasetZTNB. If NULL, parameters are estimated from the dataset provided.",
      function_name = "NBsimSingleCell_zinbwaveZinger"
    ),
    param_numeric(
      id = "min.dispersion",
      default = 0.1,
      lower = 0,
      description = "The minimum dispersion value to use for simulation. $0.1$ by default.",
      function_name = "NBsimSingleCell_zinbwaveZinger"
    ),
    param_numeric(
      id = "max.dipserion",
      default = 400,
      lower = 0,
      description = "The maximum dispersion value to use for simulation. $400$ by default.",
      function_name = "NBsimSingleCell_zinbwaveZinger"
    ),
    param_Boolean(
      id = "normalizeLambda",
      default = FALSE,
      function_name = "NBsimSingleCell_zinbwaveZinger"
    ),
    param_Boolean(
      id = "drop.low.lambda",
      default = TRUE,
      function_name = "NBsimSingleCell_zinbwaveZinger"
    )
  )

  zinbwaveZinger_method <- method_definition(
    method = "zinbwaveZinger",
    programming = "R",
    url = "https://github.com/statOmics/zinbwaveZinger",
    authors = authors_definition(
      first = "Koen",
      last = "Van den Berge",
      email = "Koen.VanDenBerge@UGent.be",
      github = "https://github.com/statOmics/zinbwaveZinger",
      orcid = "0000-0003-2214-0947"
    ),
    manuscript = manuscript_definition(
      title = "Observation weights unlock bulk RNA-seq tools for zero inflation and single-cell applications",
      doi = "10.1186/s13059-018-1406-4",
      journal = "Genome Biology",
      date = "2018",
      peer_review = TRUE
    ),
    description = "Integrating zingeR with ZINB-WaVE weights.")

  list(zinbwaveZinger_method = zinbwaveZinger_method,
       zinbwaveZinger_parameters = zinbwaveZinger_parameters)
}
