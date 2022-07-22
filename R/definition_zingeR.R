#' Get Information of zingeR
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @importFrom splatter newSimpleParams
#' @export
#'
#' @examples
#' zingeR_method_definition <- zingeR_method_definition()

zingeR_method_definition <- function(...){
  zingeR_parameters <- parameter_sets(
    param_reference(
      id = "counts",
      type = c("matrix"),
      default = NULL,
      force = TRUE,
      process = "estimation",
      description = "A numeric matrix containing gene expression counts. Note that every gene in this matrix must have at least $p+1$ positive counts, with $p$ the number of columns in the design matrix.",
      function_name = "getDatasetZTNB"
    ),
    param_others(
      id = "design",
      type = "matrix",
      default = NULL,
      force = TRUE,
      process = "estimation",
      description = "The design of the experiments with rows corresponding to samples and columns corresponding to coefficients.",
      function_name = "getDatasetZTNB"
    ),
    param_others(
      id = "drop.extreme.dispersion",
      type = c("numeric", "Boolean"),
      default = FALSE,
      process = "estimation",
      description = "Either a numeric value between $0$ and $1$, stating the proportion of genes with extreme (high) dispersions to remove for simulation, or FALSE (default), if no dispersions should be removed for the analysis.",
      function_name = c("getDatasetZTNB", "NBsimSingleCell")
    ),
    param_others(
      id = "offset",
      type = "numeric",
      default = NULL,
      description = "The offset to use (typically the sequencing depth) when estimating gene-wise means and dispersions in the zero-truncated negative binomial model. These parameters will be used as a basis for the simulation.",
      function_name = "getDatasetZTNB"
    ),
    param_reference(
      id = "dataset",
      type = "matrix",
      force = TRUE,
      description = "An expression matrix representing the dataset on which the simulation is based.",
      function_name = "NBsimSingleCell"
    ),
    param_vector(
      id = "group",
      force = TRUE,
      description = "Group indicator specifying the attribution of the samples to the different conditions of interest that are being simulated.",
      function_name = "NBsimSingleCell"
    ),
    param_integer(
      id = "nTags",
      default = 10000L,
      lower = 1,
      description = "The number of features (genes) to simulate. $1000$ by default",
      function_name = "NBsimSingleCell"
    ),
    param_others(
      id = "nlibs",
      default = "length(group)",
      type = "integer",
      description = "The number of samples to simulate. Defaults to length(group).",
      function_name = "NBsimSingleCell"
    ),
    param_others(
      id = "lib.size",
      type = "numeric",
      description = "The library sizes for the simulated samples. If NULL (default), library sizes are resampled from the original datset.",
      function_name = "NBsimSingleCell"
    ),
    param_numeric(
      id = "pUp",
      default = 0.5,
      lower = 0,
      upper = 1,
      description = "Numeric value between $0$ and $1$ ($0.5$ by default) specifying the proportion of differentially expressed genes that show an upregulation in the second condition.",
      function_name = "NBsimSingleCell"
    ),
    param_numeric(
      id = "foldDiff",
      default = 2,
      description = "The fold changes used in simulating the differentially expressed genes. Either one numeric value for specifying the same fold change for all DE genes, or a vector of the same length as ind to specify fold changes for all differentially expressed genes. Note that fold changes above $1$ should be used as input of which a fraction will be inversed (i.e. simulation downregulation) according to 'pUp'. Defaults to $3$.",
      function_name = "NBsimSingleCell"
    ),
    param_Boolean(
      id = "verbose",
      default = TRUE,
      description = "Logical, stating whether progress be printed.",
      function_name = "NBsimSingleCell"
    ),
    param_vector(
      id = "ind",
      description = "Integer vector specifying the rows of the count matrix that represent differential features.",
      function_name = "NBsimSingleCell"
    ),
    param_others(
      id = "params",
      type = "created by getDatasetZTNB",
      description = "An object containing feature-wise parameters used for simulation as created by getDatasetZTNB. If NULL, parameters are estimated from the dataset provided.",
      function_name = "NBsimSingleCell"
    ),
    param_numeric(
      id = "randomZero",
      default = 0,
      upper = 1,
      lower = 0,
      description = "A numeric value between $0$ and $1$ specifying the random fraction of cells that are set to zero after simulating the expression count matrix. Defaults to $0$.",
      function_name = "NBsimSingleCell"
    ),
    param_numeric(
      id = "min.dispersion",
      default = 0.1,
      lower = 0,
      description = "The minimum dispersion value to use for simulation. $0.1$ by default.",
      function_name = "NBsimSingleCell"
    ),
    param_numeric(
      id = "max.dipserion",
      default = 400,
      lower = 0,
      description = "The maximum dispersion value to use for simulation. $400$ by default.",
      function_name = "NBsimSingleCell"
    )
  )

  zingeR_method <- method_definition(
    method = "zingeR",
    programming = "R",
    url = "https://github.com/statOmics/zingeR",
    authors = authors_definition(
      first = "Koen",
      last = "Van den Berge",
      email = "Koen.VanDenBerge@UGent.be",
      github = "https://github.com/statOmics/zingeR",
      orcid = "0000-0003-2214-0947"
    ),
    manuscript = manuscript_definition(
      title = "Observation weights unlock bulk RNA-seq tools for zero inflation and single-cell applications",
      doi = "10.1186/s13059-018-1406-4",
      journal = "Genome Biology",
      date = "2018",
      peer_review = TRUE
    ),
    description = "Zero-inflated negative binomial gene expression in R.")

  list(zingeR_method = zingeR_method,
       zingeR_parameters = zingeR_parameters)
}
