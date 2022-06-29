#' Get Information of SparseDC
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' SparseDC_method_definition <- SparseDC_method_definition()
#'
SparseDC_method_definition <- function(...){

  SparseDC_parameters <- parameter_sets(
    param_reference(
      id = "counts",
      type = "matrix",
      default = NULL,
      process = "estimation",
      force = TRUE,
      description = "Either a counts matrix or an SingleCellExperiment object containing count data to estimate parameters from.",
      function_name = "sparseDCEstimate"
    ),
    param_vector(
      id = "conditions",
      default = FALSE,
      force = TRUE,
      process = "estimation",
      description = "Numeric vector giving the condition each cell belongs to.",
      function_name = "sparseDCEstimate"
    ),
    param_others(
      id = "nclusters",
      type = "integer",
      force = TRUE,
      process = "estimation",
      description = "Number of cluster present in the dataset.",
      function_name = "sparseDCEstimate"
    ),
    param_Boolean(
      id = "norm",
      default = TRUE,
      process = "estimation",
      description = "Logical, whether to library size normalise counts before estimation. Set this to FALSE if counts is already normalised.",
      function_name = "sparseDCEstimate"
    ),
    param_others(
      id = "params",
      type = "SparseDCParams",
      default = "newSparseDCParams()",
      process = "estimation",
      description = "PhenoParams object to store estimated values in.",
      function_name = "sparseDCEstimate"
    ),
    param_Boolean(
      id = "sparsify",
      default = TRUE,
      process = "estimation",
      description = "Logical. Whether to automatically convert assays to sparse matrices if there will be a size reduction.",
      function_name = "sparseDCSimulate"
    ),
    param_Boolean(
      id = "verbose",
      default = TRUE,
      process = "estimation",
      description = "Logical. Whether to print progress messages.",
      function_name = "sparseDCSimulate"
    ),
    param_integer(
      id = "nGenes",
      process = "simulation",
      default = 10000L,
      lower = 1L,
      function_name = "sparseDCSimulate"
    ),
    param_integer(
      id = "nCells",
      process = "simulation",
      default = 100L,
      lower = 1L,
      function_name = "sparseDCSimulate"
    ),
    param_integer(
      id = "seed",
      process = "simulation",
      default = 431492L,
      function_name = "sparseDCSimulate"
    ),
    param_integer(
      id = "markers.n",
      default = 0L,
      process = "simulation",
      description = "Number of marker genes to simulate for each cluster.",
      function_name = "sparseDCSimulate"
    ),
    param_integer(
      id = "markers.shared",
      default = 0L,
      process = "simulation",
      description = "Number of marker genes for each cluster shared between conditions. Must be less than or equal to markers.n",
      function_name = "sparseDCSimulate"
    ),
    param_Boolean(
      id = "markers.same",
      default = FALSE,
      process = "simulation",
      description = "Logical. Whether to print progress messages.",
      function_name = "sparseDCSimulate"
    ),
    param_vector(
      id = "clusts.c1",
      default = 1,
      process = "simulation",
      description = "Numeric vector of clusters present in condition 1. The number of times a cluster is repeated controls the proportion of cells from that cluster.",
      function_name = "sparseDCSimulate"
    ),
    param_vector(
      id = "clusts.c2",
      default = 1,
      process = "simulation",
      description = "Numeric vector of clusters present in condition 2. The number of times a cluster is repeated controls the proportion of cells from that cluster.",
      function_name = "sparseDCSimulate"
    ),
    param_numeric(
      id = "mean.lower",
      default = 1,
      lower = 0,
      description = "Lower bound for cluster gene means.",
      process = "simulation",
      function_name = "sparseDCSimulate"
    ),
    param_numeric(
      id = "mean.upper",
      default = 2,
      lower = 0,
      description = "Upper bound for cluster gene means.",
      process = "simulation",
      function_name = "sparseDCSimulate"
    )
  )

  SparseDC_method <- method_definition(
    method = "SparseDC",
    programming = "R",
    url = "https://cran.rstudio.com/web/packages/SparseDC/index.html",
    authors = authors_definition(
      first = "Martin",
      last = "Barron",
      email = "jun.li@nd.edu",
      github = "benjamin.schuster-boeckler@ludwig.ox.ac.uk",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "A sparse differential clustering algorithm for tracing cell type changes via single-cell RNA-sequencing data",
      doi = "10.1093/nar/gkx1113",
      journal = "Nucleic Acids Research",
      date = "2018",
      peer_review = TRUE
    ),
    description = "This algorithm clusters samples from two different populations, links the clusters across the conditions and identifies marker genes for these changes. The package was designed for scRNA-Seq data but is also applicable to many other data types, just replace cells with samples and genes with variables.")

  list(SparseDC_method = SparseDC_method,
       SparseDC_parameters = SparseDC_parameters)
}
