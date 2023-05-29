#' Get Information of scDesign3
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' muscat_method_definition <- muscat_method_definition()
#'
scDesign3_method_definition <- function(...){

  scDesign3_parameters <- parameter_sets(
    param_reference(
      id = "sce",
      type = "SingleCellExperiment",
      default = NULL,
      process = "estimation",
      force = TRUE,
      description = "A SingleCellExperiment.",
      function_name = c("construct_data", "simData")
    ),
    param_vector(
      id = "spatial",
      default = NULL,
      description = "A length two string vector of the names of spatial coordinates. Defualt is NULL.",
      process = "estimation",
      function_name = "construct_data"
    ),
    param_vector(
      id = "other_covariates",
      default = NULL,
      description = "A string or a string vector of the other covaraites you want to include in the data.",
      process = "estimation",
      function_name = "construct_data"
    ),
    param_vector(
      id = "ncell",
      default = "dim(sce)[2]",
      process = "estimation",
      description = "The number of cell you want to simulate. Default is dim(sce)[2] (the same number as the input data). If an arbitrary number is provided, the fucntion will use Vine Copula to simulate a new covaraite matrix.",
      function_name = "construct_data"
    ),
    param_vector(
      id = "predictor",
      default = "gene",
      description = "A string of the predictor for the gam/gamlss model. Default is gene. This is essentially just a name.",
      process = "estimation",
      function_name = "fit_marginal"
    ),
    param_vector(
      id = "mu_formula",
      default = "1",
      description = "A string of the mu parameter formula",
      process = "estimation",
      function_name = "fit_marginal"
    ),
    param_vector(
      id = "sigma_formula",
      default = "1",
      description = "A string of the sigma parameter formula",
      process = "estimation",
      function_name = "fit_marginal"
    ),
    param_character(
      id = "family_use",
      default = "nb",
      alternatives = c("nb", "binomial", "poisson", "zip", "zinb", "guassian"),
      description = "A string or a vector of strings of the marginal distribution. Must be one of 'binomial', 'poisson', 'nb', 'zip', 'zinb' or 'gaussian', which represent 'poisson distribution', 'negative binomial distribution', 'zero-inflated poisson distribution', 'zero-inflated negative binomail distribution', and 'gaussian distribution' respectively.",
      process = "estimation",
      function_name = "fit_marginal"
    ),
    param_integer(
      id = "n_cores",
      default = 1L,
      lower = 1L,
      description = "An integer. The number of cores to use.",
      process = "estimation",
      function_name = "fit_marginal"
    ),
    param_vector(
      id = "assay_use",
      default = "counts",
      description = "A string which indicates the assay you will use in the sce. Default is 'counts'.",
      process = "estimation",
      function_name = c("fit_copula", "simu_new")
    ),
    param_character(
      id = "copula",
      default = "gaussian",
      alternatives = c("gaussian", "vine"),
      description = "A string of the copula choice. Must be one of 'gaussian' or 'vine'. Default is 'gaussian'. Note that vine copula may have better modeling of high-dimensions, but can be very slow when features are >1000.",
      process = "estimation",
      function_name = "fit_copula"
    ),
    param_others(
      id = "mean_mat",
      type = "matrix",
      process = "simulation",
      description = "A cell by feature matrix of the mean parameter.",
      function_name = "simu_new"
    ),
    param_others(
      id = "sigma_mat",
      type = "matrix",
      process = "simulation",
      description = "A cell by feature matrix of the sigma parameter.",
      function_name = "simu_new"
    ),
    param_others(
      id = "zero_mat",
      type = "matrix",
      process = "simulation",
      description = "A cell by feature matrix of the zero-inflation parameter.",
      function_name = "simu_new"
    ),
    param_others(
      id = "quantile_mat",
      type = "matrix",
      process = "simulation",
      description = "A cell by feature matrix of the multivariate quantile.",
      function_name = "simu_new"
    ),
    param_others(
      id = "copula_list",
      type = "list",
      process = "simulation",
      description = "A list of copulas for generating the multivariate quantile matrix. If provided, the quantile_mat must be NULL.",
      function_name = "simu_new"
    ),
    param_Boolean(
      id = "fastmvn",
      default = FALSE,
      description = "An logical variable. If TRUE, the sampling of multivariate Gaussian is done by mvnfast, otherwise by mvtnorm. Default is FALSE.",
      process = "simulation",
      function_name = "simu_new"
    ),
    param_Boolean(
      id = "nonnegative",
      default = TRUE,
      description = "A logical variable. If TRUE, values < 0 in the synthetic data will be converted to 0. Default is TRUE (since the expression matrix is nonnegative).",
      process = "simulation",
      function_name = "simu_new"
    ),
    param_Boolean(
      id = "nonzerovar",
      default = TRUE,
      description = "A logical variable. If TRUE, for any gene with zero variance, a cell will be replaced with 1. This is designed for avoiding potential errors, for example, PCA.",
      process = "simulation",
      function_name = "simu_new"
    )
  )

  scDesign3_method <- method_definition(
    method = "scDesign3",
    programming = "R",
    url = "https://github.com/SONGDONGYUAN1994/scDesign3/tree/main",
    authors = authors_definition(
      first = "Dongyuan",
      last = "Song",
      email = "dongyuansong@ucla.edu",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "scDesign3 generates realistic in silico data for multimodal single-cell and spatial omics",
      doi = "10.1038/s41587-023-01772-1",
      journal = "Nature Biotechnology",
      date = "2023",
      peer_review = TRUE
    ),
    description = "The R package scDesign3 is an all-in-one single-cell data simulation tool by using reference datasets with different cell states (cell types, trajectories or and spatial coordinates), different modalities (gene expression, chromatin accessibility, protein abundance, DNA methylation, etc), and complex experimental designs.")

  list(scDesign3_method = scDesign3_method,
       scDesign3_parameters = scDesign3_parameters)
}
