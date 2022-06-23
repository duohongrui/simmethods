#' Get Information of BEARscc
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' BEARscc_method_definition <- BEARscc_method_definition()
#'
BEARscc_method_definition <- function(...){

  BEARscc_parameters <- parameter_sets(
    param_reference(
      id = "SCEList",
      type = "SingleCellExperiment",
      default = NULL,
      process = c("estimation", "simulation"),
      force = TRUE,
      description = "A SingleCellExperiment object that must contain the observed counts matrix as 'observed_expression' in assays, and must have the relevant spike-in samples identified using isSpike() as well as contain the expected actual concentrations of these spike-ins as spikeConcentrations in metadata. Please see the vignette for more detail about constructing the appropriate SCEList.",
      function_name = c("estimate_noiseparameters", "simulate_replicates")
    ),
    param_Boolean(
      id = "plot",
      default = FALSE,
      process = "estimation",
      description = "When plot=TRUE produces plots to investigate quality of data fits with root file name set by file option.",
      function_name = "estimate_noiseparameters"
    ),
    param_numeric(
      id = "sd_inflate",
      default = 0,
      process = "estimation",
      description = "An optional parameter to modulate the estimated noise. The estimated standard deviation of spike-ins can be scaled by this factor. We recommend leaving the value at the default of 0.",
      function_name = "estimate_noiseparameters"
    ),
    param_numeric(
      id = "bins",
      default = 10,
      process = "estimation",
      description = "The parameter determines the number of bins for comparison of the quality of fit between the mixed-model and observed data for each spike-in alpha in order to calculate the relationship between alpha and mean in the noise model. This should be set lower for small datasets and higher for datasets with more observations.",
      function_name = "estimate_noiseparameters"
    ),
    param_numeric(
      id = "max_cumprob",
      force = TRUE,
      default = 0.9999,
      process = c("estimation", "simulation"),
      description = "Because a cumulative distribution will range from n=0 to a countable infinity, the event space needs to be set to cover a reasonable fraction of the probability density. This parameter determines the the fraction of probability density covered by the event space, which in turn defines the highes count number in the event space. We recommend users use the default value of 0.9999.",
      function_name = c("estimate_noiseparameters", "simulate_replicates")
    ),
    param_Boolean(
      id = "write.noise.model",
      process = "estimation",
      default = TRUE,
      description = "When write.noise.model=TRUE outputs two tab-delimited files containing the dropout effects and noise model parameters; this allows users to apply the noise generation on a seperate high compute node. The root file name is set by file option.",
      function_name = "estimate_noiseparameters"
    ),
    param_character(
      id = "file",
      default = "noise_estimation",
      alternatives = c("noise_estimation"),
      process = "estimation",
      description = "Describes the root name for files written out by write.noise.model and plot options.",
      function_name = "estimate_noiseparameters"
    ),
    param_numeric(
      id = "dropout_inflate",
      default = 1,
      process = "estimation",
      description = "A scaling parameter for increasing explicitly the number of drop-outs present beyond those estimated by spike-ins. The value must be greater than 0 or an error will occur. Values below one will diminish drop-outs in simulated replicates, and values above one will increase drop-outs in simulated replicates. We recommend users use the default value of 1.",
      function_name = "estimate_noiseparameters"
    ),
    param_character(
      id = "model_view",
      default = c("Observed", "Optimized"),
      alternatives = c("Observed", "Optimized", "Poisson", "Neg. Binomial"),
      process = "estimation",
      description = "model_view=c('Observed', 'Optimized', 'Poisson', 'Neg. Binomial', determines the statistical distributions that should be plotted for the ERCC plots output by plot=TRUE.",
      function_name = "estimate_noiseparameters"
    ),
    param_numeric(
      id = "alpha_resolution",
      default = 0.005,
      process = "estimation",
      description = "Because the alpha parameter is enumerated discretely and empirically evaluated for each value for each spike-in, it is necessary to specify the resolution (how small the step is between each explicit alpha test); this parameter defines the resolution of alpha values tested for maximum empirical fit to spike-ins. It is recommended that users utilize the default resolution.",
      function_name = "estimate_noiseparameters"
    ),
    param_character(
      id = "tie_function",
      default = "maximum",
      alternatives = c("minimum", "maximum"),
      process = "estimation",
      description = "The parameter tie_function=c('minimum', 'maximum') tells BEARscc how to handle a tie alpha value for fitting the mixture model to an individual spike-in. If maximum, then BEARscc will chose the maximum alpha value with the best fit; conversely, if minimum is set, then BEARscc will choose the minimum alpha value with the best fit.",
      function_name = "estimate_noiseparameters"
    ),
    param_integer(
      id = "n",
      default = 3L,
      description = "The number of simulated technical replicates to generate.",
      process = "simulation",
      function_name = "simulate_replicates"
    )
  )

  BEARscc_method <- method_definition(
    method = "BEARscc",
    programming = "R",
    url = "https://www.bioconductor.org/packages/release/bioc/html/BEARscc.html",
    authors = authors_definition(
      first = "David T.",
      last = "Severson",
      email = "benjamin.schuster-boeckler@ludwig.ox.ac.uk",
      github = "https://github.com/seversond12/BEARscc",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "BEARscc determines robustness of single-cell clusters using simulated technical replicates",
      doi = "10.1038/s41467-018-03608-y",
      journal = "Nature Communications",
      date = "2018",
      peer_review = TRUE
    ),
    description = "BEARscc is a noise estimation and injection tool that is designed to assess putative single-cell RNA-seq clusters in the context of experimental noise estimated by ERCC spike-in controls.")

  list(BEARscc_method = BEARscc_method,
       BEARscc_parameters = BEARscc_parameters)
}
