#' Get Information of scGAN
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' scGAN_method_definition <- scGAN_method_definition()
#'
scGAN_method_definition <- function(...){

  scGAN_parameters <- parameter_sets(
    param_reference(
      id = "ref_data",
      type = "matrix",
      default = NULL,
      process = "estimation",
      force = TRUE,
      description = "A numeric matrix with rows representing genes and columns representing cells. Gene names are given as row names.",
      function_name = "scGAN_estimation"
    ),
    param_integer(
      id = "GPU",
      default = 1L,
      process = "estimation",
      description = "The GPU number to use and execute the training step.",
      function_name = "scGAN_estimation"
    ),
    param_numeric(
      id = "res",
      process = "estimation",
      default = 1.5,
      lower = 0,
      description = "The resolution to perform louvain clustering.",
      function_name = "scGAN_estimation"
    ),
    param_integer(
      id = "min_cells",
      default = 3L,
      lower = 0L,
      process = "estimation",
      description = "Remove the genes which have expressions below the specified number of cells.",
      function_name = "scGAN_estimation"
    ),
    param_integer(
      id = "min_genes",
      default = 10L,
      lower = 0L,
      process = "estimation",
      description = "Remove cells which have expressions below in the specified number of genes.",
      function_name = "scGAN_estimation"
    ),
    param_character(
      id = "raw_input",
      process = "estimation",
      default = "/scGAN/input_data.h5ad",
      alternatives = "/scGAN/input_data.h5ad",
      description = "The path of reference data in the docker container.",
      function_name = "scGAN_estimation"
    ),
    param_character(
      id = "scale",
      process = "estimation",
      default = "normalize_per_cell_LS_20000",
      alternatives = "normalize_per_cell_LS_20000",
      description = "The normalization method for reference data.",
      function_name = "scGAN_estimation"
    ),
    param_Boolean(
      id = "balanced_split",
      default = TRUE,
      process = "estimation",
      description = "Whether to split the test and valid cells to all clusters equally or not.",
      function_name = "scGAN_estimation"
    ),
    param_character(
      id = "split_seed",
      process = "estimation",
      default = "default",
      alternatives = "default",
      description = "The seed when spilting data into valid and test datasets.",
      function_name = "scGAN_estimation"
    ),
    param_integer(
      id = "test_cells",
      default = 0L,
      lower = 0L,
      description = "The number of cells in test dataset.",
      process = "estimation",
      function_name = "scGAN_estimation"
    ),
    param_integer(
      id = "valid_cells",
      default = 2000L,
      lower = 0L,
      description = "The number of cells in valid dataset.",
      process = "estimation",
      function_name = "scGAN_estimation"
    ),
    param_integer(
      id = "max_steps",
      default = 1000000L,
      lower = 0L,
      description = "The Number of steps in the (outer) training loop.",
      process = "estimation",
      function_name = "scGAN_estimation"
    ),
    param_Boolean(
      id = "decay",
      default = TRUE,
      process = "estimation",
      description = "If True, uses an exponential decay of the learning rate.",
      function_name = "scGAN_estimation"
    ),
    param_numeric(
      id = "alpha_0",
      process = "estimation",
      default = 0.0001,
      lower = 0,
      description = "Initial learning rate value.",
      function_name = "scGAN_estimation"
    ),
    param_numeric(
      id = "alpha_final",
      process = "estimation",
      default = 0.00001,
      lower = 0,
      description = "Final value of the learning rate if the decay is set to True.",
      function_name = "scGAN_estimation"
    ),
    param_character(
      id = "algorithm",
      process = "estimation",
      default = "AMSGrad",
      alternatives = c("AMSGrad", "Adam"),
      description = "Optimizer used in the training..",
      function_name = "scGAN_estimation"
    ),
    param_numeric(
      id = "beta1",
      process = "estimation",
      default = 0.5,
      lower = 0,
      description = "Exponential decay for the first-moment estimates.",
      function_name = "scGAN_estimation"
    ),
    param_numeric(
      id = "beta2",
      process = "estimation",
      default = 0.9,
      lower = 0,
      description = "Exponential decay for the second-moment estimates.",
      function_name = "scGAN_estimation"
    ),
    param_integer(
      id = "batch_size",
      default = 128L,
      lower = 0L,
      description = "Batch size used for the training.",
      process = "estimation",
      function_name = "scGAN_estimation"
    ),
    param_integer(
      id = "critic_iter",
      default = 5L,
      lower = 0L,
      description = "Number of training iterations of the critic (inner loop) for each iteration on the generator (outer loop).",
      process = "estimation",
      function_name = "scGAN_estimation"
    ),
    param_integer(
      id = "progress_freq",
      default = 10L,
      lower = 0L,
      description = "Period (in steps) between displays of the losses values on the standard output.",
      process = "estimation",
      function_name = "scGAN_estimation"
    ),
    param_integer(
      id = "validation_freq",
      default = 1000L,
      lower = 0L,
      description = "The frequency in steps for validation (e.g. running t-SNE plots).",
      process = "estimation",
      function_name = "scGAN_estimation"
    ),
    param_integer(
      id = "save_freq",
      default = 1000L,
      lower = 0L,
      description = "Period (in steps) between saves of the model.",
      process = "estimation",
      function_name = "scGAN_estimation"
    ),
    param_integer(
      id = "summary_freq",
      default = 50L,
      lower = 0L,
      description = "Period (in steps) between logs for the Tensorboard.",
      process = "estimation",
      function_name = "scGAN_estimation"
    ),
    param_integer(
      id = "latent_dim",
      default = 128L,
      lower = 0L,
      description = "Dimension of the latent space used from which the input noise of the generator is sampled.",
      process = "estimation",
      function_name = "scGAN_estimation"
    ),
    param_integer(
      id = "output_LSN",
      default = 20000L,
      lower = 0L,
      description = "Parameter of the LSN layer at the output of the critic (i.e. total number of counts per generated cell). If set to None, the layer won't be added in the generator.",
      process = "estimation",
      function_name = "scGAN_estimation"
    ),
    param_others(
      id = "critic_layers",
      type = "list",
      default = list(1024, 512, 256),
      description = "List of integers corresponding to the number of neurons of each layer of the critic.",
      process = "estimation",
      function_name = "scGAN_estimation"
    ),
    param_others(
      id = "gen_layers",
      type = "list",
      default = list(256, 512, 1024),
      description = "List of integers corresponding to the number of neurons of each layer of the generator.",
      process = "estimation",
      function_name = "scGAN_estimation"
    ),
    param_numeric(
      id = "lambd",
      default = 10,
      lower = 0,
      description = "Regularization hyper-parameter to be used with the gradient penalty in the WGAN loss.",
      function_name = "scGAN_estimation"
    )
  )

  scGAN_method <- method_definition(
    method = "scGAN",
    programming = "Python",
    url = "https://github.com/imsb-uke/scGAN",
    authors = authors_definition(
      first = "Mohamed",
      last = "Marouf",
      email = "fabian.hausmann@zmnh.uni-hamburg.de",
      github = "https://github.com/imsb-uke/scGAN",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "Realistic in silico generation and augmentation of single-cell RNA-seq data using generative adversarial networks",
      doi = "10.1038/s41467-019-14018-z",
      journal = "Nature Communications",
      date = "2020",
      peer_review = TRUE
    ),
    description = "A statistical simulator for rational scRNA-seq experimental design.")

  list(scGAN_method = scGAN_method,
       scGAN_parameters = scGAN_parameters)
}
