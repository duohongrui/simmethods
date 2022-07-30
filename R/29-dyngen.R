#' Estimate Parameters From Real Datasets by dyngen
#'
#' This function is used to estimate useful parameters from a real dataset by
#' using \code{infer_trajectory} function in dynwrap package.
#'
#' @param ref_data A count matrix. Each row represents a gene and each column
#' represents a cell.
#' @param verbose Logical.
#' @param seed An integer of a random seed.
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @importFrom dynwrap infer_trajectory wrap_expression add_grouping
#' @importFrom NbClust NbClust
#' @importFrom tislingshot ti_slingshot
#' @return A list contains the estimated parameters and the results of execution
#' detection.
#' @export
#' @details
#' In dyngen, users can input cell group information if it is available. If cell
#' group information is not provided, the procedure will detect cell groups by
#' kmeans automatically.
#' See `Examples` for more instructions.
#'
#' @references
#' Cannoodt R, Saelens W, Deconinck L, et al. Spearheading future omics analyses using dyngen, a multi-modal simulator of single cells. Nature Communications, 2021, 12(1): 1-9. <https://doi.org/10.1038/s41467-021-24152-2>
#'
#' CRAN URL: <https://cran.r-project.org/web/packages/dyngen/index.html>
#'
#' Github URL: <https://github.com/dynverse/dyngen>
#' @examples
#' ref_data <- simmethods::data
#'
#' estimate_result <- simmethods::dyngen_estimation(
#'   ref_data = ref_data,
#'   other_prior = NULL,
#'   verbose = TRUE,
#'   seed = 111
#' )
#'
#' ## estimation with cell group information
#' group_condition <- paste0("Group", as.numeric(simmethods::group_condition))
#' estimate_result <- simmethods::dyngen_estimation(
#'   ref_data = ref_data,
#'   other_prior = list(group.condition = group_condition),
#'   verbose = TRUE,
#'   seed = 111
#' )
dyngen_estimation <- function(ref_data,
                              verbose = FALSE,
                              other_prior,
                              seed){
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  if(!is.matrix(ref_data)){
    ref_data <- as.matrix(ref_data)
  }
  if(is.null(other_prior[["group.condition"]])){
    message("Performing k-means and determin the best number of clusters...")
    clust <- NbClust::NbClust(data = t(ref_data),
                              distance = 'euclidean',
                              min.nc = 2,
                              max.nc = sqrt(nrow(t(ref_data))),
                              method = "kmeans",
                              index = "dunn")
    other_prior[["group.condition"]] <- paste0("Group", clust[["Best.partition"]])
  }
  ref_data <- dynwrap::wrap_expression(counts = t(ref_data),
                                       expression = log2(t(ref_data) + 1))
  ## Add group
  if(verbose){
    message("Add grouping to data...")
  }
  ref_data <- dynwrap::add_grouping(dataset = ref_data,
                                    grouping = other_prior[["group.condition"]])
  ##############################################################################
  ####                            Estimation                                 ###
  ##############################################################################
  if(verbose){
    cat("Estimating parameters using dyngen\n")
  }
  # Estimation
  tryCatch({
    # Estimate parameters from real data and return parameters and detection results
    estimate_detection <- peakRAM::peakRAM(
      estimate_result <- dynwrap::infer_trajectory(dataset = ref_data,
                                                   method = tislingshot::ti_slingshot(),
                                                   parameters = NULL,
                                                   give_priors = NULL,
                                                   seed = seed,
                                                   verbose = verbose)
    )
  }, error = function(e){
    as.character(e)
  })
  estimate_result <- list(estimate_result = estimate_result,
                          data_dim = c(length(ref_data[["feature_ids"]]),
                                       length(ref_data[["cell_ids"]])))
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  estimate_output <- list(estimate_result = estimate_result,
                          estimate_detection = estimate_detection)
  return(estimate_output)
}


#' Simulate Datasets by dyngen
#'
#' @param parameters A object generated by \code{\link[dynwrap]{infer_trajectory}}
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @param return_format A character. Alternative choices: list, SingleCellExperiment,
#' Seurat, h5ad. If you select `h5ad`, you will get a path where the .h5ad file saves to.
#' @param verbose Logical. Whether to return messages or not.
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#' @import dyngen
#' @importFrom tools R_user_dir
#' @export
#' @details
#' In dyngen, users can only set `nCells` and `nGenes` to specify the number of cells and genes in the
#' simulated dataset. See `Examples` for instructions.
#'
#' @references
#' Cannoodt R, Saelens W, Deconinck L, et al. Spearheading future omics analyses using dyngen, a multi-modal simulator of single cells. Nature Communications, 2021, 12(1): 1-9. <https://doi.org/10.1038/s41467-021-24152-2>
#'
#' CRAN URL: <https://cran.r-project.org/web/packages/dyngen/index.html>
#'
#' Github URL: <https://github.com/dynverse/dyngen>
#'
#' @examples
#' ref_data <- simmethods::data
#'
#' ## estimation with cell group information
#' group_condition <- paste0("Group", as.numeric(simmethods::group_condition))
#' estimate_result <- simmethods::dyngen_estimation(
#'   ref_data = ref_data,
#'   other_prior = list(group.condition = group_condition),
#'   verbose = TRUE,
#'   seed = 111
#' )
#'
#' # # 1) Simulate with default parameters (need a lot of memory)
#' # simulate_result <- simmethods::dyngen_simulation(
#' #   parameters = estimate_result[["estimate_result"]],
#' #   other_prior = NULL,
#' #   return_format = "list",
#' #   verbose = TRUE,
#' #   seed = 111
#' # )
#' # ## counts
#' # counts <- simulate_result[["simulate_result"]][["count_data"]]
#' # dim(counts)
#'
#' # 2) 100 cells and 100 genes
#' # simulate_result <- simmethods::dyngen_simulation(
#' #   parameters = estimate_result[["estimate_result"]],
#' #   other_prior = list(nCells = 100,
#' #                      nGenes = 100),
#' #   return_format = "list",
#' #   verbose = TRUE,
#' #   seed = 111
#' # )
#'
#' ## counts
#' # counts <- simulate_result[["simulate_result"]][["count_data"]]
#' # dim(counts)
dyngen_simulation <- function(parameters,
                              other_prior = NULL,
                              return_format,
                              verbose = FALSE,
                              seed
){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("dyntoy", quietly = TRUE)){
    cat("dyntoy is not installed on your device\n")
    cat("Installing dyntoy...\n")
    devtools::install_github("dynverse/dyngen")
  }
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  true_topology <- parameters[["estimate_result"]][["trajectory_type"]]
  if(true_topology=='linear'){
    backbone <- dyngen::backbone_linear()
  }
  if(true_topology=='bifurcation'){
    backbone <- dyngen::backbone_bifurcating()
  }
  if(true_topology=='cycle'){
    backbone <- dyngen::backbone_cycle()
  }
  if(true_topology=='tree'){
    backbone <- dyngen::backbone_binary_tree()
  }
  if(true_topology=='multifurcation'){
    backbone <- dyngen::backbone_trifurcating()
  }
  # nCells
  if(!is.null(other_prior[["nCells"]])){
    nCells <- other_prior[["nCells"]]
  }else{
    nCells <- parameters[["data_dim"]][2]
  }
  # nGenes
  if(!is.null(other_prior[["nGenes"]])){
    nGenes <- other_prior[["nGenes"]]
  }else{
    nGenes <- parameters[["data_dim"]][1]
  }

  # Return to users
  cat(glue::glue("nCells: {nCells}"), "\n")
  cat(glue::glue("nGenes: {nGenes}"), "\n")
  # TFs and HKs
  num_tfs <- nrow(backbone$module_info)
  num_targets <- round((nGenes - num_tfs) / 2)
  num_hks <- nGenes - num_targets - num_tfs
  # Seed
  set.seed(seed)
  # Preparation
  init <- dyngen::initialise_model(backbone = backbone,
                                   num_cells = nCells,
                                   num_tfs = num_tfs,
                                   num_targets = num_targets,
                                   num_hks = num_hks,
                                   download_cache_dir = tools::R_user_dir("dyngen", "data"),
                                   simulation_params = dyngen::simulation_default(
                                     census_interval = 0.01,
                                     experiment_params = dyngen::simulation_type_wild_type(num_simulations = 1)
                                     ),
                                   verbose = TRUE)

  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    cat("Simulating datasets using dyngen\n")
  }
  # Seed
  set.seed(seed)
  # Estimation
  tryCatch({
    simulate_detection <- peakRAM::peakRAM(
      simulate_result <- dyngen::generate_dataset(init, make_plots = FALSE)
    )
  }, error = function(e){
    as.character(e)
  })
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  counts <- t(as.matrix(simulate_result[["dataset"]][["counts"]]))
  colnames(counts) <- paste0("Cell", 1:ncol(counts))
  rownames(counts) <- paste0("Gene", 1:nrow(counts))
  ## col_data
  col_data <- data.frame("cell_name" = colnames(counts))
  ## row_data
  row_data <- data.frame("gene_name" = rownames(counts))
  # Establish SingleCellExperiment
  simulate_result <- SingleCellExperiment::SingleCellExperiment(list(counts = counts),
                                                                colData = col_data,
                                                                rowData = row_data)
  simulate_result <- simutils::data_conversion(SCE_object = simulate_result,
                                               return_format = return_format)

  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  simulate_output <- list(simulate_result = simulate_result,
                          simulate_detection = simulate_detection)
  return(simulate_output)
}
