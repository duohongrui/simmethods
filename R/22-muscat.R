#' Estimate Parameters From Real Datasets by muscat
#'
#' This function is used to estimate useful parameters from a real dataset by
#' using `prepSim` function in muscat package.
#'
#' @param ref_data A count matrix. Each row represents a gene and each column
#' represents a cell.
#' @param verbose Logical.
#' @param seed An integer of a random seed.
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs and other variables are usually customed, so before simulating a dataset
#' you must point it out.
#' @importFrom peakRAM peakRAM
#' @return A list contains the estimated parameters and the results of execution
#' detection.
#' @export
#' @details
#' In muscat, cell group information is not neccessary but users can input it when
#' it is available by `other_prior = list(group.condition = xxx)`.
#'
#' For more information, see `Examples` and [muscat::prepSim()]
#'
#' @references
#' Crowell H L, Soneson C, Germain P L, et al. Muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. Nature communications, 2020, 11(1): 1-12. <https://doi.org/10.1038/s41467-020-19894-4>
#'
#' Github URL: <https://github.com/HelenaLC/muscat>
#'
#' @examples
#' ref_data <- simmethods::data
#' ## estimation
#' estimate_result <- simmethods::muscat_estimation(
#'   ref_data = ref_data,
#'   other_prior = NULL,
#'   verbose = TRUE,
#'   seed = 111
#' )
#'
#' ## cell groups
#' group_condition <- as.numeric(simmethods::group_condition)
#' ## estimation
#' estimate_result <- simmethods::muscat_estimation(
#'   ref_data = ref_data,
#'   other_prior = list(group.condition = group_condition),
#'   verbose = TRUE,
#'   seed = 111
#' )
muscat_estimation <- function(ref_data,
                              verbose = FALSE,
                              other_prior = NULL,
                              seed
){
  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("muscat", quietly = TRUE)){
    message("muscat is not installed on your device...")
    message("Installing muscat...")
    devtools::install_github("HelenaLC/muscat")
  }
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  if(!is.matrix(ref_data)){
    ref_data <- as.matrix(ref_data)
  }
  ## group
  if(is.null(other_prior[["group.condition"]])){
    other_prior[["group_id"]] <- rep("A", ncol(ref_data))
  }else{
    other_prior[["group_id"]] <- other_prior[["group.condition"]]
  }
  ## cluster
  if(is.null(other_prior[["cluster_id"]])){
    other_prior[["cluster_id"]] <- rep("A", ncol(ref_data))
  }

  col_data <- data.frame("sample_id" = other_prior[["cluster_id"]],
                         "group_id" = other_prior[["group_id"]],
                         "cluster_id" = other_prior[["cluster_id"]])

  ref_data <- SingleCellExperiment::SingleCellExperiment(list(counts = ref_data),
                                                         colData = col_data)

  ref_data <- muscat::prepSCE(x = ref_data)

  other_prior[["x"]] <- ref_data
  if(is.null(other_prior[["min_count"]])){
    other_prior[["min_count"]] <- 0
  }
  if(is.null(other_prior[["min_cells"]])){
    other_prior[["min_cells"]] <- 0
  }
  if(is.null(other_prior[["min_genes"]])){
    other_prior[["min_genes"]] <- 0
  }
  estimate_formals <- simutils::change_parameters(function_expr = "muscat::prepSim",
                                                  other_prior = other_prior,
                                                  step = "estimation")
  ##############################################################################
  ####                            Estimation                                 ###
  ##############################################################################
  if(verbose){
    message("Estimating parameters using muscat")
  }
  # Seed
  set.seed(seed)
  # Estimation
  estimate_detection <- peakRAM::peakRAM(
    estimate_result <- muscat::prepSim(x = estimate_formals[["x"]],
                                       min_count = estimate_formals[["min_count"]],
                                       min_cells = estimate_formals[["min_cells"]],
                                       min_genes = estimate_formals[["min_genes"]],
                                       min_size = NULL,
                                       group_keep = as.character(unique(other_prior[["group_id"]])),
                                       verbose = estimate_formals[["verbose"]])
  )
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  estimate_output <- list(estimate_result = estimate_result,
                          estimate_detection = estimate_detection)
  return(estimate_output)
}



#' Simulate Datasets by muscat
#'
#' This function is used to simulate datasets from learned parameters by `simData`
#' function in muscat package.
#'
#' @param parameters A object generated by [muscat::prepSim()]
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @param return_format A character. Alternative choices: list, SingleCellExperiment,
#' Seurat, h5ad. If you select `h5ad`, you will get a path where the .h5ad file saves to.
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#' @importFrom SummarizedExperiment rowData<-
#' @export
#' @details
#' In addtion to simulate datasets with default parameters, users want to simulate
#' other kinds of datasets, e.g. a counts matrix with 2 or more cell groups. In
#' muscat, you can set extra parameters to simulate datasets.
#'
#' The customed parameters you can set are below:
#' 1. nCells. In muscat, you can set nCells directly. For example, if you want to simulate 1000 cells, you can type `other_prior = list(nCells = 1000)`.
#' 2. nGenes. You can directly set `other_prior = list(nGenes = 5000)` to simulate 5000 genes.
#' 3. nGroups. In muscat, `nGroups` can be 1 or 2 because muscat can only simulate two cell groups.
#' 4. de.prob. You can directly set `other_prior = list(de.prob = 0.2)` to simulate DEGs that account for 20 percent of all genes.
#' 5. fc.group. You can directly set `other_prior = list(fc.group = 2)` to specify the minimum fold change of DEGs.
#'
#' For more customed parameters in muscat, please check [muscat::simData()].
#'
#' @references
#' Crowell H L, Soneson C, Germain P L, et al. Muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. Nature communications, 2020, 11(1): 1-12. <https://doi.org/10.1038/s41467-020-19894-4>
#'
#' Github URL: <https://github.com/HelenaLC/muscat>
#'
#' @examples
#' ref_data <- simmethods::data
#' ## cell groups
#' group_condition <- as.numeric(simmethods::group_condition)
#' ## estimation
#' estimate_result <- simmethods::muscat_estimation(
#'   ref_data = ref_data,
#'   other_prior = list(group.condition = group_condition),
#'   verbose = TRUE,
#'   seed = 111
#' )
#'
#' # 1) Simulate with default parameters
#' simulate_result <- simmethods::muscat_simulation(
#'   parameters = estimate_result[["estimate_result"]],
#'   other_prior = NULL,
#'   return_format = "list",
#'   verbose = TRUE,
#'   seed = 111
#' )
#' ## counts
#' counts <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(counts)
#'
#'
#' # 2) Simulate 1000 cells and 2000 genes
#' simulate_result <- simmethods::muscat_simulation(
#'   parameters = estimate_result[["estimate_result"]],
#'   other_prior = list(nCells = 1000,
#'                      nGenes = 2000),
#'   return_format = "list",
#'   verbose = TRUE,
#'   seed = 111
#' )
#'
#' ## counts
#' counts <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(counts)
#'
#'
#' # 3) Simulate 2 groups (20% proportion of DEGs, 4 fold change)
#' simulate_result <- simmethods::muscat_simulation(
#'   parameters = estimate_result[["estimate_result"]],
#'   other_prior = list(nCells = 1000,
#'                      nGenes = 2000,
#'                      nGroups = 2,
#'                      de.prob = 0.2,
#'                      fc.group = 4),
#'   return_format = "list",
#'   verbose = TRUE,
#'   seed = 111
#' )
#' ## cell information
#' col_data <- simulate_result[["simulate_result"]][["col_meta"]]
#' table(col_data$group)/1000
#' ## gene information
#' row_data <- simulate_result[["simulate_result"]][["row_meta"]]
#' table(row_data$de_gene)[2]/2000
muscat_simulation <- function(parameters,
                              other_prior = NULL,
                              return_format,
                              verbose = FALSE,
                              seed
){
  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("muscat", quietly = TRUE)){
    message("muscat is not installed on your device...")
    message("Installing muscat...")
    devtools::install_github("HelenaLC/muscat")
  }
  other_prior[["x"]] <- parameters
  ## nCells
  if(!is.null(other_prior[["nCells"]])){
    other_prior[["nc"]] <- other_prior[["nCells"]]
  }else{
    other_prior[["nc"]] <- ncol(parameters)
  }
  ## nGenes
  if(!is.null(other_prior[["nGenes"]])){
    other_prior[["ng"]] <- other_prior[["nGenes"]]
  }else{
    other_prior[["ng"]] <- nrow(parameters)
  }
  ## fc.group
  if(!is.null(other_prior[["fc.group"]])){
    other_prior[["lfc"]] <- log2(other_prior[["fc.group"]])
  }else{
    other_prior[["lfc"]] <- 1
  }
  ## de.prob
  if(!is.null(other_prior[["de.prob"]])){
    other_prior[["p_dd"]] <- c(1-other_prior[["de.prob"]],
                               0, other_prior[["de.prob"]], 0, 0, 0)
  }else{
    other_prior[["p_dd"]] <- c(0.9, 0, 0.1, 0, 0, 0)
  }
  ## nGroups
  if(is.null(other_prior[["nGroups"]])){
    other_prior[["nGroups"]] <- 1
  }
  if(other_prior[["nGroups"]] == 1){
    other_prior[["dd"]] <- FALSE
  }else{
    other_prior[["ns"]] <- 2
    other_prior[["dd"]] <- TRUE
    other_prior[["paired"]] <- TRUE
  }
  ## phylo_pars
  other_prior[["phylo_pars"]] <- c(0, 3)
  other_prior[["force"]] <- TRUE
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  simulate_formals <- simutils::change_parameters(function_expr = "muscat::simData",
                                                  other_prior = other_prior,
                                                  step = "simulation")
  nGroups <- ifelse(simulate_formals[["dd"]], 2, 1)
  de.group <- simulate_formals[["p_dd"]][3]/sum(simulate_formals[["p_dd"]])
  # Return to users
  message(glue::glue("nCells: {simulate_formals[['nc']]}"))
  message(glue::glue("nGenes: {simulate_formals[['ng']]}"))
  message(glue::glue("nGroups: {nGroups}"))
  message(glue::glue("de.group: {de.group}"))
  message(glue::glue("fc.group: {2^simulate_formals[['lfc']]}"))
  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    message("Simulating datasets using muscat")
  }
  # Seed
  set.seed(seed)
  # Simulation
  simulate_detection <- peakRAM::peakRAM(
    simulate_result <- do.call(muscat::simData, simulate_formals)
  )
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  ## counts
  counts <- SingleCellExperiment::counts(simulate_result)
  colnames(counts) <- paste0("Cell", 1:ncol(counts))
  rownames(counts) <- paste0("Gene", 1:nrow(counts))
  ## gene information
  gene_info <- metadata(simulate_result)$gene_info
  row_data <- gene_info %>%
    dplyr::transmute("gene_name" = rownames(counts),
                     "de_gene" = case_when(
                       gene_info$"category" == "ee" ~ "no",
                       gene_info$"category" == "de" ~ "yes"
                     ),
                     "fc_gene" = 2^gene_info$"logFC")
  rownames(row_data) <- row_data$gene_name
  ## cell information
  if(is.null(simulate_result$group_id)){
    group <- rep("Group1", ncol(counts))
  }else{
    group <- paste0("Group", as.numeric(simulate_result$group_id))
  }
  col_data <- data.frame("cell_gene" = colnames(counts),
                         "group" = group)
  rownames(col_data) <- col_data$cell_gene
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

