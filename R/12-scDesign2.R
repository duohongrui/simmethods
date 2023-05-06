#' Estimate Parameters From Real Datasets by scDesign2
#'
#' This function is used to estimate useful parameters from a real dataset by
#' using `fit_model_scDesign2` function in scDesign2 package.
#'
#' @param ref_data A count matrix. Each row represents a gene and each column
#' represents a cell.
#' @param verbose Logical.
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @param seed An integer of a random seed.
#'
#' @return A list contains the estimated parameters and the results of execution
#' detection.
#' @export
#' @details
#' In scDesign2, users can input cell group or cell type information before estimating
#' parameters from real datasets. For more instructions or information, see `Examples`
#' or [scDesign2::fit_model_scDesign2()]
#' @references
#' Sun T, Song D, Li W V, et al. scDesign2: a transparent simulator that generates high-fidelity single-cell gene expression count data with gene correlations captured. Genome biology, 2021, 22(1): 1-37. <https://doi.org/10.1186/s13059-021-02367-2>
#'
#' Github URL: <https://github.com/JSB-UCLA/scDesign2>
#' @examples
#' \dontrun{
#' ref_data <- simmethods::data
#'
#' ## cell groups
#' group_condition <- as.numeric(simmethods::group_condition)
#' ## In scDesign2, cell group information is not neccessary which indicates the type
#' ## that each cell belongs to.
#' estimate_result <- simmethods::scDesign2_estimation(
#'   ref_data = ref_data,
#'   other_prior = list(group.condition = group_condition),
#'   verbose = TRUE,
#'   seed = 111
#' )
#' ## Or you can input information of cell types via cell_type_sel parameter described
#' ## in scDesign2::fit_model_scDesign2 function
#' estimate_result <- simmethods::scDesign2_estimation(
#'   ref_data = ref_data,
#'   other_prior = list(cell_type_sel = paste0("cell_type",
#'                                             group_condition)),
#'   verbose = TRUE,
#'   seed = 111
#' )
#' }
#'
scDesign2_estimation <- function(ref_data,
                                 verbose = FALSE,
                                 other_prior = NULL,
                                 seed
){
  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("scDesign2", quietly = TRUE)){
    message("scDesign2 is not installed on your device...")
    message("Installing scDesign2...")
    devtools::install_github("JSB-UCLA/scDesign2")
  }
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  if(!is.matrix(ref_data)){
    ref_data <- as.matrix(ref_data)
  }
  if(!is.null(other_prior[["group.condition"]])){
    fake_cell_name <- paste0(LETTERS[other_prior[["group.condition"]]], "cell")
    colnames(ref_data) <- fake_cell_name
    other_prior[["cell_type_sel"]] <- unique(fake_cell_name)
  }else{
    if(!is.null(other_prior[["cell_type_sel"]])){
      colnames(ref_data) <- other_prior[["cell_type_sel"]]
      other_prior[["cell_type_sel"]] <- unique(other_prior[["cell_type_sel"]])
    }else{
      colnames(ref_data) <- rep("Acell", ncol(ref_data))
      other_prior[["cell_type_sel"]] <- "Acell"
    }
  }
  other_prior[["data_mat"]] <- ref_data
  if(is.null(other_prior[["sim_method"]])){
    other_prior[["sim_method"]] <- "copula"
  }
  if(is.null(other_prior[["marginal"]])){
    other_prior[["marginal"]] <- "auto_choose"
  }
  estimate_formals <- simutils::change_parameters(function_expr = "scDesign2::fit_model_scDesign2",
                                                  other_prior = other_prior,
                                                  step = "estimation")
  ##############################################################################
  ####                            Estimation                                 ###
  ##############################################################################
  if(verbose){
    message("Estimating parameters using scDesign2")
  }
  # Seed
  set.seed(seed)
  # Estimation
  estimate_detection <- peakRAM::peakRAM(
    estimate_result <- do.call(scDesign2::fit_model_scDesign2, estimate_formals)
  )
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  estimate_output <- list(estimate_result = estimate_result,
                          estimate_detection = estimate_detection)
  return(estimate_output)
}



#' Simulate Datasets by scDesign2
#'
#' This function is used to simulate datasets from learned parameters by `simulate_count_scDesign2`
#' function in scDesign2 package.
#'
#' @param parameters A object generated by [scDesign2::fit_model_scDesign2()]
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @param return_format A character. Alternative choices: list, SingleCellExperiment,
#' Seurat, h5ad. If you select `h5ad`, you will get a path where the .h5ad file saves to.
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#' @importFrom BiocGenerics sapply grep
#' @export
#' @details
#' In addtion to simulate datasets with default parameters, users want to simulate
#' other kinds of datasets, e.g. a counts matrix with 2 or more cell groups. In
#' scDesign2, you can set extra parameters to simulate datasets.
#'
#' The customed parameters you can set are below:
#' 1. nCells. In scDesign2, You can directly set `other_prior = list(nCells = 5000)` to simulate 5000 cells.
#' 2. prob.group. You can directly set `other_prior = list(prob.group = c(0.4, 0.6))` to assign two proportions of cell groups. Note that the the length of the vector must equal to the groups or types that users input in estimation step.
#'
#' For more customed parameters in scDesign2, please check [scDesign2::simulate_count_scDesign2()].
#'
#' @references
#' Sun T, Song D, Li W V, et al. scDesign2: a transparent simulator that generates high-fidelity single-cell gene expression count data with gene correlations captured. Genome biology, 2021, 22(1): 1-37. <https://doi.org/10.1186/s13059-021-02367-2>
#'
#' Github URL: <https://github.com/JSB-UCLA/scDesign2>
#'
#' @examples
#' \dontrun{
#' ref_data <- simmethods::data
#'
#' ## cell groups
#' group_condition <- as.numeric(simmethods::group_condition)
#' ## In scDesign2, cell group information is not neccessary which indicates the type
#' ## that each cell belongs to.
#' estimate_result <- simmethods::scDesign2_estimation(
#'   ref_data = ref_data,
#'   other_prior = list(group.condition = group_condition),
#'   verbose = TRUE,
#'   seed = 111
#' )
#'
#' ## Simulate 1000 cells (40% in Group1, 60% in Group2)
#' simulate_result <- simmethods::scDesign2_simulation(
#'   parameters = estimate_result[["estimate_result"]],
#'   other_prior = list(nCells = 1000,
#'                      prob.group = c(0.4, 0.6)),
#'   return_format = "list",
#'   verbose = TRUE,
#'   seed = 111
#' )
#' ## counts
#' counts <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(counts)
#' ## cell information
#' col_data <- simulate_result[["simulate_result"]][["col_meta"]]
#' table(col_data$group)
#' }
#'
scDesign2_simulation <- function(parameters,
                                 other_prior = NULL,
                                 return_format,
                                 verbose = FALSE,
                                 seed
){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("scDesign2", quietly = TRUE)){
    message("scDesign2 is not installed on your device...")
    message("Installing scDesign2...")
    devtools::install_github("JSB-UCLA/scDesign2")
  }
  other_prior[["model_params"]] <- parameters
  ## nCells
  if(!is.null(other_prior[["nCells"]])){
    other_prior[["n_cell_new"]] <- other_prior[["nCells"]]
  }else{
    other_prior[["n_cell_new"]] <- sum(sapply(parameters, function(x){x[["n_cell"]]}))
  }
  ## prob.group
  if(is.null(other_prior[["prob.group"]])){
    other_prior[["cell_type_prop"]] <- rep(1/length(parameters), length(parameters))
  }else{
    other_prior[["cell_type_prop"]] <- other_prior[["prob.group"]]
  }

  if(is.null(other_prior[["sim_method"]])){
    other_prior[["sim_method"]] <- "copula"
  }
  if(is.null(other_prior[["reseq_method"]])){
    other_prior[["reseq_method"]] <- "mean_scale"
  }

  simulate_formals <- simutils::change_parameters(function_expr = "scDesign2::simulate_count_scDesign2",
                                                  other_prior = other_prior,
                                                  step = "simulation")

  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  # Return to users
  nGenes <- sum(BiocGenerics::sapply(BiocGenerics::grep(names(parameters[[1]]), pattern = "^gene_sel"), function(x){length(parameters[[1]][[x]])}))
  message(paste0("nCells: ", simulate_formals[['n_cell_new']]))
  message(paste0("nGenes: ", nGenes))
  message(paste0("nGroups: ", length(parameters)))
  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    message("Simulating datasets using scDesign2")
  }
  # Seed
  set.seed(seed)
  # Estimation
  simulate_detection <- peakRAM::peakRAM(
    simulate_result <- BiocGenerics::do.call(scDesign2::simulate_count_scDesign2, simulate_formals)
  )
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  counts <- simulate_result
  # Rename
  rownames(counts) <- paste0("Gene", 1:nrow(counts))
  colnames(counts) <- paste0("Cell", 1:ncol(counts))
  # col_data
  col_data <- data.frame("cell_name" = colnames(counts),
                         "group" = paste0("Group", rep(1:length(table(colnames(simulate_result))),
                                                       table(colnames(simulate_result)))))
  # Row data
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

