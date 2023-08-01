#' Estimate Parameters From Real Datasets by scDesign3
#'
#' This function is used to estimate useful parameters from a real dataset by
#' using `fit_marginal` and `fit_copula` function in scDesign3 package.
#'
#' @param ref_data A count matrix. Each row represents a gene and each column
#' represents a cell.
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @param verbose Logical.
#' @param seed An integer of a random seed.
#' @importFrom peakRAM peakRAM
#' @importFrom splatter splatEstimate
#' @importFrom SingleCellExperiment counts
#' @return A list contains the estimated parameters and the results of execution
#' detection.
#' @details
#' When you use scDesign3 to estimate parameters from a real dataset, you can optionally input
#' a numeric vector to specify the groups or plates that each cell comes from,
#' like `other_prior = list(group.condition = the numeric vector)`.
#'
#' You can also optionally input a numeric vector to specify the batches that each cell comes from,
#' like `other_prior = list(batch.condition = the numeric vector)`.
#'
#' See `Examples` and learn from it.
#' @export
#' @references
#' Song, D., Wang, Q., Yan, G. et al. scDesign3 generates realistic in silico data for multimodal single-cell and spatial omics. Nat Biotechnol (2023). https://doi.org/10.1038/s41587-023-01772-1
#'
#' Github URL: <https://github.com/SONGDONGYUAN1994/scDesign3/tree/main>
#' @examples
#' \dontrun{
#' ref_data <- simmethods::data
#' estimate_result <- simmethods::Splat_estimation(ref_data = data,
#'                                                 verbose = TRUE,
#'                                                 seed = 10)
#' estimate_result <- estimate_result[["estimate_result"]]
#' ## Check the class
#' class(estimate_result) == "SplatParams"
#' }
#'
scDesign3_estimation <- function(ref_data,
                                 other_prior = NULL,
                                 verbose = FALSE,
                                 seed
){
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  if(!requireNamespace("scDesign3", quietly = TRUE)){
    message("scDesign3 is not installed on your device...")
    message("Installing scDesign3...")
    devtools::install_github("SONGDONGYUAN1994/scDesign3")
  }

  ### group
  if(!is.null(other_prior[["group.condition"]])){
    col_data <- data.frame("cell_type" = paste0("Group", as.numeric(other_prior[["group.condition"]])))
  }else{
    col_data <- data.frame("cell_type" = rep("A", ncol(ref_data)))
  }
  rownames(col_data) <- colnames(ref_data)
  col_data$"cell_name" <- colnames(ref_data)
  ### batch
  if(!is.null(other_prior[["batch.condition"]])){
    col_data$batch.condition <- paste0("Batch", as.numeric(other_prior[["batch.condition"]]))
  }
  if(!is.null(other_prior[["spatial.x"]]) & !is.null(other_prior[["spatial.y"]])){
    col_data$spatial.x <- other_prior[["spatial.x"]]
    col_data$spatial.y <- other_prior[["spatial.y"]]
    spatial <- c("spatial.x", "spatial.y")
  }else{
    spatial <- NULL
  }
  ### trajectory
  if(!is.null(other_prior[["traj"]])){
    if(!requireNamespace("dyndimred", quietly = TRUE)){
      devtools::install_github("dynverse/dyndimred")
    }
    cat("Constructing lineages for the data...\n")
    traj_info <- pseudotime_info(ref_data = ref_data,
                                 other_prior = other_prior,
                                 col_data = col_data,
                                 seed = seed)
    col_data <- traj_info$col_data
    mu_formula <- traj_info$mu_formula
    pseudotime <- traj_info$pseudotime
  }else{
    pseudotime <- NULL
  }

  # Establish SingleCellExperiment
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = ref_data),
                                                         colData = col_data)
  ### step by step
  #### other_covariates
  other_covariates <- NULL
  if(!is.null(other_prior[["batch.condition"]])){
    other_covariates <- "batch.condition"
  }
  scDeisgn3_data <- scDesign3::construct_data(
    sce = sce,
    assay_use = "counts",
    celltype = "cell_type",
    pseudotime = pseudotime,
    spatial = spatial,
    other_covariates = other_covariates,
    corr_by = "cell_type",
    ncell = ifelse(!is.null(other_prior[["nCells"]]), other_prior[["nCells"]], dim(sce)[2])
  )

  ### mu_formula
  if(!is.null(other_prior[["mu_formula"]])){
    mu_formula <- other_prior[["mu_formula"]]
  }else{
    mu_formula <- NULL
    for(i in c("cell_type", "batch_condition")){
      if(!is.null(other_prior[[i]])){
        if(is.null(mu_formula)){
          mu_formula <- append(mu_formula, i)
        }else{
          mu_formula <- paste0(mu_formula, " + ", i)
        }
      }
    }
    if(!is.null(other_prior[["spatial.x"]]) & !is.null(other_prior[["spatial.y"]])){
      mu_formula <- "s(spatial.x, spatial.y, bs = 'gp')"
    }
    if(!is.null(other_prior[["traj"]])){
      mu_formula <- traj_info$mu_formula
    }
  }
  scDesign3_est <- function(){
    scDeisgn3_marginal <- scDesign3::fit_marginal(
      data = scDeisgn3_data,
      predictor = ifelse(!is.null(other_prior[["predictor"]]), other_prior[["predictor"]], "gene"),
      mu_formula = ifelse(!is.null(mu_formula), mu_formula, "1"),
      sigma_formula = ifelse(!is.null(other_prior[["sigma_formula"]]), other_prior[["sigma_formula"]], "1"),
      family_use = ifelse(!is.null(other_prior[["family_use"]]), other_prior[["family_use"]], "nb"),
      n_cores = ifelse(!is.null(other_prior[["n_cores"]]), other_prior[["n_cores"]], 1),
      usebam = FALSE
    )
    set.seed(seed)
    scDeisgn3_copula <- scDesign3::fit_copula(
      sce = sce,
      assay_use = "counts",
      marginal_list = scDeisgn3_marginal,
      family_use = ifelse(!is.null(other_prior[["family_use"]]), other_prior[["family_use"]], "nb"),
      copula = ifelse(!is.null(other_prior[["copula"]]), other_prior[["copula"]], "gaussian"),
      n_cores = ifelse(!is.null(other_prior[["n_cores"]]), other_prior[["n_cores"]], 1),
      new_covariate = NULL,
      input_data = scDeisgn3_data$dat
    )

    scDeisgn3_para <- scDesign3::extract_para(
      sce = sce,
      marginal_list = scDeisgn3_marginal,
      n_cores = 1,
      family_use = ifelse(!is.null(other_prior[["family_use"]]), other_prior[["family_use"]], "nb"),
      new_covariate = NULL,
      data = scDeisgn3_data$dat
    )
    estimate_result <- list(scDeisgn3_data = scDeisgn3_data,
                            scDeisgn3_marginal = scDeisgn3_marginal,
                            scDeisgn3_copula = scDeisgn3_copula,
                            scDeisgn3_para = scDeisgn3_para)
    return(estimate_result)

  }
  ##############################################################################
  ####                            Estimation                                 ###
  ##############################################################################
  if(verbose){
    message("Estimating parameters using scDesign3")
  }
  # Seed
  set.seed(seed)
  # Estimation
  estimate_detection <- peakRAM::peakRAM(
    estimate_result <- scDesign3_est()
  )
  estimate_result <- list(sce = sce,
                          scDeisgn3_data = estimate_result$scDeisgn3_data,
                          scDeisgn3_marginal = estimate_result$scDeisgn3_marginal,
                          scDeisgn3_copula = estimate_result$scDeisgn3_copula,
                          scDeisgn3_para = estimate_result$scDeisgn3_para,
                          group.condition = other_prior[["group.condition"]],
                          batch.condition = other_prior[["batch.condition"]])
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  estimate_output <- list(estimate_result = estimate_result,
                          estimate_detection = estimate_detection)
  return(estimate_output)
}



#' Simulate Datasets by scDesign3
#'
#' This function is used to simulate datasets by `simu_new` function in scDesign3 package.
#'
#' @param parameters A object generated by [scDesign3::fit_marginal()] and [scDesign3::fit_copula()]
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @param return_format A character. Alternatives choices: list, SingleCellExperiment,
#' Seurat, h5ad. If you select `h5ad`, you will get a path where the .h5ad file saves to.
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#' @export
#' @details
#' In addtion to simulate datasets with default parameters, users want to simulate
#' other kinds of datasets, e.g. a counts matrix with 2 or more cell groups. In
#' scDesign3, you can set extra parameters to simulate datasets.
#'
#' The customed parameters you can set are below:
#' 1. nCells. In scDesign3, you can set nCells directly `other_prior = list(nCells = 1000)` to simulate 1000 cells.
#' 2. nGroups. You can not directly set `other_prior = list(nGroups = 3)` to simulate 3 groups. Instead, you should specify the group labels for cells in the **estimation** step in `scDesign3_estimation` function.
#' 3. nBatches You can not directly set `other_prior = list(nBatches = 3)` to simulate 3 groups. Instead, you should specify the batch labels for cells in the **estimation** step in `scDesign3_estimation` function. Note that, if you customed another simulated cell number which is not equal to the one of real data, the batch information for simulated cells is not returned.
#'
#' For more customed parameters in scDesign3, please check [scDesign3::simu_new()].
#' @references
#' Song, D., Wang, Q., Yan, G. et al. scDesign3 generates realistic in silico data for multimodal single-cell and spatial omics. Nat Biotechnol (2023). https://doi.org/10.1038/s41587-023-01772-1
#'
#' Github URL: <https://github.com/SONGDONGYUAN1994/scDesign3/tree/main>
#' @examples
#' \dontrun{
#' ref_data <- simmethods::data
#'
#' ## Simulate datasets with default parameters
#' simulate_result <- scDesign3_simulation(ref_data = ref_data,
#'                                        return_format = "list",
#'                                        verbose = TRUE,
#'                                        seed = 111)
#' counts <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(counts)
#'
#'
#' ## Simulate two groups with 20% proportion of DEGs and 2 fold change. Note that
#' ## scDesign3 does not provide fold changes for genes so users would better set
#' ## fc.group parameter in simulation function.
#' simulate_result <- scDesign3_simulation(ref_data = ref_data,
#'                                        other_prior = list(nCells = 1000,
#'                                                           nGroups = 2,
#'                                                           de.prob = 0.2,
#'                                                           fc.group = 2),
#'                                        return_format = "list",
#'                                        verbose = TRUE,
#'                                        seed = 111)
#' counts <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(counts)
#' ## cell information
#' col_data <- simulate_result[["simulate_result"]][["col_meta"]]
#' table(col_data$group)
#' ## gene information
#' row_data <- simulate_result[["simulate_result"]][["row_meta"]]
#' table(row_data$de_gene)/4000
#' table(row_data$up_down)
#'
#'
#' ## Simulate three groups with 20% proportion of DEGs and 4 fold change. 20%, 40%
#' ## and 40% of cells belong to Group1, Group2 and Group3, respectively.
#' simulate_result <- scDesign3_simulation(ref_data = ref_data,
#'                                        other_prior = list(nCells = 1000,
#'                                                           nGroups = 3,
#'                                                           prob.group = c(0.2, 0.4, 0.4),
#'                                                           de.prob = 0.2,
#'                                                           fc.group = 4),
#'                                        return_format = "list",
#'                                        verbose = TRUE,
#'                                        seed = 111)
#' counts <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(counts)
#' ## cell information
#' col_data <- simulate_result[["simulate_result"]][["col_meta"]]
#' table(col_data$group)
#' }
#'
scDesign3_simulation <- function(parameters,
                                 other_prior = NULL,
                                 return_format,
                                 verbose = FALSE,
                                 seed
){
  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("scDesign3", quietly = TRUE)){
    message("scDesign3 is not installed on your device...")
    message("Installing scDesign3...")
    devtools::install_github("SONGDONGYUAN1994/scDesign3")
  }
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  sce <- parameters[["sce"]]
  scDeisgn3_data <-  parameters[["scDeisgn3_data"]]
  scDeisgn3_marginal <- parameters[["scDeisgn3_marginal"]]
  scDeisgn3_copula <- parameters[["scDeisgn3_copula"]]
  scDeisgn3_para <- parameters[["scDeisgn3_para"]]
  group.condition <- parameters[["group.condition"]]
  batch.condition <- parameters[["batch.condition"]]

  if(is.null(scDeisgn3_data[["newCovariate"]])){
    nCells <- dim(scDeisgn3_data[["count_mat"]])[1]
  }else{
    nCells <- dim(scDeisgn3_data[["newCovariate"]])[1]
  }
  if(!is.null(other_prior[["group.condition"]])){
    nGroups <- length(unique(as.numeric(other_prior[["group.condition"]])))
  }else{
    nGroups <- ifelse(is.null(group.condition), 1, length(unique(group.condition)))
  }
  if(!is.null(other_prior[["batch.condition"]])){
    nBatches <- length(unique(as.numeric(other_prior[["batch.condition"]])))
  }else{
    nBatches <- ifelse(is.null(batch.condition), 1, length(unique(batch.condition)))
  }
  # Return to users
  message(paste0("nCells: ", nCells))
  message(paste0("nGenes: ", dim(scDeisgn3_data[["count_mat"]])[2]))
  message(paste0("nGroups: ", nGroups))
  message(paste0("nBatches: ", nBatches))
  if(nBatches != 1){
    if(!is.null(scDeisgn3_data[["newCovariate"]])){
      cat("The number of simulated cells is not equal to the origin one and batch information can not be returned. \n")
    }
  }
  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    message("Simulating datasets using scDesign3")
  }
  # Seed
  set.seed(seed)
  # Estimation
  simulate_detection <- peakRAM::peakRAM(
    simulate_result <- scDesign3::simu_new(
      sce = sce,
      mean_mat = scDeisgn3_para$mean_mat,
      sigma_mat = scDeisgn3_para$sigma_mat,
      zero_mat = scDeisgn3_para$zero_mat,
      quantile_mat = NULL,
      copula_list = scDeisgn3_copula$copula_list,
      n_cores = ifelse(!is.null(other_prior[["n_cores"]]), other_prior[["n_cores"]], 1),
      family_use = ifelse(!is.null(other_prior[["family_use"]]), other_prior[["family_use"]], "nb"),
      input_data = scDeisgn3_data$dat,
      new_covariate = scDeisgn3_data$newCovariate,
      important_feature = scDeisgn3_copula$important_feature
    )
  )
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  counts <- as.matrix(simulate_result)
  col_data <- data.frame("cell_name" = colnames(counts))
  ### group
  if(nGroups != 1){
    if(!is.null(scDeisgn3_data[["newCovariate"]])){
      group <- as.character(scDeisgn3_data[["newCovariate"]][["corr_group"]])
      col_data$group <- group
    }else{
      group <- scDeisgn3_data[["dat"]][["corr_group"]]
      col_data$group <- group
    }
  }
  ### batch
  if(nBatches != 1){
    if(is.null(scDeisgn3_data[["newCovariate"]])){
      batch <- scDeisgn3_data[["dat"]][["batch.condition"]]
      col_data$batch <- batch
    }
  }
  ### row_data
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

