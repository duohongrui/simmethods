#' Estimate Parameters From Real Datasets by SPARSim
#'
#' This function is used to estimate useful parameters from a real dataset by
#' using \code{SPARSim_estimate_parameter_from_data} function in SPARSim package.
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
#' @importFrom SPARSim SPARSim_estimate_parameter_from_data scran_normalization
#' @importFrom stats model.matrix
#' @importFrom scater normalizeCounts
#' @return A list contains the estimated parameters and the results of execution
#' detection.
#' @export
#' @details
#' In SPARSim, the information of cell group condition can be input if neccessary
#' by `other_prior = list(group.condition = xxx)`. Note that the cell group condition
#' must be an integer vactor (e.g. 1, 2, 3, ...) to specify which condition that
#' each cell belongs to. See `Examples` below for more.
#'
#' @references
#' Baruzzo G, Patuzzi I, Di Camillo B. SPARSim single cell: a count data simulator for scRNA-seq data. Bioinformatics, 2020, 36(5): 1468-1475. <https://doi.org/10.1093/bioinformatics/btz752>
#'
#' Gitlab URL: <https://gitlab.com/sysbiobig/sparsim>
#'
#' @examples
#' ref_data <- simmethods::data
#'
#' # 1) Estimation without cell group information
#' estimate_result <- simmethods::SPARSim_estimation(
#'   ref_data = ref_data,
#'   other_prior = NULL,
#'   verbose = TRUE,
#'   seed = 111
#' )
#'
#' # 2) Estimation with cell group information (Note that an integer vector to specify
#' # which condition that each cell belongs to)
#' group_condition <- as.numeric(simmethods::group_condition)
#' estimate_result <- simmethods::SPARSim_estimation(
#'   ref_data = ref_data,
#'   other_prior = list(group.condition = group_condition),
#'   verbose = TRUE,
#'   seed = 111
#' )
#'
#' # 3) Users can also utilize spike-in genes to estimate parameters. In this case, users
#' ## must input dilution.factor and volume (nanoliter) parameters. Note that the
#' ## reference matrix must contain spike-in gene counts.
#' ref_data <- simmethods::data
#'
#' group_condition <- as.numeric(simmethods::group_condition)
#' estimate_result <- simmethods::SPARSim_estimation(
#'   ref_data = ref_data,
#'   other_prior = list(group.condition = group_condition,
#'                      dilution.factor = 50000,
#'                      volume = 0.01),
#'   verbose = TRUE,
#'   seed = 111
#' )
#' ## check spike-in parameters
#' spikein_params <- estimate_result[["estimate_result"]][["SPARSim_spikein_parameter"]]
SPARSim_estimation <- function(ref_data,
                               verbose = FALSE,
                               other_prior = NULL,
                               seed
){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("SPARSim", quietly = TRUE)){
    cat("SPARSim is not installed on your device\n")
    cat("Installing SPARSim...\n")
    devtools::install_gitlab("sysbiobig/sparsim")
  }
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  if(!is.matrix(ref_data)){
    ref_data <- as.matrix(ref_data)
  }
  colnames(ref_data) <- paste0("Cell", 1:ncol(ref_data))
  other_prior[["raw_data"]] <- ref_data
  other_prior[["norm_data"]] <- scater::normalizeCounts(ref_data)

  if(is.null(other_prior[["group.condition"]])){
    count_matrix_conditions <- list(conditionA = 1:ncol(ref_data))
    other_prior[["conditions"]] <- count_matrix_conditions
  }else{
    count_matrix_conditions <- list()
    condition_length <- length(unique(other_prior[["group.condition"]]))
    for(i in 1:condition_length){
      index <- which(other_prior[["group.condition"]] == i)
      count_matrix_conditions[[paste0("cond_", LETTERS[i])]] <- index
    }
    other_prior[["conditions"]] <- count_matrix_conditions
  }

  if(!is.null(other_prior[["dilution.factor"]]) & !is.null(other_prior[["volume"]])){
    if(S4Vectors::isEmpty(grep(rownames(ref_data), pattern = "^ERCC"))){
      stop("Reference data does not contain ERCC spike-in")
    }
    ERCC_index <- grep(rownames(ref_data), pattern = "^ERCC")
    ERCC_counts <- ref_data[ERCC_index, ]
    ref_data <- ref_data[-ERCC_index, ]
    ref_data <- rbind(ref_data, ERCC_counts)
    ERCC_info <- simmethods::ERCC_info
    ERCC_info <- ERCC_info[which(rownames(ERCC_counts) %in% ERCC_info$ERCC_id), ]
    concentration <- ERCC_info$con_Mix1_attomoles_ul
    spikein_abund <- concentration*10^-18*6.022*10^23*other_prior[["volume"]]/other_prior[["dilution.factor"]]
    spikein <- SPARSim::SPARSim_create_spikein_mix(mix_name= "spikein",
                                                   abundance = spikein_abund)
    spikein_set <- SPARSim::SPARSim_create_spikein_set(spikein_mixes = list(spikein = spikein))
    spikein_sample_association <- c(rep("spikein", ncol(ref_data)))
    spikein_abundance <- 0.05
    SPARSim_spikein_parameter <- SPARSim_create_spikein_parameter(spikein_set = spikein_set,
                                                                  spikein_sample = spikein_sample_association,
                                                                  spikein_proportion = spikein_abundance)

  }else{
    SPARSim_spikein_parameter <- NULL
  }

  estimate_formals <- simutils::change_parameters(function_expr = "SPARSim::SPARSim_estimate_parameter_from_data",
                                                  other_prior = other_prior,
                                                  step = "estimation")

  ##############################################################################
  ####                            Estimation                                 ###
  ##############################################################################
  if(verbose){
    cat("Estimating parameters using SPARSim\n")
  }
  # Seed
  set.seed(seed)
  # Estimation
  tryCatch({
    # Estimate parameters from real data and return parameters and detection results
    estimate_detection <- peakRAM::peakRAM(
      estimate_result <- do.call(SPARSim::SPARSim_estimate_parameter_from_data, estimate_formals)
    )
  }, error = function(e){
    as.character(e)
  })
  estimate_result <- list(estimate_result = estimate_result,
                          SPARSim_spikein_parameter = SPARSim_spikein_parameter)
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  estimate_output <- list(estimate_result = estimate_result,
                          estimate_detection = estimate_detection)
  return(estimate_output)
}


#' Simulate Datasets by SPARSim
#'
#' @param parameters A object generated by [SPARSim::SPARSim_simulation()]
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @param return_format A character. Alternative choices: list, SingleCellExperiment,
#' Seurat, h5ad. If you select `h5ad`, you will get a path where the .h5ad file saves to.
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#' @import SPARSim
#' @importFrom stats runif
#' @importFrom BiocGenerics get
#' @export
#' @details
#' In addtion to simulate datasets with default parameters, users want to simulate
#' other kinds of datasets, e.g. a counts matrix with 2 or more cell groups. In
#' SPARSim, you can set extra parameters to simulate datasets.
#'
#' The customed parameters you can set are below:
#' 1. de.prob. You can directly set `other_prior = list(de.prob = 0.2)` to simulate DEGs that account for 20 percent of all genes.
#' 2. fc.group. You can directly set `other_prior = list(fc.group = 3)` to specify the fold change of DEGs.
#' 3. batch.condition. You can input the batch vector that each cell belongs to and set `other_prior = list(batch.condition = xxxx)`. This parameter also determine the number of batches.
#'
#' If users want to simulate groups, they should estimate group parameters by inputting `group.condition` parameter previously. Otherwise, thay can not simulate groups.
#'
#' For more customed parameters in SPARSim, please check [SPARSim::SPARSim_simulation()].
#'
#' @references
#' Baruzzo G, Patuzzi I, Di Camillo B. SPARSim single cell: a count data simulator for scRNA-seq data. Bioinformatics, 2020, 36(5): 1468-1475. <https://doi.org/10.1093/bioinformatics/btz752>
#'
#' Gitlab URL: <https://gitlab.com/sysbiobig/sparsim>
#'
#' @examples
#' ref_data <- simmethods::data
#' ## Estimation
#' group_condition <- as.numeric(simmethods::group_condition)
#' estimate_result <- simmethods::SPARSim_estimation(
#'   ref_data = ref_data,
#'   other_prior = list(group.condition = group_condition),
#'   verbose = TRUE,
#'   seed = 111
#' )
#' ## 1) Simulation (20% proportion of DEGs, fold change 3)
#' simulate_result <- simmethods::SPARSim_simulation(
#'   parameters = estimate_result[["estimate_result"]],
#'   other_prior = list(de.prob = 0.2,
#'                      fc.group = 3),
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
#' ## gene information
#' row_data <- simulate_result[["simulate_result"]][["row_meta"]]
#' table(row_data$de_genes)[2]/4000
#'
#' ## 2) In SPARSim, users can simulate batches when batch.condition parameter is available
#' simulate_result <- simmethods::SPARSim_simulation(
#'   parameters = estimate_result[["estimate_result"]],
#'   other_prior = list(de.prob = 0.2,
#'                      fc.group = 3,
#'                      batch.condition = sample(1:3, 160, replace = TRUE)),
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
#' table(col_data$batch)
#'
#' ## 3) Users can also utilize spike-in genes to simulate datasets. In this case, users
#' ## must input dilution.factor and volume (nanoliter) parameters. Note that the
#' ## reference matrix must contain spike-in gene counts.
#' ref_data <- simmethods::data
#'
#' group_condition <- as.numeric(simmethods::group_condition)
#' estimate_result <- simmethods::SPARSim_estimation(
#'   ref_data = ref_data,
#'   other_prior = list(group.condition = group_condition,
#'                      dilution.factor = 50000,
#'                      volume = 0.01),
#'   verbose = TRUE,
#'   seed = 111
#' )
#' ## check spike-in parameters
#' spikein_params <- estimate_result[["estimate_result"]][["SPARSim_spikein_parameter"]]
#' ## simulate
#' simulate_result <- simmethods::SPARSim_simulation(
#'   parameters = estimate_result[["estimate_result"]],
#'   other_prior = list(de.prob = 0.2,
#'                      fc.group = 3),
#'   return_format = "list",
#'   verbose = TRUE,
#'   seed = 111
#' )
SPARSim_simulation <- function(parameters,
                               other_prior = NULL,
                               return_format,
                               verbose = FALSE,
                               seed
){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("SPARSim", quietly = TRUE)){
    cat("SPARSim is not installed on your device\n")
    cat("Installing SPARSim...\n")
    devtools::install_gitlab("sysbiobig/sparsim")
  }
  other_prior[["dataset_parameter"]] <- parameters[["estimate_result"]]
  other_prior[["spikein_parameter"]] <- parameters[["SPARSim_spikein_parameter"]]
  parameters <- parameters[["estimate_result"]]
  ## Batch
  if(!is.null(other_prior[["batch.condition"]])){
    other_prior[["output_batch_matrix"]] <- TRUE
    if(is.null(other_prior[["distribution"]])){
      other_prior[["distribution"]] <- "normal"
    }
    if(!length(unique(other_prior[["batch.condition"]])) == length(other_prior[["distribution"]])){
      other_prior[["distribution"]] <- rep(other_prior[["distribution"]][1],
                                           length(unique(other_prior[["batch.condition"]])))
    }
    assertthat::assert_that(length(unique(other_prior[["batch.condition"]])) == length(other_prior[["distribution"]]),
                            msg = "The length of batch.condition must equal to the vector of distributions")
    ## Default param A
    if(is.null(other_prior[["param_A"]])){
      for(w in 1:length(other_prior[["distribution"]])){
        dis <- other_prior[["distribution"]][w]
        if(dis == "normal"){
          other_prior[["param_A"]] <- BiocGenerics::append(other_prior[["param_A"]], w-1)
        }
        if(dis == "gamma"){
          other_prior[["param_A"]] <- BiocGenerics::append(other_prior[["param_A"]], w)
        }
      }
    }
    ## Default param B
    if(is.null(other_prior[["param_B"]])){
      for(w in 1:length(other_prior[["distribution"]])){
        dis <- other_prior[["distribution"]][w]
        if(dis == "normal"){
          other_prior[["param_B"]] <- BiocGenerics::append(other_prior[["param_B"]], w)
        }
        if(dis == "gamma"){
          other_prior[["param_B"]] <- BiocGenerics::append(other_prior[["param_B"]], w)
        }
      }
    }
    ## Check param length
    assertthat::assert_that(length(other_prior[["param_A"]]) == length(other_prior[["distribution"]]),
                                   msg = "The length of param_A must equal to the vector of distributions \n")
    assertthat::assert_that(length(other_prior[["param_B"]]) == length(other_prior[["distribution"]]),
                                   msg = "The length of param_B must equal to the vector of distributions \n")


    batch_set_tmp <- purrr::map(.x = seq_len(length(other_prior[["distribution"]])),
                                .f = function(id){
                                  SPARSim::SPARSim_create_batch(name = paste0("Batch", id),
                                                                distribution = other_prior[["distribution"]][id],
                                                                param_A = other_prior[["param_A"]][id],
                                                                param_B = other_prior[["param_B"]][id])
                                }) %>% stats::setNames(paste0("Batch_", seq_len(length(other_prior[["distribution"]]))))
    batch_set <- SPARSim::SPARSim_create_batch_set(batch_list = batch_set_tmp)
    batch_sample_association <- paste0("Batch", other_prior[["batch.condition"]])
    SPARSim_batch_parameter <-SPARSim::SPARSim_create_batch_parameter(batch_set = batch_set,
                                                                      batch_sample = batch_sample_association)
    other_prior[["batch_parameter"]] <- SPARSim_batch_parameter
  }
  ## DEGs
  if(length(parameters) != 1){
    ## gene number
    n_genes <- length(parameters[[1]][['intensity']])
    ## fold change
    if(is.null(other_prior[["fc.group"]])){
      cat("You do not point the fold change of DEGs between two groups, we will set it as 2 \n")
      other_prior[["fc.group"]] <- 2
    }
    ## proportion of DEGs
    if(is.null(other_prior[["de.prob"]])){
      cat("You do not point the percent of DEGs, we will set it as 0.1 \n")
      other_prior[["de.prob"]] <- 0.1
    }

    DE_gene_number <- round(n_genes * other_prior[["de.prob"]])
    DE_group <- simutils::proportionate(number = DE_gene_number,
                                        result_sum_strict = DE_gene_number,
                                        prop = rep(1/(length(parameters)-1), (length(parameters)-1)),
                                        prop_sum_strict = 1,
                                        digits = 0)

    cond_A_param <- parameters[[1]]
    SPARSim_param_with_DE <- list(cond_A = cond_A_param)
    for(group in seq_len(length(parameters)-1)){

      set.seed(seed)
      DE_multiplier <- c(stats::runif(n = ceiling(DE_group[group]/2),
                                      min = 1/other_prior[["fc.group"]],
                                      max = 1/other_prior[["fc.group"]]),
                         stats::runif(n = floor(DE_group[group]/2),
                                      min = other_prior[["fc.group"]],
                                      max = other_prior[["fc.group"]]))
      if(group >= 2){
        pre_DE <- sum(DE_group[1:(group-1)])
        if(group < 3){
          fold_change_multiplier_record <- fold_change_multiplier
        }
        fold_change_multiplier_record[(pre_DE+1):(pre_DE+length(DE_multiplier))] <- DE_multiplier
        fold_change_multiplier <- c(rep(1, pre_DE),
                                    DE_multiplier,
                                    rep(1, n_genes-pre_DE-length(DE_multiplier)))
      }else{
        not_DE_multiplier <- rep(1, n_genes-DE_group[group])
        fold_change_multiplier <- c(DE_multiplier, not_DE_multiplier)
        fold_change_multiplier_record <- NULL
      }
      cond_name <- paste0("cond_", LETTERS[group+1], "_param")
      assign(cond_name, parameters[[group+1]])

      assign(cond_name, SPARSim::SPARSim_create_DE_genes_parameter(
        sim_param = cond_A_param,
        fc_multiplier = fold_change_multiplier,
        N_cells = length(BiocGenerics::get(cond_name)[["lib_size"]]),
        lib_size_DE = BiocGenerics::get(cond_name)[["lib_size"]],
        condition_name = paste0("cond_", LETTERS[group+1]))
      )
      SPARSim_param_with_DE[[paste0("cond_", LETTERS[group+1])]] <- BiocGenerics::get(cond_name)
    }
    other_prior[["dataset_parameter"]] <- SPARSim_param_with_DE
  }

  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################

  simulate_formals <- simutils::change_parameters(function_expr = "SPARSim::SPARSim_simulation",
                                                  other_prior = other_prior,
                                                  step = "simulation")

  condition_num <- length(parameters)
  cell_num <- c()
  for(i in 1:condition_num){
    tmp <- length(parameters[[i]][["lib_size"]])
    cell_num <- BiocGenerics::append(cell_num, tmp)
  }
  cell_num <- sum(cell_num)

  ## Batches
  if(is.null(simulate_formals[["batch_parameter"]])){
    nBatches <- 0
  }else{
    nBatches <- length(simulate_formals[["batch_parameter"]][["batch_set"]])
  }

  ## Group
  if(!is.null(other_prior[["fc.group"]])){
    fc.group <- other_prior[["fc.group"]]
  }else{
    fc.group <- 0
  }
  if(!is.null(other_prior[["de.prob"]])){
    de.prob <- other_prior[["de.prob"]]
  }else{
    de.prob <- 0
  }
  groups <- length(simulate_formals[["dataset_parameter"]])
  # Return to users
  cat(glue::glue("nCells: {cell_num}"), "\n")
  cat(glue::glue("nGenes: {length(parameters[[1]][['intensity']])}"), "\n")
  cat(glue::glue("nGroups: {groups}"), "\n")
  cat(glue::glue("fc.group: {fc.group}"), "\n")
  cat(glue::glue("de.prob: {de.prob}"), "\n")
  cat(glue::glue("nBatches: {nBatches}"), "\n")

  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    cat("Simulating datasets using SPARSim\n")
  }
  # Seed
  set.seed(seed)
  # Estimation
  tryCatch({
    simulate_detection <- peakRAM::peakRAM(
      simulate_result <- do.call(SPARSim::SPARSim_simulation, simulate_formals))
  }, error = function(e){
    as.character(e)
  })
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  ## counts
  counts <- simulate_result[["count_matrix"]]
  if(!is.null(other_prior[["spikein_parameter"]])){
    spikein_index <- grep(rownames(counts), pattern = "^spikein")
    filter_index <- (nrow(counts)-2*length(spikein_index)+1):(nrow(counts)-length(spikein_index))
    counts <- counts[-filter_index, ]
  }

  if(length(parameters) != 1){
    group_name_tmp <- stringr::str_split(colnames(counts),
                                         pattern = "_",
                                         simplify = TRUE)[, 2]
    group_name <- unique(group_name_tmp)
    group <- rep("Group1", ncol(counts))
    for(i in 1:length(group_name)){
      index <- which(group_name[i] == group_name_tmp)
      group[index] <- rep(paste0("Group", i), length(index))
    }
  }

  colnames(counts) <- paste0("Cell", 1:ncol(counts))
  rownames(counts) <- paste0("Gene", 1:nrow(counts))
  ## col_data
  if(!is.null(simulate_formals[["batch_parameter"]])){
    batch <- paste0("Batch", other_prior[["batch.condition"]])
    col_data <- data.frame("cell_name" = colnames(counts),
                           "group" = group,
                           "batch" = batch)
  }else{
    col_data <- data.frame("cell_name" = colnames(counts))
  }

  ## row_data
  if(length(parameters) != 1){
    if(is.null(fold_change_multiplier_record)){
      de_genes <- ifelse(fold_change_multiplier == 1, "no", "yes")
      fc <- fold_change_multiplier
    }else{
      de_genes <- ifelse(fold_change_multiplier_record == 1, "no", "yes")
      fc <- fold_change_multiplier_record
    }
    row_data <- data.frame("gene_name" = rownames(counts),
                           "de_genes" = de_genes,
                           "fc_gene" = fc,
                           "de_genes_group" = c(rep(paste0("Group1_Group",
                                                           2:(length(DE_group)+1)),
                                                    DE_group),
                                                rep("Group1", nrow(counts)-DE_gene_number)))
  }else{
    row_data <- data.frame("gene_name" = rownames(counts))
  }

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

