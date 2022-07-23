#' Estimate Parameters From Real Datasets by powsimR
#'
#' This function is used to estimate useful parameters from a real dataset by
#' using \code{estimateParam} or \code{estimateSpike} function in powsimR package.
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
#' @importFrom powsimR estimateParam estimateSpike
#' @details
#' powsimR provides some choices for users to select suitable parameters according
#' to different types of data, platforms, normalization methods, distributions and
#' so on.
#' 1. RNAseq. "bulk" or "singlecell" (default).
#' 2. Protocol. Options are "UMI" (default) (e.g. 10X Genomics, CEL-seq2) or "Read" (e.g. Smart-seq2).
#' 3. Distribution. "NB" (default) for negative binomial or "ZINB" for zero-inflated negative binomial distribution fitting.
#' 4. Normalisation. "TMM" (default), "MR", "PosCounts", "UQ", "scran", "Linnorm", "SCnorm", "Census", "depth", "none".
#'
#' powsimR also provides an another choice to estimate parameters (not neccessary)
#' via spike-ins. If users want to use this, make sure that the reference data
#' must contain ERCC spike-in counts. In addtion, users must set dilution.factor and
#' volume information by `other_prior = list(dilution.factor = xxx, volume = xxx)`.
#' For more instructions, see `Examples`.
#' @return A list contains the estimated parameters and the results of execution
#' detection.
#' @export
#' @references
#' Vieth B, Ziegenhain C, Parekh S, et al. powsimR: power analysis for bulk and single cell RNA-seq experiments. Bioinformatics, 2017, 33(21): 3486-3488. <https://doi.org/10.1093/bioinformatics/btx435>
#'
#' Github URL: <https://github.com/bvieth/powsimR>
#'
#' @examples
#' ref_data <- simmethods::data
#'
#' ## Estimate parameters without ERCC spike-in
#' estimate_result <- powsimR_estimation(ref_data = ref_data,
#'                                       other_prior = list(RNAseq = "singlecell",
#'                                                          Protocol = "UMI",
#'                                                          Normalisation = "scran"),
#'                                       verbose = TRUE,
#'                                       seed = 111)
#'
#' ## Estimate parameters with ERCC spike-in
#' ### Make sure there are ERCC names in reference data
#' rownames(ref_data)[grep(rownames(ref_data), pattern = "^ERCC")]
#' ### Users must input the dilution.factor and volume (nanoliter) to determine the ERCC molecules
#' estimate_result <- powsimR_estimation(ref_data = ref_data,
#'                                       other_prior = list(RNAseq = "singlecell",
#'                                                          Protocol = "UMI",
#'                                                          Normalisation = "scran",
#'                                                          dilution.factor = 50000,
#'                                                          volume = 1),
#'                                       verbose = TRUE,
#'                                       seed = 111)
powsimR_estimation <- function(ref_data,
                               verbose = FALSE,
                               other_prior = NULL,
                               seed){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("powsimR", quietly = TRUE)){
    cat("powsimR is not installed on your device\n")
    cat("Installing powsimR...\n")
    devtools::install_github("bvieth/powsimR", build_vignettes = FALSE, dependencies = TRUE)
  }
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  if(!is.matrix(ref_data)){
    ref_data <- as.matrix(ref_data)
  }
  # Add ref_data to other_prior
  other_prior[["countData"]] <- ref_data
  other_prior[["verbose"]] <- verbose
  # Check
  if(is.null(other_prior[["RNAseq"]])){
    other_prior[["RNAseq"]] <- "singlecell"
  }
  if(is.null(other_prior[["Protocol"]])){
    other_prior[["Protocol"]] <- "UMI"
  }
  # Sub function
  if(is.null(other_prior[["Distribution"]])){
    other_prior[["Distribution"]] <- "NB"
  }
  if(is.null(other_prior[["Normalisation"]])){
    other_prior[["Normalisation"]] <- "TMM"
  }

  ## Spike-in
  if(!is.null(other_prior[["volume"]]) & !is.null(other_prior[["dilution.factor"]])){
    spikeData <- ref_data[grep(rownames(ref_data), pattern = "^ERCC"), ]
    concentration <- simmethods::ERCC_info$con_Mix1_attomoles_ul
    spikeInfo <- data.frame(SpikeID = simmethods::ERCC_info$ERCC_id,
                            SpikeInput = concentration*10^-18*6.022*10^23*other_prior[["volume"]]/other_prior[["dilution.factor"]],
                            row.names = simmethods::ERCC_info$ERCC_id)
    spikeInfo <- spikeInfo[rownames(spikeData), ]
    other_prior[["spikeData"]] <- spikeData
    other_prior[["spikeInfo"]] <- spikeInfo
    spike_in_formals <- simutils::change_parameters(function_expr = "powsimR::estimateSpike",
                                                    other_prior = other_prior,
                                                    step = "estimation")
    spike_in_formals$Normalisation <- "depth"
  }else spike_in_formals <- NULL
  estimate_formals <- simutils::change_parameters(function_expr = "powsimR::estimateParam",
                                                  other_prior = other_prior,
                                                  step = "estimation")

  ##############################################################################
  ####                            Estimation                                 ###
  ##############################################################################
  if(verbose){
    cat("Estimating parameters using powsimR\n")
  }
  # Seed
  set.seed(seed)
  # Estimation
  tryCatch({
    # Estimate parameters from real data and return parameters and detection results
    cat("Estimating parameters using estimateParam function\n")
    estimate_detection <- peakRAM::peakRAM(
      estimate_result <- do.call(powsimR::estimateParam, estimate_formals)
    )
    if(!is.null(spike_in_formals)){
      estSpikeRes <- do.call(powsimR::estimateSpike, spike_in_formals)
    }else estSpikeRes <- NULL
  }, error = function(e){
    as.character(e)
  })
  estimate_result <- list(estimate_result = estimate_result,
                          estSpikeRes = estSpikeRes)
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  estimate_output <- list(estimate_result = estimate_result,
                          estimate_detection = estimate_detection)
  return(estimate_output)
}


#' Simulate Datasets by powsimR
#'
#' @param parameters A object generated by [powsimR::estimateParam()] or [powsimR::estimateSpike()].
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @param return_format A character. Alternative choices: list, SingleCellExperiment,
#' Seurat, h5ad. If you select `h5ad`, you will get a path where the .h5ad file saves to.
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#' @importFrom powsimR Setup simulateDE
#' @importFrom simutils proportionate
#' @export
#' @details
#' powsimR provides some choices for users to select suitable parameters according
#' to different normalization methods, and methods for differential expressed analysis.
#' 1. Normalisation. "TMM" (default), "MR", "PosCounts", "UQ", "scran", "Linnorm", "sctransform", "SCnorm", "Census", "depth".
#' 2. DEmethod. "T-Test", "edgeR-LRT", "edgeR-QL", "edgeR-zingeR", "edgeR-ZINB-WaVE", "limma-voom", "limma-trend" (default), "DESeq2", "DESeq2-zingeR", "DESeq2-ZINB-WaVE", "ROTS", "baySeq", "NOISeq", "EBSeq", "MAST", "BPSC", "scDD", "DECENT".
#'
#' In addtion to simulate datasets with default parameters, users want to simulate
#' other kinds of datasets, e.g. a counts matrix with 2 or more cell groups. In
#' powsimR, you can set extra parameters to simulate datasets.
#'
#' The customed parameters you can set are below:
#' 1. nCells. In powsimR, You can directly set `other_prior = list(nCells = 5000)` to simulate 5000 cells.
#' 2. nGenes. You can directly set `other_prior = list(nGenes = 5000)` to simulate 5000 genes.
#' 3. de.prob. You can directly set `other_prior = list(de.prob = 0.2)` to simulate DEGs that account for 20 percent of all genes.
#' 4. prob.group. You can directly set `other_prior = list(prob.group = c(0.4, 0.6))` to assign two proportions of cell groups. Note that the the length of the vector must be **2**.
#' 5. fc.group. You can directly set `other_prior = list(fc.group = 2)` to specify the fold change of DEGs.
#' 6. nBatches. You can not directly set `other_prior = list(nBatches = 3)` to simulate 3 batches. Instead, you should set `other_prior = list(prob.batch = c(0.3, 0.4, 0.3))` to reach the goal.
#' 7. fc.batch. You can directly set `other_prior = list(fc.batch = 2)` to specify the fold change of genes between batches.
#'
#' For more customed parameters in powsimR, please check [powsimR::simulateDE()].
#' @references
#' Vieth B, Ziegenhain C, Parekh S, et al. powsimR: power analysis for bulk and single cell RNA-seq experiments. Bioinformatics, 2017, 33(21): 3486-3488. <https://doi.org/10.1093/bioinformatics/btx435>
#'
#' Github URL: <https://github.com/bvieth/powsimR>
#'
#' @examples
#' ref_data <- simmethods::data
#'
#' ## Estimate parameters with ERCC spike-in
#' ### Make sure there are ERCC names in reference data
#' rownames(ref_data)[grep(rownames(ref_data), pattern = "^ERCC")]
#' ### Users must input the dilution.factor and volume (nanoliter) to determine the ERCC molecules
#' estimate_result <- powsimR_estimation(ref_data = ref_data,
#'                                       other_prior = list(RNAseq = "singlecell",
#'                                                          Protocol = "UMI",
#'                                                          Normalisation = "scran",
#'                                                          dilution.factor = 50000,
#'                                                          volume = 1),
#'                                       verbose = TRUE,
#'                                       seed = 111)
#'
#' # (1) Simulate 500 cells and 2000 genes
#' simulate_result <- simmethods::powsimR_simulation(parameters = estimate_result[["estimate_result"]],
#'                                                   other_prior = list(nCells = 500,
#'                                                                      nGenes = 2000),
#'                                                   return_format = "list",
#'                                                   verbose = TRUE,
#'                                                   seed = 111)
#' count_data <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(count_data)
#' ## In powsimR, users always get two groups of cells, and the numbers of cells in
#' ## different groups are the same (default, prob.group = c(0.5, 0.5)).
#' col_data <- simulate_result[["simulate_result"]][["col_meta"]]
#' table(col_data$group)
#' ## row_data
#' row_data <- simulate_result[["simulate_result"]][["row_meta"]]
#' table(row_data$de_genes)[2]/2000
#'
#' # (2) Simulate two groups with 20% of DEGs and 4 fold change.
#' simulate_result <- simmethods::powsimR_simulation(parameters = estimate_result[["estimate_result"]],
#'                                                   other_prior = list(prob.group = c(0.3, 0.7),
#'                                                                      de.prob = 0.2,
#'                                                                      fc.group = 4),
#'                                                   return_format = "list",
#'                                                   verbose = TRUE,
#'                                                   seed = 111)
#' count_data <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(count_data)
#' ## cell information
#' col_data <- simulate_result[["simulate_result"]][["col_meta"]]
#' table(col_data$group)
#' ## gene information
#' row_data <- simulate_result[["simulate_result"]][["row_meta"]]
#' table(row_data$de_genes)[2]/4000
#'
#'
#' # (3) Simulate two batches (2 fold change)
#' simulate_result <- simmethods::powsimR_simulation(parameters = estimate_result[["estimate_result"]],
#'                                                   other_prior = list(prob.batch = c(0.4, 0.6),
#'                                                                      fc.batch = 2),
#'                                                   return_format = "list",
#'                                                   verbose = TRUE,
#'                                                   seed = 111)
#' count_data <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(count_data)
#' ## cell information
#' col_data <- simulate_result[["simulate_result"]][["col_meta"]]
#' table(col_data$group)
#' table(col_data$batch)
#' ## gene information
#' row_data <- simulate_result[["simulate_result"]][["row_meta"]]
#' head(row_data)
#' table(row_data$de_fc)
#' table(row_data$batch_genes)
powsimR_simulation <- function(parameters,
                               other_prior = NULL,
                               return_format,
                               verbose = FALSE,
                               seed
){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("powsimR", quietly = TRUE)){
    cat("powsimR is not installed on your device\n")
    cat("Installing powsimR\n")
    devtools::install_github("bvieth/powsimR", build_vignettes = FALSE, dependencies = TRUE)
  }
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################

  other_prior[["estParamRes"]] <- parameters$estimate_result
  other_prior[["setup.seed"]] <- seed
  other_prior[["verbose"]] <- verbose

  if(!is.null(parameters[["estSpikeRes"]])){
    other_prior[["estSpikeRes"]] <- parameters$estSpikeRes
  }
  ##############################################################################
  ####                                Setup                                  ###
  ##############################################################################

  ## Sim num
  other_prior[["nsims"]] <- 1

  ## genes
  if(is.null(other_prior[["nGenes"]])){
    other_prior[["ngenes"]] <- parameters$estimate_result[["totalG"]]
  }else{
    other_prior[["ngenes"]] <- other_prior[["nGenes"]]
  }
  ## cells
  ### prob.batch
  if(!is.null(other_prior[["prob.batch"]])){
    prob.batch <- other_prior[["prob.batch"]]
  }else{
    prob.batch <- 1
  }
  ### prob.group
  if(!is.null(other_prior[["prob.group"]])){
    prob.group <- other_prior[["prob.group"]]
    if(sum(prob.group) != 1){
      stop("The sum of prob.group is not euqal to 1")
    }
    if(length(prob.group) != 2){
      stop("The length of prob.group is not euqal to 2")
    }
  }else{
    prob.group <- c(0.5, 0.5)
  }
  ### Assign cells
  if(is.null(other_prior[["nCells"]])){
    nCells <- parameters$estimate_result[["totalS"]]
  }else{
    nCells <- other_prior[["nCells"]]
  }

  group1_cell <- round(nCells * prob.group[1])
  group2_cell <- nCells - group1_cell

  if(length(prob.batch) == 1){
    nbatches <- 1
    other_prior[["n1"]] <- group1_cell
    other_prior[["n2"]] <- group2_cell
  }else{
    nbatches <- length(prob.batch)
    other_prior[["n1"]] <- simutils::proportionate(number = group1_cell,
                                                   result_sum_strict = group1_cell,
                                                   prop = prob.batch,
                                                   prop_sum_strict = 1,
                                                   digits = 0)
    other_prior[["n2"]] <- simutils::proportionate(number = group2_cell,
                                                   result_sum_strict = group1_cell,
                                                   prop = prob.batch,
                                                   prop_sum_strict = 1,
                                                   digits = 0)
  }

  ## proportion of genes in groups
  if(is.null(other_prior[["de.prob"]])){
    other_prior[["p.DE"]] <- 0.1
  }else{
    other_prior[["p.DE"]] <- other_prior[["de.prob"]]
  }
  ## fold change of genes in groups
  if(is.null(other_prior[["fc.group"]])){
    other_prior[["pLFC"]] <- 1
  }else{
    other_prior[["pLFC"]] <- log2(other_prior[["fc.group"]])
  }
  ## fold change of genes in batches
  if(is.null(other_prior[["fc.batch"]])){
    if(nbatches == 1){
      other_prior[["bLFC"]] <- NULL
    }else{
      other_prior[["bLFC"]] <- 1
    }
  }else{
    other_prior[["bLFC"]] <- log2(other_prior[["fc.batch"]])
  }
  ## proportion of genes in batches
  if(nbatches == 1){
    other_prior[["p.B"]] <- NULL
  }else{
    if(is.null(other_prior[["p.B"]])){
      other_prior[["p.B"]] <- 0.1
    }
  }


  cat(glue::glue("nCells: {nCells}"), "\n")
  cat(glue::glue("nGenes: {other_prior[['ngenes']]}"), "\n")
  cat(glue::glue("nGroups: 2"), "\n")
  cat(glue::glue("de.prob: {other_prior[['p.DE']]}"), "\n")
  cat(glue::glue("fc.group: {2^other_prior[['pLFC']]}"), "\n")
  cat(glue::glue("nBatches: {nbatches}"), "\n")
  cat(glue::glue("fc.batch: {2^other_prior[['bLFC']]}"), "\n")

  setup_formals <- simutils::change_parameters(function_expr = "powsimR::Setup",
                                               other_prior = other_prior,
                                               step = "simulation")

  SetupRes <- do.call(powsimR::Setup, setup_formals)

  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  other_prior[["SetupRes"]] <- SetupRes
  other_prior[["Counts"]] <- TRUE

  if(is.null(other_prior[["Normalisation"]])){
    other_prior[["Normalisation"]] <- "TMM"
  }
  if(is.null(other_prior[['DEmethod']])){
    other_prior[['DEmethod']] <- "limma-trend"
  }
  simulate_formals <- simutils::change_parameters(function_expr = "powsimR::simulateDE",
                                                  other_prior = other_prior,
                                                  step = "simulation")

  if(verbose){
    cat("Simulating datasets using powsimR\n")
  }
  # Estimation
  tryCatch({
    simulate_detection <- peakRAM::peakRAM(
      simulate_result <- do.call(powsimR::simulateDE, simulate_formals)
      )
  }, error = function(e){
    as.character(e)
  })
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  # Extract count data
  if(length(simulate_result[["Counts"]]) != 1){
    counts <- c()
    for(name in names(simulate_result[["Counts"]])){
      counts <- cbind(counts, simulate_result[["Counts"]][[name]][[1]])
    }
  }else{
    counts <- simulate_result[["Counts"]][[1]][[1]]
  }
  # Save the old name
  cell_tmp_name <- colnames(counts)
  gene_tmp_name <- rownames(counts)
  # Rename
  colnames(counts) <- paste0("Cell", c(1:ncol(counts)))
  rownames(counts) <- paste0("Gene", c(1:nrow(counts)))
  # Extract infomation
  if(!is.null(SetupRes[["DESetup"]][["bLFC"]])){
    row_data <- data.frame("gene_name" = rownames(counts),
                           "de_genes" = ifelse(simulate_result[["DESetup"]][["pLFC"]][[1]] == 0, FALSE, TRUE),
                           "de_fc" = ifelse(simulate_result[["DESetup"]][["pLFC"]][[1]] == 0, 0,
                                            2^simulate_result[["DESetup"]][["pLFC"]][[1]]),
                           "batch_genes" = ifelse(simulate_result[["DESetup"]][["bLFC"]][[1]] == 0, FALSE, TRUE),
                           "batch_fc" = ifelse(simulate_result[["DESetup"]][["bLFC"]][[1]] ==0, 0,
                                               2^simulate_result[["DESetup"]][["bLFC"]][[1]]))
  }else{
    row_data <- data.frame("gene_name" = rownames(counts),
                           "de_genes" = ifelse(simulate_result[["DESetup"]][["pLFC"]][[1]] == 0, FALSE, TRUE),
                           "de_fc" = ifelse(simulate_result[["DESetup"]][["pLFC"]][[1]] == 0, 0,
                                            2^simulate_result[["DESetup"]][["pLFC"]][[1]]))
  }
  if(length(simulate_result[["Counts"]]) != 1){
    col_data <- data.frame("cell_name" = colnames(counts),
                           "group" = c(rep("Group1", other_prior[["n1"]][1]),
                                       rep("Group2", other_prior[["n2"]][1]),
                                       rep("Group1", other_prior[["n1"]][2]),
                                       rep("Group2", other_prior[["n2"]][2])),
                           "batch" = rep(paste0("Batch", 1:nbatches), (other_prior[["n1"]] + other_prior[["n2"]])))
  }else{
    col_data <- data.frame("cell_name" = colnames(counts),
                           "group" = paste0("Group",
                                            stringr::str_split(cell_tmp_name, pattern = "_n", simplify = T)[, 2]))
  }
  # SingleCellExperiment object
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

