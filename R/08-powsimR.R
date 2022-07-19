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
#' @return A list contains the estimated parameters and the results of execution
#' detection.
#' @export
#' @references
#' Vieth B, Ziegenhain C, Parekh S, et al. powsimR: power analysis for bulk and single cell RNA-seq experiments[J]. Bioinformatics, 2017, 33(21): 3486-3488. <https://doi.org/10.1093/bioinformatics/btx435>
#'
#' Github URL: <https://github.com/bvieth/powsimR>
#'
powsimR_estimation <- function(ref_data,
                               verbose = FALSE,
                               other_prior,
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
  # Check
  if(is.null(other_prior[["RNAseq"]])){
    other_prior[["RNAseq"]] <- "singlecell"
    # stop("Please input the parameter of RNAseq: one of 'bulk' and 'singlecell'")
  }
  if(is.null(other_prior[["Protocol"]])){
    other_prior[["Protocol"]] <- "UMI"
    # stop("Please input the parameter of Protocol: one of 'UMI' and 'Read'")
  }
  # Sub function
  if(is.null(other_prior[["Distribution"]])){
    other_prior[["Distribution"]] <- "NB"
    # stop("Please input the parameter of Distribution: one of 'NB' and 'ZINB'")
  }
  if(is.null(other_prior[["Normalisation"]])){
    other_prior[["Normalisation"]] <- "TMM"
    # stop("Please input the parameter of Normalisation: one of 'TMM', 'MR', 'PosCounts', 'UQ', 'scran', 'Linnorm', 'sctransform', 'SCnorm', 'Census', 'depth' and 'none'")
  }

  if(is.null(other_prior[["subfunction"]])){
    other_prior[["subfunction"]] <- "estimateParam"
  }

  estimate_formals <- simutils::change_parameters(function_expr = "powsimR::estimateParam",
                                                  other_prior = other_prior,
                                                  step = "estimation")
  ## Spike-in
  if(other_prior[["subfunction"]] == "estimateSpike"){
    if(is.null(other_prior[["spikeData"]])){
      stop("Please input the parameter of spikeData")
    }
    if(is.null(other_prior[["spikeInfo"]])){
      stop("Please input the parameter of spikeInfo")
    }
    if(is.null(other_prior[["Protocol"]])){
      stop("Please input the parameter of Protocol: one of 'UMI' and 'Read'")
    }
    sub_method_formals <- as.list(formals(powsimR::estimateSpike))
    for(param in names(sub_method_formals)){
      names_wait_check <- names(other_prior)
      if(param %in% names_wait_check){
        sub_method_formals[[param]] <- other_prior[[param]]
      }
    }
  }

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
    estSpikeRes <- NULL
    if(other_prior[["subfunction"]] == "estimateSpike"){
      cat("Estimating parameters using estimateSpike function\n")
      estSpikeRes <- do.call(powsimR::estimateSpike, sub_method_formals)
    }
  }, error = function(e){
    as.character(e)
  })
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  estimate_output <- list(estimate_result = estimate_result,
                          estimate_detection = estimate_detection,
                          estSpikeRes = estSpikeRes)
  return(estimate_output)
}


#' Simulate Datasets by powsimR
#'
#' @param parameters A object generated by [powsimR::estimateParam()] or [powsimR::estimateSpike()].
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs and other variables are usually customed, so before simulating a dataset
#' you must point it out.
#' @param return_format A character. Alternatives choices: list, SingleCellExperiment,
#' Seurat, h5ad
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#' @importFrom powsimR Setup simulateDE
#' @importFrom simutils proportionate
#'
#' @export
#'
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

  other_prior[["estParamRes"]] <- parameters
  other_prior[["setup.seed"]] <- seed
  other_prior[["verbose"]] <- verbose

  if(!is.null(parameters[["estSpikeRes"]])){
    other_prior[["estSpikeRes"]] <- parameters[["estSpikeRes"]]
  }
  ##############################################################################
  ####                                Setup                                  ###
  ##############################################################################

  ## Sim num
  other_prior[["nsims"]] <- 1

  ## genes
  if(is.null(other_prior[["nGenes"]])){
    other_prior[["ngenes"]] <- parameters[["totalG"]]
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
  }else{
    prob.group <- c(0.5, 0.5)
  }
  ### Assign cells
  if(is.null(other_prior[["nCells"]])){
    nCells <- parameters[["totalS"]]
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
    # stop("Please input the parameter of Normalisation: one of 'TMM', 'MR', 'PosCounts', 'UQ', 'scran', 'Linnorm', 'sctransform', 'SCnorm', 'Census', and 'depth'")
  }
  if(is.null(other_prior[['DEmethod']])){
    other_prior[['DEmethod']] <- "limma-trend"
    # stop("Please input the parameter of DEmethod: one of 'T-Test', 'edgeR-LRT', 'edgeR-QL', 'edgeR-zingeR', 'edgeR-ZINB-WaVE', 'limma-voom', 'limma-trend', 'DESeq2', 'DESeq2-zingeR', 'DESeq2-ZINB-WaVE', 'ROTS', 'baySeq', 'NOISeq', 'EBSeq', 'MAST', 'BPSC', 'scDD', 'DECENT'")
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
                           "group" = paste0("Group",
                                            c(rep(1, other_prior[["n1"]][1]),
                                              rep(2, other_prior[["n1"]][2]),
                                              rep(1, other_prior[["n2"]][1]),
                                              rep(2, other_prior[["n2"]][2]))),
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

