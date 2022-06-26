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
#' DEGs and other variables are usually customed, so before simulating a dataset
#' you must point it out.
#' @param seed An integer of a random seed.
#' @importFrom peakRAM peakRAM
#' @importFrom powsimR estimateParam estimateSpike
#'
#' @return A list contains the estimated parameters and the results of execution
#' detection.
#' @export
#'
powsimR_estimation <- function(ref_data, verbose = FALSE, other_prior, seed){

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



  method_formals <- as.list(formals(powsimR::estimateParam))
  for(param in names(method_formals)){
    names_wait_check <- names(other_prior)
    if(param %in% names_wait_check){
      method_formals[[param]] <- other_prior[[param]]
    }
  }
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
      estimate_result <- do.call(powsimR::estimateParam, method_formals)
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
#' @param parameters A object generated by \code{\link[powsimR]{estimateParam}} or
#' \code{\link[powsimR]{estimateSpike}}
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs and other variables are usually customed, so before simulating a dataset
#' you must point it out.
#' @param return_format A character. Alternatives choices: list, SingleCellExperiment,
#' Seurat, h5ad
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#'
#' @importFrom powsimR Setup simulateDE
#' @importFrom assertthat assert_that
#' @importFrom glue glue
#' @importFrom SingleCellExperiment counts colData rowData
#' @importFrom Seurat as.Seurat
#' @importFrom SeuratDisk SaveH5Seurat Convert
#' @importFrom stringr str_replace
#'
#' @export
#'
powsimR_simulation <- function(parameters,
                               other_prior,
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
  if(!is.null(other_prior[["prob.batch"]]) | !is.null(other_prior[["fc.batch"]])){
    if(is.null(other_prior[["nCells"]])){
      other_prior[["nCells"]] <- parameters[["totalS"]]
    }
    if(length(other_prior[["nCells"]]) == 2){
      other_prior[["n1"]] <- c(ceiling(other_prior[["nCells"]][1]/2), floor(other_prior[["nCells"]][1]/2))
      other_prior[["n2"]] <- c(ceiling(other_prior[["nCells"]][2]/2), floor(other_prior[["nCells"]][2]/2))
    }else{
      other_prior[["n1"]] <- c(ceiling(ceiling(other_prior[["nCells"]]/2)/2),
                               floor(ceiling(other_prior[["nCells"]]/2)/2))
      other_prior[["n2"]] <- c(ceiling(floor(other_prior[["nCells"]]/2)/2),
                               floor(floor(other_prior[["nCells"]]/2)/2))
    }
  }else{
    if(is.null(other_prior[["nCells"]])){
      other_prior[["nCells"]] <- parameters[["totalS"]]
    }
    if(length(other_prior[["nCells"]]) != 1){
      other_prior[["n1"]] <- other_prior[["nCells"]][1]
      other_prior[["n2"]] <- other_prior[["nCells"]][2]
    }else{
      other_prior[["n1"]] <- ceiling(other_prior[["nCells"]]/2)
      other_prior[["n2"]] <- floor(other_prior[["nCells"]]/2)
    }
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
    other_prior[["pLFC"]] <- other_prior[["fc.group"]]
  }
  ## fold change of genes in batches
  if(is.null(other_prior[["fc.batch"]])){
    other_prior[["bLFC"]] <- NULL
  }else{
    other_prior[["bLFC"]] <- other_prior[["fc.batch"]]
  }
  ## proportion of genes in batches
  if(is.null(other_prior[["bLFC"]])){
    other_prior[["p.B"]] <- NULL
  }else{
    other_prior[["p.B"]] <- 0.1
  }


  cat(glue::glue("nCells: {other_prior[['nCells']]}"), "\n")
  cat(glue::glue("nGenes: {other_prior[['ngenes']]}"), "\n")
  cat(glue::glue("nGroups: 2"), "\n")
  cat(glue::glue("de.prob: {other_prior[['p.DE']]}"), "\n")
  cat(glue::glue("fc.group: {other_prior[['pLFC']]}"), "\n")
  cat(glue::glue("nBatches: 2 \n"), "\n")
  cat(glue::glue("fc.batch: {other_prior[['bLFC']]}"), "\n")


  setup_formals <- as.list(formals(powsimR::Setup))
  for(param in names(setup_formals)){
    names_wait_check <- names(other_prior)
    if(param %in% names_wait_check){
      setup_formals[[param]] <- other_prior[[param]]
    }
  }

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

  simulate_formals <- as.list(formals(powsimR::simulateDE))
  for(param in names(simulate_formals)){
    names_wait_check <- names(other_prior)
    if(param %in% names_wait_check){
      simulate_formals[[param]] <- other_prior[[param]]
    }
  }

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
  if(length(simulate_result[["Counts"]]) == 2){

    counts_1 <- simulate_result[["Counts"]][[1]][[1]]
    counts_2 <- simulate_result[["Counts"]][[2]][[1]]
    counts <- cbind(counts_1, counts_2)

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
                           "de_fc" = simulate_result[["DESetup"]][["pLFC"]][[1]],
                           "batch_genes" = ifelse(simulate_result[["DESetup"]][["bLFC"]][[1]] == 0, FALSE, TRUE),
                           "batch_fc" = simulate_result[["DESetup"]][["bLFC"]][[1]])
  }else{
    row_data <- data.frame("gene_name" = rownames(counts),
                           "de_genes" = ifelse(simulate_result[["DESetup"]][["pLFC"]][[1]] == 0, FALSE, TRUE),
                           "de_fc" = simulate_result[["DESetup"]][["pLFC"]][[1]])
  }
  if(length(simulate_result[["Counts"]]) == 2){
    col_data <- data.frame("cell_name" = colnames(counts),
                           "group" = paste0("Group",
                                            c(rep(1, other_prior[["n1"]][1]),
                                              rep(2, other_prior[["n1"]][2]),
                                              rep(1, other_prior[["n2"]][1]),
                                              rep(2, other_prior[["n2"]][2]))),
                           "batch" = paste0("Batch",
                                            c(rep(1, ncol(counts_1)), rep(2, ncol(counts_2)))))
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

