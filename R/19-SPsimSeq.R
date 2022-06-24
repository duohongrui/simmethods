#' Simulate Datasets by SPsimSeq
#'
#' @param ref_data A matrix for one dataset or a list of datasets with their own
#' names. This is usually unused except for some methods, e.g. SCRIP, scDesign,
#' zingeR, SPsimSeq
#' @param return_format A character. Alternatives choices: list, SingleCellExperiment,
#' Seurat, h5ad
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#'
#' @importFrom SPsimSeq SPsimSeq
#' @importFrom assertthat assert_that
#' @importFrom glue glue
#' @importFrom SingleCellExperiment counts colData rowData
#'
#' @export
#'
SPsimSeq_simulation <- function(ref_data,
                                return_format,
                                verbose = FALSE,
                                seed
){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("SPsimSeq", quietly = TRUE)){
    cat("SPsimSeq is not installed on your device\n")
    cat("Installing SPsimSeq\n")
    BiocManager::install("SPsimSeq")
  }
  other_prior[["s.data"]] <- ref_data
  ## Batch
  if(is.null(other_prior[["batch.condition"]])){
    other_prior[["batch"]] <- rep(1, ncol(ref_data))
  }else{
    other_prior[["batch"]] <- other_prior[["batch.condition"]]
    other_prior[["batch.config"]] <- table(other_prior[["batch"]])/length(other_prior[["batch"]])
  }
  ## Group
  if(is.null(other_prior[["group.condition"]])){
    other_prior[["group"]] <- rep(1, ncol(ref_data))
  }else{
    other_prior[["group"]] <- other_prior[["group.condition"]]
    other_prior[["group.config"]] <- table(other_prior[["group"]])/length(other_prior[["group"]])
  }
  if(is.null(other_prior[["prob.group"]])){
    other_prior[["pDE"]] <- 0.1
  }else{
    other_prior[["pDE"]] <- other_prior[["prob.group"]]
  }
  ## gene number
  if(is.null(other_prior[["nGenes"]])){
    other_prior[["n.genes"]] <- nrow(ref_data)
  }else{
    other_prior[["n.genes"]] <- other_prior[["nGenes"]]
  }
  ## cell number
  if(is.null(other_prior[["nCells"]])){
    other_prior[["tot.samples"]] <- ncol(ref_data)
  }else{
    other_prior[["tot.samples"]] <- other_prior[["nCells"]]
  }
  ## fc
  if(is.null(other_prior[["fc.group"]])){
    other_prior[["lfc.thrld"]] <- 0.5
  }else{
    other_prior[["lfc.thrld"]] <- other_prior[["fc.group"]]
  }
  ## verbose
  other_prior[["verbose"]] <- verbose
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################

  simulate_formals <- simutils::change_parameters(function_expr = "SPsimSeq::SPsimSeq",
                                                  other_prior = other_prior,
                                                  step = "simulation")
  # Return to users
  cat(glue::glue("nCells: {simulate_formals[['tot.samples']]}"), "\n")
  cat(glue::glue("nGenes: {simulate_formals[['n.genes']]}"), "\n")
  cat(glue::glue("nGroups: {length(unique(simulate_formals[['group']]))}"), "\n")
  cat(glue::glue("prob.group: {simulate_formals[['pDE']]}"), "\n")
  cat(glue::glue("nBatches: {length(unique(simulate_formals[['batch']]))}"), "\n")

  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    cat("Simulating datasets using SPsimSeq\n")
  }
  # Seed
  set.seed(seed)
  # Estimation
  tryCatch({
    simulate_detection <- peakRAM::peakRAM(
      simulate_result <- BiocGenerics::do.call(SPsimSeq::SPsimSeq, simulate_formals))
  }, error = function(e){
    as.character(e)
  })
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  ## counts
  colnames(simulate_result[[1]]) <- paste0("Cell", 1:ncol(simulate_result[[1]]))
  rownames(simulate_result[[1]]) <- paste0("Gene", 1:nrow(simulate_result[[1]]))
  ## col_data
  simulate_result[[1]]$Batch <- paste0("Batch", simulate_result[[1]]$Batch)
  simulate_result[[1]]$Group <- paste0("Group", simulate_result[[1]]$Group)
  ## Data format conversion
  simulate_result <- simutils::data_conversion(SCE_object = simulate_result,
                                               return_format = return_format)

  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  simulate_output <- list(simulate_result = simulate_result,
                          simulate_detection = simulate_detection)
  return(simulate_output)
}

