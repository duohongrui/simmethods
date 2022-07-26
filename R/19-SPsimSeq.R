#' Simulate Datasets by SPsimSeq
#'
#' @param ref_data A matrix for one dataset or a list of datasets with their own
#' names. This is usually unused except for some methods, e.g. SCRIP, scDesign,
#' zingeR, SPsimSeq.
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @param return_format A character. Alternative choices: list, SingleCellExperiment,
#' Seurat, h5ad. If you select `h5ad`, you will get a path where the .h5ad file saves to.
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#' @export
#'
#' @references
#' Assefa A T, Vandesompele J, Thas O. SPsimSeq: semi-parametric simulation of bulk and single-cell RNA-sequencing data[J]. Bioinformatics, 2020, 36(10): 3276-3278. <https://doi.org/10.1093/bioinformatics/btaa105>
#'
#' Bioconductor URL: <https://www.bioconductor.org/packages/release/bioc/html/SPsimSeq.html>
#'
#' Github URL: <https://github.com/CenterForStatistics-UGent/SPsimSeq>
#'
SPsimSeq_simulation <- function(ref_data,
                                other_prior = NULL,
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
  if(is.null(other_prior[["de.prob"]])){
    other_prior[["pDE"]] <- 0.1
  }else{
    other_prior[["pDE"]] <- other_prior[["de.prob"]]
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
    other_prior[["lfc.thrld"]] <- 2
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
  cat(glue::glue("de.prob: {simulate_formals[['pDE']]}"), "\n")
  cat(glue::glue("fc.group: {other_prior[['lfc.thrld']]}"), "\n")
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
  counts <- counts(simulate_result[[1]])
  colnames(counts) <- paste0("Cell", 1:ncol(counts))
  rownames(counts) <- paste0("Gene", 1:nrow(counts))
  ## cell information
  col_data <- data.frame("cell_name" = colnames(counts),
                         "group" = paste0("Group", simulate_result[[1]]$Group),
                         "batch" = paste0("Batch", simulate_result[[1]]$Batch))
  rownames(col_data) <- col_data$cell_name
  ## gene information
  row_data <- as.data.frame(SummarizedExperiment::rowData(simulate_result[[1]]))
  if(is.null(other_prior[["group.condition"]])){
    row_data <- data.frame("gene_name" = rownames(counts))
  }else{
    row_data <- data.frame("gene_name" = rownames(counts),
                           "de_gene" = ifelse(row_data$DE.ind, "yes", "no"))
  }
  rownames(row_data) <- row_data$gene_name
  # Establish SingleCellExperiment
  simulate_result <- SingleCellExperiment::SingleCellExperiment(list(counts = counts),
                                                                colData = col_data,
                                                                rowData = row_data)
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

