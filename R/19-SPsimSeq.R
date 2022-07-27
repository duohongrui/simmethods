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
#' @details
#' In addtion to simulate datasets with default parameters, users want to simulate
#' other kinds of datasets, e.g. a counts matrix with 2 or more cell groups. In
#' SPsimSeq, you can set extra parameters to simulate datasets.
#'
#' The customed parameters you can set are below:
#' 1. nCells. In SPsimSeq, you can set nCells directly. For example, if you want to simulate 1000 cells, you can type `other_prior = list(nCells = 1000)`.
#' 2. nGenes. You can directly set `other_prior = list(nGenes = 5000)` to simulate 5000 genes.
#' 3. group.condition. You can input cell group information as an integer vector to specify which group that each cell belongs to. See `Examples`.
#' 4. de.prob. You can directly set `other_prior = list(de.prob = 0.2)` to simulate DEGs that account for 20 percent of all genes.
#' 5. fc.group. You can directly set `other_prior = list(fc.group = 2)` to specify the minimum fold change of DEGs.
#' 6. batch.condition. You can input cell batch information as an integer vector to specify which batch that each cell belongs to. See `Examples`.
#'
#' For more customed parameters in SPsimSeq, please check [SPsimSeq::SPsimSeq()].
#'
#' @references
#' Assefa A T, Vandesompele J, Thas O. SPsimSeq: semi-parametric simulation of bulk and single-cell RNA-sequencing data. Bioinformatics, 2020, 36(10): 3276-3278. <https://doi.org/10.1093/bioinformatics/btaa105>
#'
#' Bioconductor URL: <https://www.bioconductor.org/packages/release/bioc/html/SPsimSeq.html>
#'
#' Github URL: <https://github.com/CenterForStatistics-UGent/SPsimSeq>
#'
#' @examples
#' # SPsimSeq can simulate datasets directly without estimation step.
#' ref_data <- simmethods::data
#'
#' # 1) Simulate with default parameters
#' simulate_result <- simmethods::SPsimSeq_simulation(
#'   ref_data = ref_data,
#'   other_prior = NULL,
#'   return_format = "list",
#'   verbose = TRUE,
#'   seed = 111
#' )
#'
#' ## counts
#' counts <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(counts)
#' ## cell information
#' col_data <- simulate_result[["simulate_result"]][["col_meta"]]
#' table(col_data$group)
#' table(col_data$batch)
#'
#'
#' # 2) Simulate two groups (20% proportion of DEGs, minimum 2 fold change)
#' group_condition <- as.numeric(simmethods::group_condition)
#' simulate_result <- simmethods::SPsimSeq_simulation(
#'   ref_data = ref_data,
#'   other_prior = list(nCells = 1000,
#'                      nGenes = 2000,
#'                      group.condition = group_condition,
#'                      de.prob = 0.2,
#'                      fc.group = 2),
#'   return_format = "list",
#'   verbose = TRUE,
#'   seed = 111
#' )
#'
#' ## counts
#' counts <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(counts)
#' ## cell information
#' col_data <- simulate_result[["simulate_result"]][["col_meta"]]
#' table(col_data$group)
#' table(col_data$batch)
#' ## gene information
#' row_data <- simulate_result[["simulate_result"]][["row_meta"]]
#' table(row_data$de_gene)[2]/2000
#'
#'
#' # 3) Simulate two batches
#' group_condition <- as.numeric(simmethods::group_condition)
#' simulate_result <- simmethods::SPsimSeq_simulation(
#'   ref_data = ref_data,
#'   other_prior = list(nCells = 1000,
#'                      nGenes = 2000,
#'                      group.condition = group_condition,
#'                      de.prob = 0.2,
#'                      fc.group = 2,
#'                      batch.condition = sample(1:2, ncol(ref_data), replace = TRUE)),
#'   return_format = "list",
#'   verbose = TRUE,
#'   seed = 111
#' )
#'
#' ## counts
#' counts <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(counts)
#' ## cell information
#' col_data <- simulate_result[["simulate_result"]][["col_meta"]]
#' table(col_data$group)
#' table(col_data$batch)
#' ## gene information
#' row_data <- simulate_result[["simulate_result"]][["row_meta"]]
#' table(row_data$de_gene)[2]/2000
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

