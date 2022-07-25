#' Estimate Parameters From Real Datasets by BEARscc
#'
#' This function is used to estimate useful parameters from a real dataset by
#' using \code{estimate_noiseparameters} function in BEARscc package.
#'
#' @param ref_data A count matrix. Each row represents a gene and each column
#' represents a cell.
#' @param verbose Logical.
#' @param seed An integer of a random seed.
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. When
#' you use BEARscc, you must input the \code{dilution.factor} and the \code{volume}
#' information to calculate the number of molecules of spike-ins.
#' @import BEARscc
#' @importFrom assertthat not_empty
#' @importFrom S4Vectors metadata
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment altExp
#' @importFrom methods as
#' @importFrom BiocGenerics do.call
#' @importFrom stats na.omit
#' @return A list contains the estimated parameters and the results of execution
#' detection.
#' @export
#'
#' @references
#' Severson D T, Owen R P, White M J, et al. BEARscc determines robustness of single-cell clusters using simulated technical replicates. Nature communications, 2018, 9(1): 1-7. <https://doi.org/10.1038/s41467-018-03608-y>
#'
#' Bioconductor URL: <https://www.bioconductor.org/packages/release/bioc/html/BEARscc.html>
#'
#' Github URL: <https://github.com/seversond12/BEARscc>
BEARscc_estimation <- function(ref_data,
                               verbose = FALSE,
                               other_prior = NULL,
                               seed){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("BEARscc", quietly = TRUE)){
    cat("BEARscc is not installed on your device\n")
    cat("Installing BEARscc...\n")
    BiocManager::install("BEARscc")
  }
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  if(!is.matrix(ref_data)){
    ref_data <- as.matrix(ref_data)
  }
  rownames(ref_data) <- stringr::str_replace_all(rownames(ref_data),
                                                 pattern = '[_]|[;]|[.]|[__]|[/]|[,]',
                                                 replacement = '--')
  colnames(ref_data) <- paste0("Cell", 1:ncol(ref_data))
  ## Calculate the number of molecules of spike-ins
  if(is.null(other_prior[["dilution.factor"]])){
    stop("Please input dilution.factor.")
  }
  if(is.null(other_prior[["volume"]])){
    stop("Please input volume")
  }
  if(!assertthat::not_empty(grep(pattern = '^ERCC-', rownames(ref_data)))){
    stop("ref_data must contain the counts of ERCC spike-in and the names of the molecule must be formal.")
  }
  ## ERCC
  ERCC_count <- ref_data[grep(pattern = '^ERCC-', rownames(ref_data)), ]
  concentration <- simmethods::ERCC_info$con_Mix1_attomoles_ul
  ERCC_meta <- data.frame(Transcripts = concentration*10^-18*6.022*10^23*other_prior[["volume"]]/other_prior[["dilution.factor"]],
                          row.names = simmethods::ERCC_info$ERCC_id)
  ERCC_meta <- as.data.frame(ERCC_meta[match(rownames(ERCC_count), rownames(ERCC_meta)), ])
  rownames(ERCC_meta) <- rownames(ERCC_count)
  colnames(ERCC_meta) <- "Transcripts"
  ref_data <- ref_data[-grep(pattern = '^ERCC-', rownames(ref_data)), ]
  ## Gene name transformation
  ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

  id_convert <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                              "transcript_biotype",
                                              "external_gene_name"), mart = ensembl) %>%
    dplyr::filter("external_gene_name" != "")
  gene_filter <- id_convert$external_gene_name[stats::na.omit(match(rownames(ref_data),
                                                                    id_convert$external_gene_name))]
  ref_data <- ref_data[gene_filter, ]
  rownames(ref_data) <- id_convert$ensembl_gene_id[stats::na.omit(match(rownames(ref_data),
                                                                  id_convert$external_gene_name))]
  ref_data <- rbind(ref_data, ERCC_count)
  ## Data format
  ref_data <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = as.matrix(ref_data)))
  ref_data <- methods::as(ref_data, "SingleCellExperiment")
  S4Vectors::metadata(ref_data) <- list(spikeConcentrations = ERCC_meta)
  SummarizedExperiment::assay(ref_data, "observed_expression") <- SingleCellExperiment::counts(ref_data)
  SingleCellExperiment::altExp(ref_data, "ERCC_spikes") <- ref_data[grepl("^ERCC-", rownames(ref_data)), ]
  ## SCEList
  other_prior[["SCEList"]] <- ref_data
  if(is.null(other_prior[["write.noise.model"]])){
    other_prior[["write.noise.model"]] <- FALSE
  }
  estimate_formals <- simutils::change_parameters(function_expr = "BEARscc::estimate_noiseparameters",
                                                  other_prior = other_prior,
                                                  step = "estimation")
  ##############################################################################
  ####                            Estimation                                 ###
  ##############################################################################
  if(verbose){
    cat("Estimating parameters using BEARscc\n")
  }
  # Seed
  set.seed(seed)
  # Estimation
  tryCatch({
    # Estimate parameters from real data and return parameters and detection results
    estimate_detection <- peakRAM::peakRAM(
      estimate_result <- BiocGenerics::do.call(BEARscc::estimate_noiseparameters, estimate_formals)
    )
  }, error = function(e){
    as.character(e)
  })
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  estimate_output <- list(estimate_result = estimate_result,
                          estimate_detection = estimate_detection)
  return(estimate_output)
}


#' Simulate Datasets by BEARscc
#'
#' @param parameters A object generated by [BEARscc::simulate_replicates()]
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @param return_format A character. Alternative choices: list, SingleCellExperiment,
#' Seurat, h5ad. If you select `h5ad`, you will get a path where the .h5ad file saves to.
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#' @importFrom BEARscc simulate_replicates
#' @export
#'
#' @references
#' Severson D T, Owen R P, White M J, et al. BEARscc determines robustness of single-cell clusters using simulated technical replicates. Nature communications, 2018, 9(1): 1-7. <https://doi.org/10.1038/s41467-018-03608-y>
#'
#' Bioconductor URL: <https://www.bioconductor.org/packages/release/bioc/html/BEARscc.html>
#'
#' Github URL: <https://github.com/seversond12/BEARscc>
BEARscc_simulation <- function(parameters,
                               other_prior = NULL,
                               return_format,
                               verbose = FALSE,
                               seed
){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("BEARscc", quietly = TRUE)){
    cat("BEARscc is not installed on your device\n")
    cat("Installing BEARscc...\n")
    BiocManager::install("BEARscc")
  }
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  # Return to users
  other_prior[["n"]] <- 1
  other_prior[["SCEList"]] <- parameters
  cat(glue::glue("nCells: {ncol(parameters)}"), "\n")
  cat(glue::glue("nGenes: {nrow(parameters)}"), "\n")

  simulate_formals <- simutils::change_parameters(function_expr = "BEARscc::simulate_replicates",
                                                  other_prior = other_prior,
                                                  step = "simulation")
  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    cat("Simulating datasets using BEARscc\n")
  }
  # Seed
  set.seed(seed)
  # Estimation
  tryCatch({
    simulate_detection <- peakRAM::peakRAM(
      simulate_result <- do.call(BEARscc::simulate_replicates, simulate_formals)
    )
  }, error = function(e){
    as.character(e)
  })
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  counts <- simulate_result@metadata[["simulated_replicates"]][["Iteration_1"]]
  col_data <- data.frame("cell_name" = colnames(counts))
  row_data <- data.frame("gene_name" = rownames(counts))
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

