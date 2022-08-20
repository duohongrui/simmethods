#' Estimate Parameters From Real Datasets by BASiCS
#'
#' This function is used to estimate useful parameters from a real dataset by
#' using `BASiCSEstimate` function in Splatter package.
#'
#' @param ref_data A count matrix. Each row represents a gene and each column
#' represents a cell.
#' @param verbose Logical.
#' @param seed An integer of a random seed.
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @importFrom splatter BASiCSEstimate
#' @return A list contains the estimated parameters and the results of execution
#' detection.
#' @export
#'
#' @references
#' Zappia L, Phipson B, Oshlack A. Splatter: simulation of single-cell RNA sequencing data. Genome biology, 2017, 18(1): 1-15. <https://doi.org/10.1186/s13059-017-1305-0>
#'
#' Bioconductor URL: <https://bioconductor.org/packages/release/bioc/html/splatter.html>
#'
#' Github URL: <https://github.com/Oshlack/splatter>
BASiCS_estimation <- function(ref_data,
                              verbose = FALSE,
                              other_prior = NULL,
                              seed
){
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  if(!is.matrix(ref_data)){
    ref_data <- as.matrix(ref_data)
  }

  other_prior[["params"]] <- splatter::newBASiCSParams()
  if(is.null(other_prior[["dilution.factor"]]) & is.null(other_prior[["batch.condition"]])){
    stop("At least one of spike.info and batch must be provided.")
  }
  if(!is.null(other_prior[["dilution.factor"]]) & !is.null(other_prior[["volume"]])){
    ERCC_index <- grep(rownames(ref_data), pattern = "^ERCC")
    ERCC_counts <- ref_data[ERCC_index, ]
    cell_remove_index <- colnames(ERCC_counts)[colSums(ERCC_counts) == 0]
    warning(glue::glue("These cells have zero counts in spike-in genes and will be moved: {paste0(cell_remove_index, collapse = ', ')}"))
    spikeData <- ref_data[grep(rownames(ref_data), pattern = "^ERCC"), ]
    concentration <- simmethods::ERCC_info$con_Mix1_attomoles_ul
    spikeInfo <- data.frame(Name = simmethods::ERCC_info$ERCC_id,
                            Input = concentration*10^-18*6.022*10^23*other_prior[["volume"]]/other_prior[["dilution.factor"]],
                            row.names = simmethods::ERCC_info$ERCC_id)
    spikeInfo <- spikeInfo[rownames(spikeData), ]
    other_prior[["spike.info"]] <- spikeInfo
    other_prior[["counts"]] <- ref_data[, -which(colnames(ref_data) %in% cell_remove_index)]
  }
  if(!is.null(other_prior[["batch.condition"]])){
    other_prior[["batch"]] <- other_prior[["batch.condition"]]
  }
  estimate_formals <- simutils::change_parameters(function_expr = "splatter::BASiCSEstimate",
                                                  other_prior = other_prior,
                                                  step = "estimation")
  ##############################################################################
  ####                            Estimation                                 ###
  ##############################################################################
  if(verbose){
    message("Estimating parameters using BASiCS")
  }
  # Seed
  set.seed(seed)
  # Estimation
  tryCatch({
    # Estimate parameters from real data and return parameters and detection results
    estimate_detection <- peakRAM::peakRAM(
      estimate_result <- splatter::BASiCSEstimate(counts = estimate_formals[["counts"]],
                                                  batch = estimate_formals[["batch"]],
                                                  spike.info = estimate_formals[["spike.info"]],
                                                  n = estimate_formals[["n"]],
                                                  thin = estimate_formals[["thin"]],
                                                  burn = estimate_formals[["burn"]],
                                                  regression = estimate_formals[["regression"]],
                                                  params = estimate_formals[["params"]],
                                                  verbose = estimate_formals[["verbose"]],
                                                  progress = estimate_formals[["progress"]])
    )
  }, error = function(e){
    print(e)
  })
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  estimate_output <- list(estimate_result = estimate_result,
                          estimate_detection = estimate_detection)
  return(estimate_output)
}



#' Simulate Datasets by BASiCS
#'
#' This function is used to simulate datasets from learned parameters by `BASiCSSimulate`
#' function in Splatter package.
#'
#' @param parameters A object generated by [splatter::BASiCSEstimate()]
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @param return_format A character. Alternative choices: list, SingleCellExperiment,
#' Seurat, h5ad. If you select `h5ad`, you will get a path where the .h5ad file saves to.
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#' @importFrom splatter BASiCSSimulate
#' @export
#'
#' @references
#' Zappia L, Phipson B, Oshlack A. Splatter: simulation of single-cell RNA sequencing data. Genome biology, 2017, 18(1): 1-15. <https://doi.org/10.1186/s13059-017-1305-0>
#'
#' Bioconductor URL: <https://bioconductor.org/packages/release/bioc/html/splatter.html>
#'
#' Github URL: <https://github.com/Oshlack/splatter>
BASiCS_simulation <- function(parameters,
                              other_prior = NULL,
                              return_format,
                              verbose = FALSE,
                              seed
){
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  assertthat::assert_that(class(parameters) == "BASiCSParams")
  if(!is.null(other_prior)){
    parameters <- simutils::set_parameters(parameters = parameters,
                                           other_prior = other_prior,
                                           method = "BASiCS")
  }
  # batchCells
  if(!is.null(other_prior[["batchCells"]])){
    if(length(other_prior[["batchCells"]]) != length(parameters@theta)){
      stop("The length of batchCells must equal to the theta in parameters")
    }
    parameters <- splatter::setParam(parameters, name = "batchCells", value = other_prior[["batchCells"]])
  }
  # nGenes
  if(!is.null(other_prior[["nGenes"]])){
    parameters <- splatter::setParam(parameters, name = "nGenes", value = other_prior[["nGenes"]])
  }
  # Get params to check
  params_check <- splatter::getParams(parameters, c("nCells",
                                                    "nGenes",
                                                    "nBatches"))
  # Return to users
  message(glue::glue("nCells: {params_check[['nCells']]}"))
  message(glue::glue("nGenes: {params_check[['nGenes']]}"))
  message(glue::glue("nBatches: {params_check[['nBatches']]}"))
  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    message("Simulating datasets using BASiCS")
  }
  # Seed
  parameters <- splatter::setParam(parameters, name = "seed", value = seed)
  # Estimation
  tryCatch({
    simulate_detection <- peakRAM::peakRAM(
      simulate_result <- splatter::BASiCSSimulate(parameters, verbose = verbose)
    )
  }, error = function(e){
    print(e)
  })
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  # counts
  counts <- as.matrix(SingleCellExperiment::counts(simulate_result))
  # col_data
  col_data <- as.data.frame(SummarizedExperiment::colData(simulate_result))
  col_data <- col_data[, c("Cell", "Batch")]
  colnames(col_data) <- c("cell_name", "batch")
  col_data$batch <- paste0("Batch", col_data$batch)
  # row_data
  row_data <- data.frame("gene_name" = paste0("Gene", 1:nrow(counts)))
  rownames(row_data) <- row_data$gene_name
  simulate_result <- simutils::data_conversion(SCE_object = simulate_result,
                                               return_format = return_format)

  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  simulate_output <- list(simulate_result = simulate_result,
                          simulate_detection = simulate_detection)
  return(simulate_output)
}

