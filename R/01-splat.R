#' Estimate Parameters From Real Datasets by Splat
#'
#' This function is used to estimate useful parameters from a real dataset by
#' using \code{splatEstimate} function in Splatter package.
#'
#' @param ref_data A count matrix. Each row represents a gene and each column
#' represents a cell.
#' @param verbose Logical.
#' @param seed An integer of a random seed.
#' @importFrom peakRAM peakRAM
#' @importFrom splatter splatEstimate
#'
#' @return A list contains the estimated parameters and the results of execution
#' detection.
#' @export
#'
splat_estimation <- function(ref_data, verbose = FALSE, seed){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("splatter", quietly = TRUE)){
    cat("Splatter is not installed on your device\n")
    cat("Installing splatter...\n")
    BiocManager::install("splatter")
  }
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  if(!is.matrix(ref_data)){
    ref_data <- as.matrix(ref_data)
  }
  ##############################################################################
  ####                            Estimation                                 ###
  ##############################################################################
  if(verbose){
    cat("Estimating parameters using splat\n")
  }
  # Seed
  set.seed(seed)
  # Estimation
  tryCatch({
    # Estimate parameters from real data and return parameters and detection results
    estimate_detection <- peakRAM::peakRAM(estimate_result <- splatter::splatEstimate(ref_data))
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
