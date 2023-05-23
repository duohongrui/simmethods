#' Simulate Datasets by SimBPDD
#'
#' This function is used to simulate datasets `design_data` function in SimBPDD package.
#'
#' @param ref_data A matrix for one dataset or a list of datasets with their own
#' names. This is usually unused except for some methods, e.g. SCRIP, SimBPDD,
#' zingeR.
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @param return_format A character. Alternatives choices: list, SingleCellExperiment,
#' Seurat, h5ad. If you select `h5ad`, you will get a path where the .h5ad file saves to.
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#' @importFrom assertthat assert_that
#' @importFrom stringr str_extract
#' @importFrom dplyr mutate case_when
#' @export
#' @details
#' In addtion to simulate datasets with default parameters, users want to simulate
#' other kinds of datasets, e.g. a counts matrix with 2 or more cell groups. In
#' SimBPDD, you can set extra parameters to simulate datasets.
#'
#' The customed parameters you can set are below:
#' 1. nCells. In SimBPDD, you can set nCells directly `other_prior = list(nCells = 1000)`
#' to simulate 1000 cells.
#' 2. prob.group. You can directly set `other_prior = list(prob.group = c(0.4, 0.6))` to assign two proportions of cell groups. Note that the the length of the vector must be **2**.
#' For more customed parameters in SimBPDD, please check [SimBPDD::bp.sim.DD()].
#' @references
#' Schefzik R. SimBPDD: Simulating differential distributions in Beta-Poisson models, in particular for single-cell RNA sequencing data, Annales Mathematicae et Informaticae, 2021, 53: 283-298. <https://doi.org/10.33039/ami.2021.03.003>
#'
#' Github URL: <https://github.com/RomanSchefzik/SimBPDD>
#'
SimBPDD_simulation <- function(ref_data,
                               other_prior = NULL,
                               return_format,
                               verbose = FALSE,
                               seed
){
  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("SimBPDD", quietly = TRUE)){
    message("SimBPDD is not installed on your device...")
    message("Installing SimBPDD...")
    devtools::install_github("RomanSchefzik/SimBPDD")
  }
  ## nCells
  if(!is.null(other_prior[["nCells"]])){
    nCells <- other_prior[["nCells"]]
  }else{
    nCells <- ncol(ref_data)
  }
  ## prob.group
  if(!is.null(other_prior[["prob.group"]])){
    prob.group <- other_prior[["prob.group"]]
  }else{
    prob.group <- c(0.5, 0.5)
  }
  N1 <- round(nCells * prob.group[1])
  N2 <- nCells - N1
  ## degree
  if(!is.null(other_prior[["degree"]])){
    degree <- other_prior[["degree"]]
  }else{
    degree <- 1/2
  }
  ## case
  if(!is.null(other_prior[["case"]])){
    case <- other_prior[["case"]]
  }else{
    case <- "DBeta"
  }

  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  # Return to users
  message(paste0("nCells: ", nCells))
  message(paste0("nGenes: ", nrow(ref_data)))
  message("nGroups: 2")
  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    message("Simulating datasets using SimBPDD")
  }
  # Seed
  set.seed(seed)
  # Estimation
  simulate_detection <- peakRAM::peakRAM(
    simulate_result <- SimBPDD::bp.sim.DD(DATA = ref_data,
                                          N1 = N1,
                                          N2 = N2,
                                          case = case,
                                          degree = degree,
                                          seedex = seed)
  )
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  counts <- as.matrix(BiocGenerics::cbind(simulate_result[["controls"]][["samples.bp.wellfit"]],
                                          simulate_result[["manipulations"]]))
  colnames(counts) <- paste0("Cell", 1:ncol(counts))
  rownames(counts) <- paste0("Gene", 1:nrow(counts))
  # col_data
  col_data <- data.frame("cell_name" = colnames(counts),
                         "group" = c(rep("Group1", N1),
                                     rep("Group2", N2)))
  row_data <- data.frame("gene_name" = rownames(counts))
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

