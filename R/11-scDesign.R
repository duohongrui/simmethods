#' Simulate Datasets by scDesign
#'
#' @param ref_data A matrix for one dataset or a list of datasets with their own
#' names. This is usually unused except for some methods, e.g. SCRIP, scDesign,
#' zingeR.
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
#' @importFrom scDesign design_data
#' @importFrom assertthat assert_that
#' @importFrom glue glue
#' @importFrom SingleCellExperiment counts colData rowData SingleCellExperiment
#' @importFrom Seurat as.Seurat
#' @importFrom SeuratDisk SaveH5Seurat Convert
#' @importFrom stringr str_extract
#' @importFrom dplyr mutate case_when
#'
#' @export
#'
#'
scDesign_simulation <- function(ref_data,
                                other_prior,
                                return_format,
                                verbose = FALSE,
                                seed
){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("scDesign", quietly = TRUE)){
    cat("Splatter is not installed on your device\n")
    cat("Installing scDesign\n")
    devtools::install_github("Vivianstats/scDesign")
  }
  other_prior[["realcount"]] <- ref_data
  ## nCells
  if(!is.null(other_prior[["nCells"]])){
    other_prior[["ncell"]] <- ncol(ref_data)
  }else{
    other_prior[["ncell"]] <- other_prior[["nCells"]]
  }

  if(!is.null(other_prior[["de.prob"]])){
    if(length(other_prior[["de.prob"]]) == 1){
      other_prior[["pUp"]] <- ceiling(other_prior[["de.prob"]]/2)
      other_prior[["pDown"]] <- floor(other_prior[["de.prob"]]/2)
    }
    if(length(other_prior[["de.prob"]]) == 2){
      other_prior[["pUp"]] <- other_prior[["de.prob"]][1]
      other_prior[["pDown"]] <- other_prior[["de.prob"]][2]
    }
  }else{
    other_prior[["pUp"]] <- 0.05
    other_prior[["pDown"]] <- 0.05
  }

  if(!is.null(other_prior[["fc.group"]])){
    if(length(other_prior[["fc.group"]]) == 1){
      stop("When you are using scDesign to simulate groups, you should input fc.group with two numbers to indicate the upper fold change and lower fold change respectively.")
    }
    if(length(other_prior[["fc.group"]]) == 2){
      other_prior[["fU"]] <- other_prior[["fc.group"]][1]
      other_prior[["fL"]] <- other_prior[["fc.group"]][2]
    }
  }else{
    other_prior[["fU"]] <- 5
    other_prior[["fL"]] <- 1.5
  }

  simulate_formals <- as.list(formals(scDesign::design_data))
  for(param in names(simulate_formals)){
    names_wait_check <- names(other_prior)
    if(param %in% names_wait_check){
      simulate_formals[[param]] <- other_prior[[param]]
    }
  }

  if(is.null(other_prior[["nGroups"]])){
    other_prior[["nGroups"]] <- 1
  }

  if(other_prior[["nGroups"]] > 1){
    simulate_formals[["ngroup"]] <- other_prior[["nGroups"]]
    if(length(simulate_formals[["ncell"]]) != simulate_formals[["ngroup"]]){
      remainder <- simulate_formals[["ncell"]] %% simulate_formals[["ngroup"]]
      if(remainder == 0){
        simulate_formals[["ncell"]] <- rep(simulate_formals[["ncell"]]/simulate_formals[["ngroup"]],
                                           simulate_formals[["ngroup"]])
      }else{
        simulate_formals[["ncell"]] <- c(rep(simulate_formals[["ncell"]]%/%simulate_formals[["ngroup"]],
                                             simulate_formals[["ngroup"]]-1),
                                         simulate_formals[["ncell"]]%/%simulate_formals[["ngroup"]] + remainder)
      }
    }
    # assertthat::assert_that(length(other_prior[["nCells"]]) == other_prior[["nGroups"]],
    #                         msg = "The length of nCells must equal to nGroups")
    if(length(simulate_formals[["S"]]) != simulate_formals[["ngroup"]]){
      if(is.null(other_prior[["S"]])){
        simulate_formals[["S"]] <- rep(simulate_formals[["S"]], simulate_formals[["ngroup"]])
      }
      assertthat::assert_that(length(simulate_formals[["S"]]) == simulate_formals[["ngroup"]],
                              msg = "The length of S must equal to nGroups")
    }
  }

  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  # Return to users
  cat(glue::glue("nCells: {sum(simulate_formals[['ncell']])}"), "\n")
  cat(glue::glue("nGenes: {nrow(simulate_formals[['realcount']])}"), "\n")
  cat(glue::glue("nGroups: {simulate_formals[['ngroup']]}"), "\n")
  cat(glue::glue("de.prob: {simulate_formals[['pUp']] + simulate_formals[['pDown']]}"), "\n")
  cat(glue::glue("fc.group: up--{simulate_formals[['fU']]}"), "\n")
  cat(glue::glue("fc.group: down--{simulate_formals[['fL']]}"), "\n")
  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    cat("Simulating datasets using scDesign\n")
  }
  # Seed
  set.seed(seed)
  # Estimation
  tryCatch({
    simulate_detection <- peakRAM::peakRAM(
      simulate_result <- do.call(scDesign::design_data, simulate_formals)
    )
  }, error = function(e){
    as.character(e)
  })
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  if(simulate_formals[["ngroup"]] == 1){
    counts <- simulate_result
    rownames(counts) <- paste0("Gene", 1:nrow(counts))
    colnames(counts) <- paste0("Cell", 1:ncol(counts))
    col_data <- data.frame("cell_name" = colnames(counts))
    row_data <- data.frame("gene_name" = rownames(counts))
  }else{
    counts_tmp <- simulate_result[["count"]]

    # col_data
    col_data <- data.frame("cell_name" = paste0("Cell", 1:sum(simulate_formals[["ncell"]])),
                           "group" = paste0("Group", unlist(purrr::map(seq_len(length(simulate_formals[["ncell"]])),
                                                                       function(x){
                                                                         rep(as.character(x), ncol(counts_tmp[[x]]))
                                                                       }))))
    # Get count data together
    counts <- c()
    for(i in 1:length(counts_tmp)){
      counts <- cbind(counts, counts_tmp[[i]])
    }
    # Rename
    rownames(counts) <- paste0("Gene", 1:nrow(counts))
    colnames(counts) <- paste0("Cell", 1:ncol(counts))
    # Up and down regulated genes
    up_gene <- simulate_result[["genesUp"]][[2]] %>%
      stringr::str_extract(pattern = "[0-9]+") %>%
      as.numeric()
    down_gene <- simulate_result[["genesDown"]][[2]] %>%
      stringr::str_extract(pattern = "[0-9]+") %>%
      as.numeric()

    row_data <- data.frame("gene_name" = rownames(counts),
                           "de_genes" = FALSE,
                           "up_down" = "no")
    # Row data
    row_data[up_gene, 3] <- "up"
    row_data[down_gene, 3] <- "down"
    row_data <- row_data %>%
      dplyr::mutate("de_genes" = dplyr::case_when(
        up_down == "no" ~ FALSE,
        TRUE ~ TRUE
      ))
  }

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

