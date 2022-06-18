#' Simulate Datasets by scDesign
#'
#' @param ref_data A matrix for one dataset or a list of datasets with their own
#' names. This is usually unused except for some methods, e.g. SCRIP, scDesign.
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
  other_prior[["ncell"]] <- other_prior[["nCells"]]

  if(!is.null(other_prior[["prob.group"]])){
    if(length(other_prior[["prob.group"]]) == 1){
      other_prior[["pUp"]] <- ceiling(other_prior[["prob.group"]]/2)
      other_prior[["pDown"]] <- floor(other_prior[["prob.group"]]/2)
    }
    if(length(other_prior[["prob.group"]]) == 2){
      other_prior[["pUp"]] <- other_prior[["prob.group"]][1]
      other_prior[["pDown"]] <- other_prior[["prob.group"]][2]
    }
  }

  if(!is.null(other_prior[["fc.group"]])){
    if(length(other_prior[["fc.group"]]) == 1){
      stop("When you are using scDesign to simulate groups, you should input fc.group with two numbers to indicate the upper fold change and lower fold change respectively.")
    }
    if(length(other_prior[["fc.group"]]) == 2){
      other_prior[["fU"]] <- other_prior[["fc.group"]][1]
      other_prior[["fL"]] <- other_prior[["fc.group"]][2]
    }
  }

  simulate_formals <- as.list(formals(scDesign::design_data))
  for(param in names(simulate_formals)){
    names_wait_check <- names(other_prior)
    if(param %in% names_wait_check){
      simulate_formals[[param]] <- other_prior[[param]]
    }
  }

  if(other_prior[["nGroups"]] > 1){
    simulate_formals[["ngroup"]] <- other_prior[["nGroups"]]
    assertthat::assert_that(length(other_prior[["nCells"]]) == other_prior[["nGroups"]],
                            msg = "The length of nCells must equal to nGroups")
    if(length(simulate_formals[["S"]]) != simulate_formals[["ngroup"]]){
      assertthat::assert_that(length(other_prior[["nCells"]] == other_prior[["nGroups"]]),
                              msg = "The length of S must equal to nGroups")
    }
  }

  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  # Return to users
  cat(glue::glue("Your simulated datasets will have {simulate_formals[['ncell']]} cells, {dim(ref_data)[1]} genes, simulate_formals[['ngroup']] group(s) and contain {other_prior[['pUp']]+other_prior[['pDown']]} percent of DEGs"), "\n")
  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    cat("Simulating datasets using scDD\n")
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
  # Establish SingleCellExperiment
  simulate_result <- SingleCellExperiment::SingleCellExperiment(list(counts = counts),
                                                                colData = col_data,
                                                                rowData = row_data)
  if(return_format == "SingleCellExperiment"){
    simulate_result <- simulate_result
  }
  if(return_format == "list"){
    count_data <- SingleCellExperiment::counts(simulate_result)
    col_meta <- as.data.frame(SingleCellExperiment::colData(simulate_result))
    row_meta <- as.data.frame(SingleCellExperiment::rowData(simulate_result))
    simulate_result <- dplyr::lst(count_data,
                                  col_meta,
                                  row_meta)
  }
  if(return_format == "Seurat"){
    simulate_result <- Seurat::as.Seurat(simulate_result,
                                         counts = "counts",
                                         data = NULL)
  }
  if(return_format == "h5ad"){
    # Convert to Seurat object
    simulate_result <- Seurat::as.Seurat(simulate_result,
                                         counts = "counts",
                                         data = NULL)
    # Tmp file
    tmp_path <- tempdir()
    data_save_name <- file.path(tmp_path, paste0(time_string(), ".h5Seurat")) %>%
      simutils::fix_path()
    SeuratDisk::SaveH5Seurat(simulate_result, filename = data_save_name)
    SeuratDisk::Convert(data_save_name, dest = "h5ad")

    # data path
    data_path <- stringr::str_replace(data_save_name,
                                      pattern = "h5Seurat",
                                      replacement = "h5ad")
    cat(glue::glue("Your data has been save to {data_path}", "\n"))
    simulate_result <- list(file_type = "h5ad",
                            save_path = data_path)
  }

  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  simulate_output <- list(simulate_result = simulate_result,
                          simulate_detection = simulate_detection)
  return(simulate_output)
}

