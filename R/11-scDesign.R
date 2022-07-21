#' Simulate Datasets by scDesign
#'
#' @param ref_data A matrix for one dataset or a list of datasets with their own
#' names. This is usually unused except for some methods, e.g. SCRIP, scDesign,
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
#'
#' @importFrom scDesign design_data
#' @importFrom stringr str_extract
#' @importFrom dplyr mutate case_when
#'
#' @export
#' @details
#' In addtion to simulate datasets with default parameters, users want to simulate
#' other kinds of datasets, e.g. a counts matrix with 2 or more cell groups. In
#' scDesign, you can set extra parameters to simulate datasets.
#'
#' The customed parameters you can set are below:
#' 1. nCells. In scDesign, you can set nCells directly `other_prior = list(nCells = 1000)` to simulate 1000 cells.
#' 2. nGroups. You can directly set `other_prior = list(nGroups = 3)` to simulate 3 groups. But the cells will be assigned to these three groups equally if you do not set `prob.group` below.
#' 3. prob.group. You can directly set `other_prior = list(prob.group = c(0.2, 0.3, 0.5))` to assign three proportions of cell groups. Note that the number of groups always equals to the length of the vector.
#' 4. de.prob. You can directly set `other_prior = list(de.prob = 0.2)` to simulate DEGs that account for 20 percent of all genes.
#' 5. fc.group. You can directly set `other_prior = list(fc.group. = 2)` to specify the fold change of DEGs. But note that, you would better set `fc.group` because scDesign dose not retrun the fold changes of DEGs in the result.
#'
#' For more customed parameters in scDesign, please check [scDesign::design_data()].
#' @references
#' Li W V, Li J J. A statistical simulator scDesign for rational scRNA-seq experimental design[J]. Bioinformatics, 2019, 35(14): i41-i50. <https://doi.org/10.1093/bioinformatics/btz321>
#'
#' Github URL: <https://github.com/Vivianstats/scDesign>
#' @examples
#' ref_data <- simmethods::data
#'
#' ## Simulate datasets with default parameters
#' simulate_result <- scDesign_simulation(ref_data = ref_data,
#'                                        return_format = "list",
#'                                        verbose = TRUE,
#'                                        seed = 111)
#' counts <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(counts)
#'
#'
#' ## Simulate two groups with 20% proportion of DEGs and 2 fold change. Note that
#' ## scDesign does not provide fold changes for genes so users would better set
#' ## fc.group parameter in simulation function.
#' simulate_result <- scDesign_simulation(ref_data = ref_data,
#'                                        other_prior = list(nCells = 1000,
#'                                                           nGroups = 2,
#'                                                           de.prob = 0.2,
#'                                                           fc.group = 2),
#'                                        return_format = "list",
#'                                        verbose = TRUE,
#'                                        seed = 111)
#' counts <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(counts)
#' ## cell information
#' col_data <- simulate_result[["simulate_result"]][["col_meta"]]
#' table(col_data$group)
#' ## gene information
#' row_data <- simulate_result[["simulate_result"]][["row_meta"]]
#' table(row_data$de_genes)/4000
#' table(row_data$up_down)
#'
#'
#' ## Simulate three groups with 20% proportion of DEGs and 4 fold change. 20%, 40%
#' ## and 40% of cells belong to Group1, Group2 and Group3, respectively.
#' simulate_result <- scDesign_simulation(ref_data = ref_data,
#'                                        other_prior = list(nCells = 1000,
#'                                                           nGroups = 3,
#'                                                           prob.group = c(0.2, 0.4, 0.4),
#'                                                           de.prob = 0.2,
#'                                                           fc.group = 4),
#'                                        return_format = "list",
#'                                        verbose = TRUE,
#'                                        seed = 111)
#' counts <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(counts)
#' ## cell information
#' col_data <- simulate_result[["simulate_result"]][["col_meta"]]
#' table(col_data$group)
 scDesign_simulation <- function(ref_data,
                                 other_prior = NULL,
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
  if(is.null(other_prior[["nCells"]])){
    other_prior[["ncell"]] <- ncol(ref_data)
  }else{
    other_prior[["ncell"]] <- other_prior[["nCells"]]
  }

  if(!is.null(other_prior[["de.prob"]])){
    if(length(other_prior[["de.prob"]]) == 1){
      other_prior[["pUp"]] <- other_prior[["de.prob"]]/2
      other_prior[["pDown"]] <- other_prior[["de.prob"]]/2
    }
    if(length(other_prior[["de.prob"]]) == 2){
      other_prior[["pUp"]] <- other_prior[["de.prob"]][1]
      other_prior[["pDown"]] <- other_prior[["de.prob"]][2]
    }
  }

  if(!is.null(other_prior[["fc.group"]])){
    if(length(other_prior[["fc.group"]]) == 1){
      other_prior[["fU"]] <- other_prior[["fc.group"]]
      other_prior[["fL"]] <- other_prior[["fc.group"]]
    }
  }

  simulate_formals <- simutils::change_parameters(function_expr = "scDesign::design_data",
                                                  other_prior = other_prior,
                                                  step = "simulation")

  if(is.null(other_prior[["nGroups"]])){
    other_prior[["nGroups"]] <- 1
  }

  if(other_prior[["nGroups"]] > 1){
    simulate_formals[["ngroup"]] <- other_prior[["nGroups"]]
    if(length(simulate_formals[["ncell"]]) != simulate_formals[["ngroup"]]){
      ## group proportions
      if(!is.null(other_prior[["prob.group"]])){
        assert_that(other_prior[["nGroups"]] == length(other_prior[["prob.group"]]))
        simulate_formals[["ncell"]] <- simutils::proportionate(number = simulate_formals[["ncell"]],
                                                               result_sum_strict = simulate_formals[["ncell"]],
                                                               prop = other_prior[["prob.group"]],
                                                               prop_sum_strict = 1,
                                                               digits = 0)
      }else{
        simulate_formals[["ncell"]] <- simutils::proportionate(number = simulate_formals[["ncell"]],
                                                               result_sum_strict = simulate_formals[["ncell"]],
                                                               prop = rep(1/simulate_formals[["ngroup"]],
                                                                          simulate_formals[["ngroup"]]),
                                                               digits = 0)
      }
    }
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
                           "de_genes" = "no",
                           "up_down" = "no")
    # Row data
    row_data[up_gene, 3] <- "up"
    row_data[down_gene, 3] <- "down"
    row_data <- row_data %>%
      dplyr::mutate("de_genes" = dplyr::case_when(
        up_down == "no" ~ "no",
        TRUE ~ "yes"
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

