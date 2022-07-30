#' Estimate Parameters From Real Datasets by ESCO
#'
#' This function is used to estimate useful parameters from a real dataset by
#' using \code{escoEstimate} function in ESCO package.
#'
#' @param ref_data A count matrix. Each row represents a gene and each column
#' represents a cell.
#' @param verbose Logical.
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @param seed An integer of a random seed.
#' @importFrom ESCO escoEstimate newescoParams
#' @return A list contains the estimated parameters and the results of execution
#' detection.
#' @export
#' @details
#' In ESCO, users can input cell group information when it is available but in this case
#' ESCO is not stable and may fail to estimate suitable distribution parameters
#' from real data.
#' If users want to estimate tree structured parameters, set `other_prior = list(tree = TRUE)`.
#' For more instructions, see `Examples`.
#' @references
#' Tian J, Wang J, Roeder K. ESCO: single cell expression simulation incorporating gene co-expression. Bioinformatics, 2021, 37(16): 2374-2381. <https://doi.org/10.1093/bioinformatics/btab116>
#'
#' Github URL: <https://github.com/JINJINT/ESCO>
#'
#' @examples
#' ref_data <- simmethods::data
#'
#' estimate_result <- simmethods::ESCO_estimation(ref_data = ref_data,
#'                                                other_prior = NULL,
#'                                                verbose = TRUE,
#'                                                seed = 111)
#' # If cell group information is available, it can be another prior information.
#' # But there is a bug in ESCO, and some datasets can not be estimated due to the
#' # failing estimation of distribution parameters.
#' # group_condition <- as.numeric(simmethods::group_condition)
#' # estimate_result <- simmethods::ESCO_estimation(
#' #   ref_data = ref_data,
#' #   other_prior = list(group.condition = group_condition),
#' #   verbose = TRUE,
#' #   seed = 111
#' # )
#'
#' # ----------------- Estimate tree or trajectory structured data -------------
#' # Load data
#' ref_data <- simmethods::data
#' # Estimate parameters
#' estimate_result <- simmethods::ESCO_estimation(ref_data = ref_data,
#'                                                other_prior = list(tree = TRUE),
#'                                                verbose = TRUE,
#'                                                seed = 10)
ESCO_estimation <- function(ref_data,
                            other_prior = NULL,
                            verbose = FALSE,
                            seed){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("ESCO", quietly = TRUE)){
    cat("ESCO is not installed on your device\n")
    cat("Installing ESCO...\n")
    devtools::install_github("JINJINT/ESCO")
  }
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  if(!is.matrix(ref_data)){
    ref_data <- as.matrix(ref_data)
  }
  other_prior[["counts"]] <- ref_data
  if(is.null(other_prior[["group.condition"]])){
    other_prior[["group"]] <- FALSE
    other_prior[["cellinfo"]] <- NULL
  }else{
    other_prior[["group"]] <- TRUE
    other_prior[["cellinfo"]] <- other_prior[["group.condition"]]
  }
  other_prior[["params"]] <- ESCO::newescoParams()
  other_prior[["dirname"]] <- tempdir()
  ## tree
  if(!is.null(other_prior[["tree"]])){
    tree <- simutils::make_trees(ref_data = ref_data,
                                 group = other_prior[["cellinfo"]],
                                 is_Newick = FALSE,
                                 is_parenthetic = TRUE,
                                 return_group = TRUE)
    group <- tree[["group"]]
    tree <- tree[["phyla"]]
  }else{
    tree <- NULL
  }
  estimate_formals <- simutils::change_parameters(function_expr = "ESCO::escoEstimate",
                                                  other_prior = other_prior,
                                                  step = "estimation")
  ##############################################################################
  ####                            Estimation                                 ###
  ##############################################################################
  if(verbose){
    cat("Estimating parameters using ESCO\n")
  }
  # Seed
  set.seed(seed)
  # Estimation
  tryCatch({
    # Estimate parameters from real data and return parameters and detection results
    estimate_detection <- peakRAM::peakRAM(
      estimate_result <- do.call(ESCO::escoEstimate, estimate_formals)
    )
  }, error = function(e){
    as.character(e)
  })
  if(!is.null(tree)){
    estimate_result <- list(estimate_result = estimate_result,
                            tree = tree,
                            group = group)
  }
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  estimate_output <- list(estimate_result = estimate_result,
                          estimate_detection = estimate_detection)
  return(estimate_output)
}


#' Simulate Datasets by ESCO
#'
#' @param parameters A object generated by [ESCO::escoEstimate()]
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @param return_format A character. Alternative choices: list, SingleCellExperiment,
#' Seurat, h5ad. If you select `h5ad`, you will get a path where the .h5ad file saves to.
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#' @importFrom ESCO escoEstimate
#' @export
#' @details
#' In addtion to simulate datasets with default parameters, users want to simulate
#' other kinds of datasets, e.g. a counts matrix with 2 or more cell groups. In
#' ESCO, you can set extra parameters to simulate datasets.
#'
#' The customed parameters you can set are below:
#' 1. nCells. You can directly set `other_prior = list(nCells = 1000)` to simulate 1000 cells.
#' 2. nGenes. You can directly set `other_prior = list(nGenes = 5000)` to simulate 5000 genes.
#' 3. nGroups. You can not directly set `other_prior = list(nGroups = 3)` to simulate 3 groups. Instead, you should set `other_prior = list(prob.group = c(0.2, 0.3, 0.5))` where the sum of group probabilities must equal to 1.
#' 4. de.prob. You can directly set `other_prior = list(de.prob = 0.2)` to simulate DEGs that account for 20 percent of all genes.
#' 5. prob.group. You can directly set `other_prior = list(prob.group = c(0.2, 0.3, 0.5))` to assign three proportions of cell groups. Note that the number of groups always equals to the length of the vector.
#' 6. If users want to simulate tree or trajectory structured data, set `other_prior = list(type = "tree")` or `other_prior = list(type = "traj")`. See `Examples`.
#'
#' For more customed parameters in ESCO, please check [ESCO::escoParams()].
#' @references
#' Tian J, Wang J, Roeder K. ESCO: single cell expression simulation incorporating gene co-expression. Bioinformatics, 2021, 37(16): 2374-2381. <https://doi.org/10.1093/bioinformatics/btab116>
#'
#' Github URL: <https://github.com/JINJINT/ESCO>
#'
#' @examples
#' ## Estimation
#' ref_data <- simmethods::data
#'
#' estimate_result <- simmethods::ESCO_estimation(ref_data = ref_data,
#'                                                other_prior = NULL,
#'                                                verbose = TRUE,
#'                                                seed = 111)
#' # 1) Simulate with default parameters
#' simulate_result <- simmethods::ESCO_simulation(
#'   parameters = estimate_result[["estimate_result"]],
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
#'
#'
#' # 2) Simulate two groups (20% proportion of DEGs)
#' simulate_result <- simmethods::ESCO_simulation(
#'   parameters = estimate_result[["estimate_result"]],
#'   other_prior = list(nCells = 1000,
#'                      nGenes = 2000,
#'                      de.prob = 0.2,
#'                      prob.group = c(0.3, 0.7)),
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
#' table(col_data$group)/1000
#' ## gene information
#' row_data <- simulate_result[["simulate_result"]][["row_meta"]]
#' table(row_data$de_gene)[2]/2000
#' ### The result of ESCO contains the factors of different groups and uses can
#' ### calculate the fold change by division. For example, the DEFactors of A gene
#' ### in Group1 and Group2 are respectively 2 and 1, and the fold change of A gene
#' ### is 2/1=2 or 1/2=0.5.
#' fc_group1_to_group2 <- row_data$DEFacGroup2/row_data$DEFacGroup1
#' table(fc_group1_to_group2 != 1)[2]/2000 ## de.prob = 0.2
#'
#' ------------------ Simulate tree or trajectory structured data --------------
#' # Load data
#' ref_data <- simmethods::data
#' # Estimate parameters
#' estimate_result <- simmethods::ESCO_estimation(ref_data = ref_data,
#'                                                other_prior = list(tree = TRUE),
#'                                                verbose = TRUE,
#'                                                seed = 10)
#'
#' # (1) Simulate tree structured cell groups
#' simulate_result <- simmethods::ESCO_simulation(parameters = estimate_result[["estimate_result"]],
#'                                                other_prior = list(batchCells = 1000,
#'                                                                   nGenes = 1000,
#'                                                                   type = "tree"),
#'                                                return_format = "SingleCellExperiment",
#'                                                verbose = TRUE,
#'                                                seed = 111)
#' ## plot
#' result <- scater::logNormCounts(simulate_result[["simulate_result"]])
#' result <- scater::runPCA(result)
#' plotPCA(result, colour_by = "group")
#'
#'
#' # (2) Simulate tree structured cell groups (specify de.prob and group.prob)
#' simulate_result <- simmethods::ESCO_simulation(parameters = estimate_result[["estimate_result"]],
#'                                                other_prior = list(batchCells = 1000,
#'                                                                   nGenes = 1000,
#'                                                                   type = "tree",
#'                                                                   group.prob = c(0.4, 0.3, 0.3),
#'                                                                   de.prob = 0.5,
#'                                                                   de.center = 2),
#'                                                return_format = "SingleCellExperiment",
#'                                                verbose = TRUE,
#'                                                seed = 111)
#' ## plot
#' result <- scater::logNormCounts(simulate_result[["simulate_result"]])
#' result <- scater::runPCA(result)
#' plotPCA(result, colour_by = "group")
#'
#'
#' # (3) Simulate continous cell trajectory
#' simulate_result <- simmethods::ESCO_simulation(parameters = estimate_result[["estimate_result"]],
#'                                                other_prior = list(batchCells = 1000,
#'                                                                   nGenes = 1000,
#'                                                                   type = "traj",
#'                                                                   group.prob = c(0.4, 0.3, 0.3),
#'                                                                   de.prob = 0.5,
#'                                                                   de.center = 2),
#'                                                return_format = "SingleCellExperiment",
#'                                                verbose = TRUE,
#'                                                seed = 111)
#' ## plot
#' result <- scater::logNormCounts(simulate_result[["simulate_result"]])
#' result <- scater::runPCA(result)
#' plotPCA(result, colour_by = "group")
ESCO_simulation <- function(parameters,
                            return_format,
                            other_prior = NULL,
                            verbose = FALSE,
                            seed
){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("ESCO", quietly = TRUE)){
    cat("ESCO is not installed on your device\n")
    cat("Installing ESCO...\n")
    devtools::install_github("JINJINT/ESCO")
  }
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  if(length(parameters) == 3){
    tree <- parameters[["tree"]]
    group <- parameters[["group"]]
    parameters <- parameters[["estimate_result"]]
    type <- other_prior[["type"]]
    if(is.null(type)){
      stop("Please input the type ('tree' or 'traj') when you want to simulate dataset with tree format information.")
    }
  }else{
    group <- NULL
    tree <- NULL
    type <- NULL
  }
  assertthat::assert_that(class(parameters) == "escoParams")
  if(!is.null(other_prior)){
    parameters <- simutils::set_parameters(parameters = parameters,
                                           other_prior = other_prior,
                                           method = "ESCO")
  }
  # nCells
  if(!is.null(other_prior[["nCells"]])){
    parameters <- splatter::setParam(parameters, name = "nCells", value = other_prior[["nCells"]])
  }
  # nGenes
  if(!is.null(other_prior[["nGenes"]])){
    parameters <- splatter::setParam(parameters, name = "nGenes", value = other_prior[["nGenes"]])
  }
  # nGroup
  if(!is.null(group)){
    parameters <- splatter::setParam(parameters, name = "nGroups", value = length(unique(group)))
  }
  if(!is.null(other_prior[["prob.group"]])){
    parameters <- splatter::setParam(parameters,
                                     name = "group.prob",
                                     value = other_prior[["prob.group"]])

  }else{
    if(!is.null(type)){
      nGroups <- splatter::getParam(parameters, name = "nGroups")
      prob.group <- c(rep(1/nGroups, nGroups-1),
                      1 - c(1/nGroups * c(nGroups-1)))
      parameters <- splatter::setParam(parameters,
                                       name = "group.prob",
                                       value = prob.group)
    }
  }
  # Get params to check
  params_check <- splatter::getParams(parameters, c("nCells",
                                                    "nGenes",
                                                    "nGroups",
                                                    "group.prob",
                                                    "de.prob"))

  # Return to users
  cat(glue::glue("nCells: {params_check[['nCells']]}"), "\n")
  cat(glue::glue("nGenes: {params_check[['nGenes']]}"), "\n")
  cat(glue::glue("nGroups: {params_check[['nGroups']]}"), "\n")
  cat(glue::glue("de.group: {params_check[['de.prob']]}"), "\n")
  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    cat("Simulating datasets using ESCO\n")
  }
  # Seed
  parameters <- splatter::setParam(parameters, name = "seed", value = seed)
  parameters <- splatter::setParam(parameters,
                                   name = "deall.prob",
                                   value = params_check[['de.prob']])
  # Estimation
  tryCatch({
    if(!is.null(type)){
      if(type == "tree"){
        cat("Simulating trajectory of trees datasets by ESCO \n")
        submethod <- "tree"
      }
      if(type == "traj"){
        cat("Simulating trajectory datasets by ESCO \n")
        submethod <- "traj"
      }
      parameters <- splatter::setParam(parameters, name = "tree", value = tree)
    }else{
      if(params_check[["nGroups"]] == 1){
        submethod <- "single"
      }else if(params_check[["nGroups"]] != 1){
        submethod <- "group"
      }
    }
    simulate_detection <- peakRAM::peakRAM(
      simulate_result <- ESCO::escoSimulate(parameters,
                                            type = submethod,
                                            verbose = verbose))
  }, error = function(e){
    as.character(e)
  })
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  # counts
  counts <- as.matrix(SingleCellExperiment::counts(simulate_result))
  # col_data
  col_data <- as.data.frame(SummarizedExperiment::colData(simulate_result))
  if(params_check[['nGroups']] == 1){
    col_data <- data.frame("cell_name" = colnames(counts),
                           "group" = rep("Group1", ncol(counts)))
  }else{
    if(type == "traj"){
      col_data <- data.frame("cell_name" = colnames(counts),
                             "group" = paste0("Group", col_data$Path))
    }else{
      col_data <- data.frame("cell_name" = colnames(counts),
                             "group" = col_data$Group)
    }
  }
  rownames(col_data) <- col_data$cell_name
  # row_data
  row_data <- as.data.frame(SummarizedExperiment::rowData(simulate_result))
  if(params_check[['nGroups']] == 1 | !is.null(parameters@tree)){
    row_data <- data.frame("gene_name" = rownames(counts))
    rownames(row_data) <- row_data$gene_name
  }else{
    group_fac <- row_data[, grep(colnames(row_data), pattern = "^DEFac")]
    total_sum <- rowSums(group_fac)
    de_gene <- ifelse(total_sum == params_check[['nGroups']], "no", "yes")
    row_data[, 2] <- de_gene
    row_data <- row_data[, -c(3:5, ncol(row_data))]
    colnames(row_data) <- c("gene_name", "de_gene", colnames(group_fac))
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

