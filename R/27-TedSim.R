.mean_exp_group <- function(ref_data){
  per_cell <- data.frame('cell_group' = ref_data[["grouping"]],
                         'mean' = apply(ref_data[["counts"]], 1, mean))
  UMI_means <- stats::aggregate(per_cell$mean, list(per_cell$cell_group), mean)
  return(UMI_means)
}

#' Bundled TedSim Estimation Function
#'
#' @param ref_data Reference dataset
#' @param phyla Tree format data
#' @param seed Random seed
#' @import TedSim
#'
#' @return A list
#' @export
#'
TedSim_est <- function(
    ref_data,
    phyla,
    seed
){
  n_samples <- 1
  ncells <- length(ref_data[["cell_ids"]])
  N_nodes <- 2*ncells-2
  ngenes <- length(ref_data[["feature_ids"]])
  p_a <- 0.4
  n_cif <- 30
  n_diff <- 20
  cif_step <- 0.4

  #Initialize data and cell meta
  states_leaves_combined <- data.frame(parent = numeric(),
                                       cluster = character(),
                                       depth = numeric(),
                                       cellID = numeric())
  counts_combined <- c()
  cell_meta_combined <- c()
  scale_s_combined <- c()

  cifs_combined <- lapply(c(1:3),function(parami){
    matrix(ncol = n_cif,
           nrow = n_samples*ncells)
  })
  # Simulate shared State Identity Factors for all datasets
  SIF_res <- TedSim::SIFGenerate(phyla,
                                 n_diff,
                                 step = cif_step)
  set.seed(seed)
  cifs <- TedSim::SimulateCIFs(ncells,
                               phyla,
                               p_a = p_a,
                               n_CIF = n_cif,
                               n_diff = n_diff,
                               step = cif_step,
                               Sigma = 0.5,
                               SIF_res = SIF_res,
                               unif_on = FALSE)
  for (parami in c(1:3)){
    cif_leaves_all <- cifs[[1]][[parami]][c(1:ncells),]
    cifs_combined[[parami]][1:ncells ,] <- cif_leaves_all
  }
  cell_meta_combined <- rbind(cell_meta_combined, cifs[[2]][1:ncells,])
  states <- cifs[[2]]
  states <- states[1:N_nodes,]
  states_leaves <- states[1:ncells,]
  means <- .mean_exp_group(ref_data = ref_data)
  states_uniq <- c(match(means$Group.1, phyla$tip.label),
                   (length(phyla$tip.label)+1):(length(phyla$tip.label)+phyla$Nnode))
  scale_s_states <- c(0.03 * means$x,rep(0.001,phyla$Nnode))
  scale_s <- states_leaves[,2]
  scale_s[scale_s %in% states_uniq] <- scale_s_states[match(scale_s, states_uniq, nomatch = 0)]
  state_intermediate <- c(match(means$Group.1, phyla$tip.label),
                          (length(phyla$tip.label)+1):(length(phyla$tip.label)+phyla$Nnode))
  state_merge_intermediate <- c(means$Group.1, rep("intermediate",phyla$Nnode))
  new_states <- states_leaves[,2]
  new_states[new_states %in% state_intermediate] <- state_merge_intermediate[match(new_states,
                                                                                   state_intermediate,
                                                                                   nomatch = 0)]
  states_leaves[, 2] <- new_states
  scale_s_combined <- c(scale_s_combined, scale_s)
  states_leaves_combined <- rbind(states_leaves_combined, states_leaves)
  cifs_combined <- list(cifs_combined, cell_meta_combined)
  return(list(cifs_combined = cifs_combined,
              scale_s_combined = scale_s_combined))
}


#' Estimate Parameters From Real Datasets by TedSim
#'
#' This function is used to estimate useful parameters from a real dataset by
#' using \code{TedSim_est} function in simmethods package.
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
#' @importFrom simutils make_trees
#' @importFrom dynwrap wrap_expression
#'
#' @return A list contains the estimated parameters and the results of execution
#' detection.
#' @export
#' @details
#' In TedSim, users can input cell group information if it is available. If cell
#' group information is not provided, the procedure will detect cell groups by
#' kmeans automatically.
#' See `Examples` for more instructions.
#'
#' Note that TedSim can only simulate the dataset whose cell number is the power of 2
#' so if the reference data does not meet the requirement, the procedure will
#' synthesize extra fake cells to achive this goal. See [simutils::synthesize_cells()]
#' to learn about the process of synthesizing fake cells.
#'
#' @references
#' Pan X, Li H, Zhang X. TedSim: temporal dynamics simulation of single-cell RNA sequencing data and cell division history. Nucleic acids research, 2022, 50(8): 4272-4288. <https://doi.org/10.1093/nar/gkac235>
#'
#' Github URL: <https://github.com/Galaxeee/TedSim>
#'
#' @examples
#' ref_data <- simmethods::data
#' ## estimation
#' estimate_result <- simmethods::TedSim_estimation(
#'   ref_data = ref_data,
#'   other_prior = NULL,
#'   verbose = TRUE,
#'   seed = 111
#' )
#'
#' ## estimation with cell group information
#' group_condition <- paste0("Group", as.numeric(simmethods::group_condition))
#' estimate_result <- simmethods::TedSim_estimation(
#'   ref_data = ref_data,
#'   other_prior = list(group.condition = group_condition),
#'   verbose = TRUE,
#'   seed = 111
#' )
TedSim_estimation <- function(ref_data,
                              other_prior = NULL,
                              verbose = FALSE,
                              seed
){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("TedSim", quietly = TRUE)){
    message("TedSim is not installed on your device")
    message("Installing TedSim...")
    devtools::install_github("Galaxeee/TedSim")
  }
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  if(!is.matrix(ref_data)){
    ref_data <- as.matrix(ref_data)
  }
  if(!is.null(other_prior[["group.condition"]])){
    group <- other_prior[["group.condition"]]
  }else{
    group <- NULL
  }
  if(log2(ncol(ref_data)) != as.integer(log2(ncol(ref_data)))){
    message("The number of cells is not the power of 2, and we will synthesize some extra cells base on your data...")
    ref_data <- simutils::synthesize_cells(ref_data,
                                           group = group,
                                           seed = seed,
                                           verbose = verbose)
  }else{
    ref_data <- dynwrap::wrap_expression(counts = t(ref_data),
                                         expression = log2(t(ref_data) + 1),
                                         grouping = group)
  }
  phyla <- simutils::make_trees(ref_data = t(as.matrix(ref_data[["expression"]])),
                                group = ref_data[["grouping"]],
                                is_Newick = FALSE,
                                is_parenthetic = FALSE)
  ##############################################################################
  ####                            Estimation                                 ###
  ##############################################################################
  if(verbose){
    message("Estimating parameters using TedSim")
  }
  # Seed
  set.seed(seed)
  # Estimation
  estimate_detection <- peakRAM::peakRAM(
    estimate_result <- simmethods::TedSim_est(ref_data = ref_data,
                                              phyla = phyla,
                                              seed = seed)
  )
  estimate_result <- list(estimate_result = estimate_result,
                          data_dim = c(length(ref_data[["feature_ids"]]),
                                       length(ref_data[["cell_ids"]])))
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  estimate_output <- list(estimate_result = estimate_result,
                          estimate_detection = estimate_detection)
  return(estimate_output)
}


#' Simulate Datasets by TedSim
#'
#' @param parameters A object generated by [simmethods::TedSim_est()]
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs are usually customed, so before simulating a dataset you must point it out.
#' See `Details` below for more information.
#' @param return_format A character. Alternative choices: list, SingleCellExperiment,
#' Seurat, h5ad. If you select `h5ad`, you will get a path where the .h5ad file saves to.
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#' @importFrom glue glue
#' @importFrom stringr str_split str_count str_extract_all str_replace
#' @importFrom reticulate source_python
#' @export
#' @details
#' In TedSim, users can only set `nGenes` to specify the number of genes in the
#' simulated dataset. See `Examples` for instructions.
#'
#' @references
#' Pan X, Li H, Zhang X. TedSim: temporal dynamics simulation of single-cell RNA sequencing data and cell division history. Nucleic acids research, 2022, 50(8): 4272-4288. <https://doi.org/10.1093/nar/gkac235>
#'
#' Github URL: <https://github.com/Galaxeee/TedSim>
#'
#' @examples
#' ref_data <- simmethods::data
#'
#' ## estimation with cell group information
#' group_condition <- paste0("Group", as.numeric(simmethods::group_condition))
#' estimate_result <- simmethods::TedSim_estimation(
#'   ref_data = ref_data,
#'   other_prior = list(group.condition = group_condition),
#'   verbose = TRUE,
#'   seed = 111
#' )
#'
#' # 1) Simulate with default parameters
#' simulate_result <- simmethods::TedSim_simulation(
#'   parameters = estimate_result[["estimate_result"]],
#'   other_prior = NULL,
#'   return_format = "list",
#'   verbose = TRUE,
#'   seed = 111
#' )
#' ## counts
#' counts <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(counts)
#'
#' # 2) 5000 genes
#' simulate_result <- simmethods::TedSim_simulation(
#'   parameters = estimate_result[["estimate_result"]],
#'   other_prior = list(nGenes = 5000),
#'   return_format = "list",
#'   verbose = TRUE,
#'   seed = 111
#' )
#'
#' ## counts
#' counts <- simulate_result[["simulate_result"]][["count_data"]]
#' dim(counts)
TedSim_simulation <- function(parameters,
                              other_prior = NULL,
                              return_format,
                              verbose = FALSE,
                              seed
){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("TedSim", quietly = TRUE)){
    message("TedSim is not installed on your device")
    message("Installing TedSim...")
    devtools::install_github("Galaxeee/TedSim")
  }
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  other_prior[["ncif"]] <- 30
  other_prior[["cif_res"]] <- parameters[["estimate_result"]][["cifs_combined"]]
  other_prior[["scale_s"]] <- parameters[["estimate_result"]][["scale_s_combined"]]
  # nCells (equal to ref data)
  other_prior[["ncells"]] <- parameters[["data_dim"]][2]
  # nGenes
  if(!is.null(other_prior[["nGenes"]])){
    other_prior[["ngenes"]] <- other_prior[["nGenes"]]
  }else{
    other_prior[["ngenes"]] <- parameters[["data_dim"]][1]
  }
  # de.prob
  if(!is.null(other_prior[["de.prob"]])){
    other_prior[["ge_prob"]] <- other_prior[["de.prob"]]
  }else{
    other_prior[["ge_prob"]] <- 0.1
  }

  simulate_formals <- simutils::change_parameters(function_expr = "TedSim::CIF2Truecounts",
                                                  other_prior = other_prior,
                                                  step = "simulation")
  # Return to users
  message(glue::glue("nCells: {other_prior[['ncells']]}"))
  message(glue::glue("nGenes: {other_prior[['ngenes']]}"))
  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    message("Simulating datasets using TedSim")
  }
  # Seed
  set.seed(seed)
  # Estimation
  simulate_detection <- peakRAM::peakRAM(
    simulate_result <- do.call(TedSim::CIF2Truecounts, simulate_formals)
  )
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  col_meta <- as.data.frame(simulate_result[["cell_meta"]])
  simulate_result <- simulate_result[["counts"]]
  colnames(simulate_result) <- paste0("Cell", 1:ncol(simulate_result))
  rownames(simulate_result) <- paste0("Gene", 1:nrow(simulate_result))
  ## col_data
  group <- as.numeric(as.factor(as.numeric(col_meta$cluster)))
  col_data <- data.frame("cell_name" = colnames(simulates_result),
                         "group" = paste0("Group", group))
  ## row_data
  row_data <- data.frame("gene_name" = rownames(simulate_result))
  # Establish SingleCellExperiment
  simulate_result <- SingleCellExperiment::SingleCellExperiment(list(counts = simulate_result),
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

