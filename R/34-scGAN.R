#' Estimate Parameters From Real Datasets by scGAN
#'
#' This function is used to estimate useful parameters in a Docker container.
#'
#' @param ref_data A count matrix. Each row represents a gene and each column
#' represents a cell.
#' @param verbose Logical.
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. In
#' simulation step, the number of cells, genes, groups, batches, the percent of
#' DEGs and other variables are usually customed, so before simulating a dataset
#' you must point it out.
#' @param seed An integer of a random seed.
#' @param ... Other attributes and new values to estimate parameters using scGAN
#' @importFrom dynwrap test_docker_installation
#' @importFrom simutils fix_path change_scGAN_parameters data_conversion
#' @importFrom babelwhale list_docker_images pull_container
#' @importFrom tidyr unite
#' @importFrom jsonlite toJSON
#'
#' @return A list contains the estimated parameters and the results of execution
#' detection.
#' @export
#'
scGAN_estimation <- function(
    ref_data,
    other_prior,
    verbose = FALSE,
    seed
){

  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  ### (1. Check docker installation
  docker_install <- dynwrap::test_docker_installation()
  if(!docker_install){
    stop("Docker has not been installed or started! Please check it!")
  }
  ### (2. Check the installation of simpipe docker image
  images <- babelwhale::list_docker_images() %>%
    tidyr::unite("Repository", "Tag", sep = ":", col = "Image") %>%
    dplyr::pull("Image")

  if(!"fhausmann/scgan:latest" %in% images){
    # If not, pull fhausmann/scgan:latest
    babelwhale::pull_container(container_id = "fhausmann/scgan:latest")
  }

  ## Local file to store .h5ad file
  local_path <- system.file("scGAN", package = "simmethods")
  tmp_path <- tempdir() %>% simutils::fix_path()
  ## Move /scGAN to a tmp dir
  file.copy(from = local_path,
            to = tmp_path,
            overwrite = TRUE,
            recursive = TRUE)
  ## Docker file path
  docker_path <- "/scGAN"
  ## change parameters in .json file
  parameters <- simutils::change_scGAN_parameters("use_scGAN",
                                                  "GPU" = 1,
                                                  "min_cells" = 0,
                                                  "min_genes" = 0,
                                                  "experiments_dir" = "/scGAN",
                                                  "raw_input" = "/scGAN/input_data.h5ad",
                                                  "max_steps" = 100000)
  ## Save to /scGAN
  param_json <- jsonlite::toJSON(parameters,
                                 pretty = 4,
                                 dataframe = "values",
                                 null = "null",
                                 auto_unbox = TRUE)
  write(param_json, file = file.path(tmp_path, "scGAN", "parameters.json"))

  ## convert data to .h5ad and save to tmp_path
  new_data <- simutils::scgan_data_conversion(data = data,
                                              data_id = "input_data",
                                              save_to_path = file.path(tmp_path, "scGAN"),
                                              verbose = verbose)
  # Prepare the input parameters-----------------------------------------------
  ## 1. docker image working directory
  wd <- c("--workdir", "/scGAN")
  ## 2. docker image directory of the mount point
  docker_path <- "/scGAN"
  ## 3. volumes
  volumes <- paste0(file.path(tmp_path, "scGAN"), ":", docker_path)
  volumes_command <- c("-v", volumes)
  ## 4. verbose
  verbose <- verbose
  ## 5. args
  args <- c("python", "main.py", "--param", "parameters.json","--process")
  ## container name
  name <- simutils::time_string()
  container_name <- c("--name", name)
  ## 7. container id
  container_id <- "fhausmann/scgan"
  ## 8. docker command
  config <- babelwhale::get_default_config()
  processx_command <- Sys.which(config$backend) %>% unname()
  ## 9. processx arguments
  processx_args <- c("run",
                     container_name,
                     wd,
                     volumes_command,
                     c("--gpus", "all"),
                     container_id,
                     args)
  ## run
  process <- processx::run(
    command = processx_command,
    args = processx_args,
    echo_cmd = verbose,
    echo = verbose,
    spinner = TRUE,
    error_on_status = FALSE,
    cleanup_tree = TRUE
  )

  ##############################################################################
  ####                            Estimation                                 ###
  ##############################################################################
  args_train <- c("python", "main.py", "--param", "parameters.json", "--train")
  ## train container name
  train_name <- simutils::time_string()
  container_train_name <- c("--name", train_name)
  ## processx args
  processx_args_train <- c("run",
                           container_train_name,
                           wd,
                           volumes_command,
                           c("--gpus", "all"),
                           container_id,
                           args_train)

  if(verbose){
    message("Estimating parameters using scGAN, This step may be long...")
  }
  # Seed
  set.seed(seed)
  # Estimation
  tryCatch({
    # Estimate parameters from real data and return parameters and detection results
    estimate_detection <- peakRAM::peakRAM(
      estimate_result <- process <- processx::run(
        command = processx_command,
        args = processx_args_train,
        echo_cmd = verbose,
        echo = verbose,
        spinner = TRUE,
        error_on_status = FALSE,
        cleanup_tree = TRUE))
  }, error = function(e){
    as.character(e)
  })
  estimate_result <- list(newick_tree = estimate_result,
                          data_dim = dim(ref_data))
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  estimate_output <- list(estimate_result = estimate_result,
                          estimate_detection = estimate_detection)
  return(estimate_output)
}


#' Simulate Datasets by PROSSTT
#'
#' @param parameters A object generated by \code{\link[simutils]{make_trees}}
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
#' @importFrom glue glue
#' @importFrom stringr str_split str_count str_extract_all str_replace
#' @importFrom reticulate source_python
#'
#' @export
#'
# PROSSTT_simulation <- function(parameters,
#                                other_prior,
#                                return_format,
#                                verbose = FALSE,
#                                seed
# ){
#   ##############################################################################
#   ####                            Environment                                ###
#   ##############################################################################
#   if(!requireNamespace("simutils", quietly = TRUE)){
#     cat("Splatter is not installed on your device\n")
#     cat("Installing simutils...\n")
#     devtools::install_github("duohongrui/simutils")
#   }
#   ##############################################################################
#   ####                               Check                                   ###
#   ##############################################################################
#   # nGenes
#   if(!is.null(other_prior[["nGenes"]])){
#     gene_num <- other_prior[["nGenes"]]
#   }else{
#     gene_num <- parameters[["data_dim"]][1]
#   }
#   if(is.null(other_prior[["newick_tree"]])){
#     data_dim <- parameters[["data_dim"]]
#     # Return to users
#     cat(glue::glue("nCells: {data_dim[2]}"), "\n")
#     cat(glue::glue("nGenes: {gene_num}"), "\n")
#     newick_tree <- parameters[["newick_tree"]]
#     inter_cell <- stringr::str_split(newick_tree, pattern = "[)]", simplify = TRUE)
#     len_inter <- length(inter_cell)-2
#     change_index <- which(as.logical(stringr::str_count(inter_cell, pattern = "^:")))
#     inter_cell[change_index] <- paste0(LETTERS[1:len_inter], inter_cell[change_index])
#     newick_tree <- paste(inter_cell, collapse = ")")
#     edge <- unlist(stringr::str_extract_all(string = newick_tree, pattern = ":[:digit:]+[.]*[:digit:]*"))
#     node <- length(edge)+1
#     ncells <- data_dim[2]
#     cell_allo <- c(rep(round(ncells/node), node-1), ncells-(round(ncells/node)*(node-1)))
#     for(o in 1:length(edge)){
#       newick_tree <- stringr::str_replace(newick_tree,
#                                           pattern = edge[o],
#                                           replacement = paste0(":", as.character(cell_allo[o])))
#     }
#     newick_tree <- paste0(stringr::str_split(newick_tree, pattern = ";",simplify = T),
#                           paste0('Z:', cell_allo[node], ";"))[1]
#   }else{
#     newick_tree <- other_prior[["newick_tree"]]
#   }
#   # alpha
#   if(!is.null(other_prior[["alpha"]])){
#     alpha <- other_prior[["alpha"]]
#   }else{
#     alpha <- 0.2
#   }
#   # beta
#   if(!is.null(other_prior[["beta"]])){
#     beta <- other_prior[["beta"]]
#   }else{
#     beta <- 3
#   }
#   # modules
#   if(!is.null(other_prior[["modules"]])){
#     modules <- other_prior[["modules"]]
#   }else{
#     modules <- 15
#   }
#   exec_text <- system.file("python", "PROSSTT_python.py", package = "simmethods")
#
#   reticulate::source_python(exec_text)
#
#   simulation_params <- list(newick_string = newick_tree,
#                             modules = modules,
#                             gene_num = gene_num,
#                             seed = as.integer(seed),
#                             alpha = alpha,
#                             beta = beta)
#   ##############################################################################
#   ####                            Simulation                                 ###
#   ##############################################################################
#   if(verbose){
#     cat("Simulating datasets using PROSSTT\n")
#   }
#   # Estimation
#   tryCatch({
#     simulate_detection <- peakRAM::peakRAM(
#       simulate_result <- do.call("PROSSTT_sim_Python", simulation_params))
#   }, error = function(e){
#     as.character(e)
#   })
#   ##############################################################################
#   ####                        Format Conversion                              ###
#   ##############################################################################
#   colnames(simulate_result) <- paste0("Cell", 1:ncol(simulate_result))
#   rownames(simulate_result) <- paste0("Gene", 1:nrow(simulate_result))
#   ## col_data
#   col_data <- data.frame("cell_name" = colnames(simulate_result))
#   ## row_data
#   row_data <- data.frame("gene_name" = rownames(simulate_result))
#   # Establish SingleCellExperiment
#   simulate_result <- SingleCellExperiment::SingleCellExperiment(list(counts = simulate_result),
#                                                                 colData = col_data,
#                                                                 rowData = row_data)
#
#   simulate_result <- simutils::data_conversion(SCE_object = simulate_result,
#                                                return_format = return_format)
#
#   ##############################################################################
#   ####                           Ouput                                       ###
#   ##############################################################################
#   simulate_output <- list(simulate_result = simulate_result,
#                           simulate_detection = simulate_detection)
#   return(simulate_output)
# }

