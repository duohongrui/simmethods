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
#' @importFrom dynwrap test_docker_installation
#' @importFrom simutils fix_path change_scGAN_parameters data_conversion scgan_data_conversion time_string
#' @importFrom babelwhale list_docker_images pull_container get_default_config
#' @importFrom tidyr unite
#' @importFrom jsonlite toJSON
#' @importFrom processx run
#' @details
#' scGAN is a novel method to simulate single-cell RNA-seq datasets using generative adversarial neural networks and users can only execute it via docker images. `scGAN_estimation` and `scGAN_simulation` functions have already implemented the codes that users can use scGAN in R environment.
#' There are some notes that users should know:
#' 1. Please install docker on you device or remote service.
#' 2. The estimation step may take a long time as scGAN trains data reference data via neural networks.
#' 3. The result of estimation will be returned as a file path which is the mounting point to connect the path in docker containers. Users can go to the mounting point to see the training result.
#'
#' There are some parameters that users may often set:
#' 1. `group.condition`. Users can input cell group information of numeric vectors. If not, clustering will be performed before the estimation step.
#' 2. `max_steps`. The max training step to train the reference data. Default is 1000000.
#' 3. `GPU`. How many GPU cores to use when training the data. This can be set as `all`. Default is 1.
#' 4. `min_cells`. Include features detected in at least this many cells when preprocessing.
#' 5. `min_genes`. Include cells where at least this many features are detected when preprocessing.
#' 6. `res`. The clustering resolution. Default is 0.15.
#' @return A list contains the estimated parameters and the results of execution
#' detection.
#' @export
#' @references
#' Marouf M, Machart P, Bansal V, et al. Realistic in silico generation and augmentation of single-cell RNA-seq data using generative adversarial networks. Nature communications, 2020, 11(1): 1-12.
scGAN_estimation <- function(
    ref_data,
    other_prior = NULL,
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
  ### (2. Check the installation of scgan docker image
  if(!requireNamespace("babelwhale", quietly = TRUE)){
    install.packages("babelwhale")
  }
  images <- babelwhale::list_docker_images() %>%
    tidyr::unite("Repository", "Tag", sep = ":", col = "Image") %>%
    dplyr::pull("Image")
  if(!"fhausmann/scgan:latest" %in% images){
    # If not, pull fhausmann/scgan:latest
    babelwhale::pull_container(container_id = "fhausmann/scgan:latest")
  }

  ## Local file to store rds file
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
  if(is.null(other_prior[["res"]])){
    res <- 0.15
  }else{
    res <- other_prior[["res"]]
  }
  if(is.null(other_prior[["GPU"]])){
    other_prior[["GPU"]] <- 1
  }
  other_prior[["experiments_dir"]] <- "/scGAN"
  other_prior[["raw_input"]] <- "/scGAN/input_data.h5ad"
  parameters <- simutils::change_scGAN_parameters("use_scGAN",
                                                  other_prior)
  ## group information
  if(is.null(other_prior[["group.condition"]])){
    group <- NULL
  }else{
    group <- other_prior[["group.condition"]]-1
  }

  ## Save to /scGAN
  if(!requireNamespace("jsonlite", quietly = TRUE)){
    install.packages("jsonlite")
  }
  param_json <- jsonlite::toJSON(parameters,
                                 pretty = 4,
                                 dataframe = "values",
                                 null = "null",
                                 auto_unbox = TRUE)
  write(param_json, file = file.path(tmp_path, "scGAN", "parameters.json"))

  ## convert data to .h5ad and save to tmp_path
  new_data <- simutils::scgan_data_conversion(data = ref_data,
                                              data_id = "input_data",
                                              res = res,
                                              group = group,
                                              save_to_path = file.path(tmp_path, "scGAN"),
                                              verbose = verbose)
  cluster_number <- table(new_data$cluster_number)
  if(other_prior[["min_cells"]] == 0){
    need_add <- ncol(ref_data)-sum(cluster_number)
    assign_result <- simutils::proportionate(number = need_add,
                                             result_sum_strict = need_add,
                                             prop = rep(1/length(cluster_number),
                                                        length(cluster_number)),
                                             prop_sum_strict = 1,
                                             digits = 0)
    cluster_number <- cluster_number + assign_result
  }

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
  args <- c("python",
            "main.py",
            "--param",
            "parameters.json",
            "--process")
  ## container name
  name <- simutils::time_string()
  container_name <- c("--name", name)
  ## 7. container id
  container_id <- "fhausmann/scgan"
  ## 8. processx arguments
  processx_args <- c("run",
                     container_name,
                     wd,
                     volumes_command,
                     c("--gpus", "all"),
                     container_id,
                     args)
  ## run preprocessing
  if(!requireNamespace("processx", quietly = TRUE)){
    install.packages("processx")
  }
  process <- processx::run(
    command = "docker",
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
    message("Estimating parameters using scGAN, this step may take a long time...")
  }
  # Seed
  set.seed(seed)
  # Estimation
  tryCatch({
    # Estimate parameters from real data and return parameters and detection results
    estimate_detection <- peakRAM::peakRAM(
      estimate_result <- process <- processx::run(
        command = "docker",
        args = processx_args_train,
        echo_cmd = verbose,
        echo = verbose,
        spinner = TRUE,
        error_on_status = FALSE,
        cleanup_tree = TRUE))
  }, error = function(e){
    as.character(e)
  })
  estimate_result <- list(local_path = file.path(tmp_path, "scGAN"),
                          cluster_number = cluster_number,
                          gpu = other_prior[["GPU"]])
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  estimate_output <- list(estimate_result = estimate_result,
                          estimate_detection = estimate_detection)
  return(estimate_output)
}



#' Simulate Datasets by scGAN
#'
#' This function is used to simulate datasets from learned parameters in the docker container.
#'
#' @param parameters A object generated by [simmethods::scGAN_estimation()]
#' @param return_format A character. Alternatives choices: list, SingleCellExperiment,
#' Seurat, h5ad. If you select `h5ad`, you will get a path where the .h5ad file saves to.
#' @param verbose Logical. Whether to return messages or not.
#' @param seed A random seed.
#' @importFrom glue glue
#' @importFrom stringr str_split str_count str_extract_all str_replace
#' @importFrom reticulate source_python
#' @export
#' @references
#' Marouf M, Machart P, Bansal V, et al. Realistic in silico generation and augmentation of single-cell RNA-seq data using generative adversarial networks. Nature communications, 2020, 11(1): 1-12.
scGAN_simulation <- function(
    parameters,
    return_format,
    verbose = FALSE,
    seed
){
  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  # Prepare the input parameters-----------------------------------------------
  ## 1. docker image working directory
  wd <- c("--workdir", "/scGAN")
  ## 2. docker image directory of the mount point
  docker_path <- "/scGAN"
  ## 3. volumes
  local_path <- parameters$local_path
  volumes <- paste0(local_path, ":", docker_path)
  volumes_command <- c("-v", volumes)
  ## 4. verbose
  verbose <- verbose
  ## 5. args
  cell_number <- parameters[["cluster_number"]]
  args <- c("python",
            "main.py",
            "--param",
            "parameters.json",
            "--generate",
            "--cells_no",
            as.numeric(cell_number),
            "--model_path",
            "/scGAN/use_scGAN/",
            "--save_path",
            "simulation_result.h5ad")
  ## container name
  name <- simutils::time_string()
  container_name <- c("--name", name)
  ## 7. container id
  container_id <- "fhausmann/scgan"
  ## 8. processx arguments
  processx_args <- c("run",
                     container_name,
                     wd,
                     volumes_command,
                     c("--gpus", parameters[["gpu"]]),
                     container_id,
                     args)
  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    message("Simulating datasets using scGAN")
  }
  # Estimation
  tryCatch({
    simulate_detection <- peakRAM::peakRAM(
      simulate_result <- processx::run(
        command = "docker",
        args = processx_args,
        echo_cmd = verbose,
        echo = verbose,
        spinner = TRUE,
        error_on_status = FALSE,
        cleanup_tree = TRUE
      ))
  }, error = function(e){
    as.character(e)
  })
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  simulation_data_path <- file.path(local_path, "simulation_result.h5ad")
  if(!requireNamespace("sceasy")){
    message("Installing sceasy...")
    devtools::install_github("cellgeni/sceasy")
  }
  if(verbose){
    message("Read h5ad file into SeuratObject in R...")
  }
  simulation_result <- sceasy::convertFormat(simulation_data_path,
                                             from="anndata",
                                             to="seurat")
  counts <- as.matrix(simulation_result@assays$RNA@counts)
  colnames(counts) <- paste0("Cell", 1:ncol(counts))
  rownames(counts) <- paste0("Gene", 1:nrow(counts))
  ## col_data
  group <- paste0("Group", as.numeric(simulation_result$cluster)+1)
  col_data <- data.frame("cell_name" = colnames(counts),
                         "group" = group,
                         row.names = colnames(counts))
  ## row_data
  row_data <- data.frame("gene_name" = rownames(counts),
                         row.names = rownames(counts))
  # Establish SingleCellExperiment
  if(verbose){
    message("Converting data...")
  }
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
  ## delete the directory
  system2(command = "rm", args = c("-rf", local_path))
  return(simulate_output)
}

