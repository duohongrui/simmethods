#' @importFrom dynwrap is_wrapper_with_expression
#' @importFrom dplyr full_join across case_when pull
#'
pseudotime_info <- function(ref_data, other_prior, col_data, seed){
  if(!requireNamespace("parallelDist", quietly = TRUE)){
    utils::install.packages("parallelDist")
  }
  if(!requireNamespace("dyndimred", quietly = TRUE)){
    devtools::install_github("dynverse/dyndimred")
  }
  ### dimension reduction
  dimred <- dyndimred::dimred_umap(t(ref_data))
  ### start cell
  if(is.null(other_prior[["start_cell_id"]])){
    set.seed(seed)
    start_cell_id <- colnames(ref_data)[sample(1:ncol(ref_data), size = 1)]
  }else{
    start_cell_id <- other_prior[["start_cell_id"]]
  }
  ### dynwrap_expression data
  if(dynwrap::is_wrapper_with_expression(other_prior[["dynwrap_data"]])){
    dynwrap_data <- other_prior[["dynwrap_data"]]
    traj_type <- dynwrap_data$trajectory_type
    milestone_network <- dynwrap_data$milestone_network
    start_milestone <- dynwrap_data$root_milestone_id
    start_cell <- names(dynwrap_data$grouping)[which(dynwrap_data$grouping %in% start_milestone)]
    set.seed(seed)
    start_cell_id <- start_cell[sample(1:length(start_cell), size = 1)]
    if(traj_type == "bifurcation"){
      inter_name <- BiocGenerics::unique(milestone_network$to)[1]
      li1_name <- c(start_milestone, milestone_network$from[-BiocGenerics::grep(start_milestone,
                                                                                milestone_network$from)][1],
                    inter_name)
      li2_name <- c(start_milestone, milestone_network$from[-BiocGenerics::grep(start_milestone,
                                                                                milestone_network$from)][2],
                    inter_name)
      dist <- parallelDist::parallelDist(dimred) %>% as.matrix()
      dist1 <- parallelDist::parallelDist(dimred[which(dynwrap_data$grouping %in% li1_name), ]) %>% as.matrix()
      dist2 <- parallelDist::parallelDist(dimred[which(dynwrap_data$grouping %in% li2_name), ]) %>% as.matrix()
      pseudotime1 <- data.frame("cell_name" = names(dist[which(dynwrap_data$grouping %in% li1_name), start_cell_id]),
                                "pseudotime1" = S4Vectors::unname(dist[which(dynwrap_data$grouping %in% li1_name), start_cell_id]))
      pseudotime2 <- data.frame("cell_name" = names(dist[which(dynwrap_data$grouping %in% li2_name), start_cell_id]),
                                "pseudotime2" = S4Vectors::unname(dist[which(dynwrap_data$grouping %in% li2_name), start_cell_id]))
      col_data <- col_data %>%
        dplyr::full_join(pseudotime1, by = "cell_name") %>%
        dplyr::full_join(pseudotime2, by = "cell_name") %>%
        mutate(
          dplyr::across(c("pseudotime1", "pseudotime2"), ~ tidyr::replace_na(.x, -1)),
          l1 = dplyr::case_when(
            pseudotime1 >= 0 & pseudotime2 < 0 ~ TRUE,
            pseudotime1 >= 0 & pseudotime2 >= 0 ~ TRUE,
            TRUE ~ FALSE
          ),
          l2 = dplyr::case_when(
            pseudotime1 < 0 & pseudotime2 >= 0 ~ TRUE,
            pseudotime1 >= 0 & pseudotime2 >= 0 ~ TRUE,
            TRUE ~ FALSE
          )
        )
      col_data$l1 <- factor(col_data$l1)
      col_data$l2 <- factor(col_data$l2)
      mu_formula = "s(pseudotime1, k = 10, by = l1, bs = 'cr') + s(pseudotime2, k = 10, by = l2, bs = 'cr')"
      pseudotime <- c("pseudotime1", "pseudotime2", "l1", "l2")
    }
    if(traj_type == "linear" | traj_type == 'cycle'){
      dist <- parallelDist::parallelDist(dimred) %>% as.matrix()
      pseudotime <- data.frame("cell_name" = names(dist[, start_cell_id]),
                               "pseudotime" = unname(dist[, start_cell_id]))
      col_data <- col_data %>%
        dplyr::full_join(pseudotime, by = "cell_name")
      col_data$l <- factor(col_data$l)
      mu_formula = "s(pseudotime, k = 10, bs = 'cr')"
      pseudotime <- c("pseudotime")
    }
    if(traj_type == "multifurcation" | traj_type == "tree"){
      dist <- parallelDist::parallelDist(dimred) %>% as.matrix()
      pseudotime <- data.frame("cell_name" = colnames(dist))

      for(i in 1:nrow(milestone_network)){
        from_cell <- milestone_network[i, ] %>% dplyr::pull("from")
        to_cell <- milestone_network[i, ] %>% dplyr::pull("to")
        pseudotime <- data.frame("cell_name" = names(dist[grep(from_cell, dynwrap_data$grouping), grep(start_cell_id, colnames(dist))]),
                                 "pseudotime" = S4Vectors::unname(dist[grep(from_cell, dynwrap_data$grouping), grep(start_cell_id, colnames(dist))]))
        col_data <- col_data %>%
          dplyr::full_join(pseudotime, by = "cell_name")
      }
      col_data <- col_data %>%
        mutate(
          dplyr::across(3:ncol(.), ~ tidyr::replace_na(.x, -1))
        )
      colnames(col_data)[3:ncol(col_data)] <- paste0("pseudotime", 1:nrow(milestone_network))
      pseudotime <- paste0("pseudotime", 1:nrow(milestone_network))
      mu_formula <- paste0(paste0("s(", pseudotime, ", k = 10, bs = 'cr')"), collapse = " + ")
    }
  }else{
    dist <- parallelDist::parallelDist(dimred) %>% as.matrix()
    pseudotime <- data.frame("cell_name" = names(dist[, start_cell_id]),
                             "pseudotime" = unname(dist[, start_cell_id]))
    col_data <- col_data %>%
      dplyr::full_join(pseudotime, by = "cell_name") %>%
      mutate(
        l = TRUE
      )
    col_data$l <- factor(col_data$l)
    mu_formula = "s(pseudotime, k = 10, bs = 'cr')"
    pseudotime <- c("pseudotime")
  }

  return(
    dplyr::lst(col_data, mu_formula, pseudotime)
  )
}



#' @importFrom scMultiSim sim_true_counts add_expr_noise
excution_function <- function(options, seed){
  simulate_result <- scMultiSim::sim_true_counts(options)
  scMultiSim::add_expr_noise(simulate_result, randseed = seed)
  return(simulate_result)
}


#' @importFrom scMultiSim divide_batches
excution_batch_function <- function(options, seed, nbatch){
  simulate_result <- scMultiSim::sim_true_counts(options)
  scMultiSim::add_expr_noise(simulate_result, randseed = seed)
  scMultiSim::divide_batches(simulate_result, nbatch = nbatch, randseed = seed)
  return(simulate_result)
}
