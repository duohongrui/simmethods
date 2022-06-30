#' Get Information of hierarchicell
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' hierarchicell_method_definition <- hierarchicell_method_definition()
#'
hierarchicell_method_definition <- function(...){

  hierarchicell_parameters <- parameter_sets(
    param_reference(
      id = "expr",
      type = "dataframe",
      default = NULL,
      process = "estimation",
      force = TRUE,
      description = "A data.frame that has been output by filter_counts where the unique cell identifier is in column one and the sample identifier is in column two with the remaining columns all being genes.",
      function_name = "compute_data_summaries"
    ),
    param_character(
      id = "type",
      default = "Raw",
      process = "estimation",
      alternative = c("Raw", "PerMillion", "Norm"),
      description = "an identifier for the type of data being submitted. If it is raw counts put 'Raw', if it is TPM or some other normalized counts per million then type 'PerMillion' or 'Norm'. The program assumes data is in one of these two formats. Other normalizations (i.e., logs) that have negative values will cause the program to malfunction.",
      function_name = "compute_data_summaries"
    ),
    param_others(
      id = "data_summaries",
      type = "SingleCellExperimrnt",
      description = "An R object that has been output by the package's compute_data_summaries function. No default.",
      function_name = "simulate_hierarchicell "
    ),
    param_others(
      id = "n_genes",
      type = "integer",
      default = "nrow(ref_data)",
      description = "An integer. The number of genes you would like to simulate for your dataset. Too large of a number may cause memory failure and may slow the simulation down tremendously. We recommend an integer less than 100,000. ",
      function_name = "simulate_hierarchicell"
    ),
    param_integer(
      id = "n_per_group",
      default = 1L,
      description = "An integer. The number of independent samples per case/control group for simulation. Use when 'binary' is specified as the outcome. Creates a balanced design, for unbalanced designs, specify n_cases and n_controls separately. If not specifying a foldchange, the number of cases and controls does not matter. ",
      function_name = "simulate_hierarchicell"
    ),
    param_integer(
      id = "n_cases",
      default = 1L,
      description = "An integer. The number of independent control samples for simulation. ",
      function_name = "simulate_hierarchicell"
    ),
    param_integer(
      id = "n_controls",
      default = 1L,
      description = "An integer. The number of independent control samples for simulation. ",
      function_name = "simulate_hierarchicell"
    ),
    param_others(
      id = "cells_per_control",
      type = "integer",
      default = "ncol(ref_data)/2",
      description = "An integer. The mean number of cells per control you would like to simulate. Too large of a number may cause memory failure and may slow the simulation down tremendously. We recommend an integer less than 1,000.",
      function_name = "simulate_hierarchicell"
    ),
    param_others(
      id = "cells_per_case",
      type = "integer",
      default = "ncol(ref_data)/2",
      description = "An integer. The mean number of cells per control you would like to simulate. Too large of a number may cause memory failure and may slow the simulation down tremendously. We recommend an integer less than 1,000.",
      function_name = "simulate_hierarchicell"
    ),
    param_character(
      id = "ncells_variation_type",
      default = "Poisson",
      alternative = c("Poisson", "NB", "Fixed"),
      description = "Either 'Poisson', 'NB', or 'Fixed'. Allows the number of cells per individual to be fixed at exactly the specified number of cells per individual, vary slightly with a poisson distribution with a lambda equal to the specified number of cells per individual, or a negative binomial with a mean equal to the specified number of cells and dispersion size equal to one.Defaults to 'Poisson'.",
      function_name = "simulate_hierarchicell"
    ),
    param_numeric(
      id = "foldchange",
      default = 2,
      lower = 0,
      upper = 10,
      description = "An integer between 1 and 10. The amount of fold change to simulate a difference in expression between case and control groups. The foldchange changes genes in either direction, so a foldchange of 2 would cause the mean expression in cases to be either twice the amount or half the amount for any particular gene. Defaults to 1.",
      function_name = "simulate_hierarchicell"
    ),
    param_numeric(
      id = "decrease_dropout",
      default = 0,
      lower = 0,
      upper = 1,
      description = "A numeric proportion between 0 and 1. The proportion by which you would like to simulate decreasing the amount of dropout in your data. For example, if you would like to simulate a decrease in the amount of dropout in your data by twenty percent, then 0.2 would be appropriate. This component of the simulation allows the user to adjust the proportion of dropout if they believe future experiments or runs will have improved calling rates (due to improved methods or improved cell viability) and thereby lower dropout rates. Defaults to 0.",
      function_name = "simulate_hierarchicell"
    ),
    param_numeric(
      id = "alter_dropout_cases",
      default = 0,
      lower = 0,
      upper = 1,
      description = "A numeric proportion between 0 and 1. The proportion by which you would like to simulate decreasing the amount of dropout between case control groups. For example, if you would like to simulate a decrease in the amount of dropout in your cases by twenty percent, then 0.2 would be appropriate. This component of the simulation allows the user to adjust the proportion of dropout if they believe the stochastic expression of a gene will differ between cases and controls. For a two-part hurdle model, like MAST implements, this will increase your ability to detect differences. Defaults to 0.",
      function_name = "simulate_hierarchicell"
    ),
    param_Boolean(
      id = "tSNE_plot",
      default = FALSE,
      description = "A TRUE/FALSE statement for the output of a tSNE plot to observe the global behavior of your simulated data. Seurat will need to be installed for this function to properly work. Defaults to FALSE.",
      function_name = "simulate_hierarchicell"
    )
  )

  hierarchicell_method <- method_definition(
    method = "hierarchicell",
    programming = "R",
    url = "https://github.com/kdzimm/hierarchicell",
    authors = authors_definition(
      first = "Kip D.",
      last = "Zimmerman",
      email = "kdzimmer@wakehealth.edu",
      github = "https://github.com/kdzimm/hierarchicell",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "Hierarchicell: an R-package for estimating power for tests of differential expression with single-cell data",
      doi = "10.1186/s12864-021-07635-w",
      journal = "BMC Genomics",
      date = "2021",
      peer_review = TRUE
    ),
    description = "An R package for simulating cell-type specific and hierarchical single-cell expression data.")

  list(hierarchicell_method = hierarchicell_method,
       hierarchicell_parameters = hierarchicell_parameters)
}
