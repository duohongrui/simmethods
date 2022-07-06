#' Get Information of PROSSTT
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' PROSSTT_method_definition <- PROSSTT_method_definition()
PROSSTT_method_definition <- function(...){

  PROSSTT_parameters <- parameter_sets(
    param_reference(
      id = "input_refernece",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      description = "Note that the input matrix for estimation is nessesary unless you can provide newick tree format string.",
      function_name = "simutils::make_trees"
    ),
    param_others(
      id = "newick_tree_format_string",
      type = "string",
      default = "(B:100,C:60)A:70;",
      description = "Usually it is learned by simutils::make_tree function. Note that the number in the string refers to the pseudotime which is similar to the length",
      function_name = "PROSSTT_python"
    ),
    param_integer(
      id = "modules",
      default = 15L,
      lower = 1L,
      description = "Total number of expression programs for the lineage tree.",
      function_name = "PROSSTT_python"
    ),
    param_others(
      id = "nGenes",
      type = "integer",
      default = "nrow(ref_data)",
      description = "Total number of genes.",
      function_name = "PROSSTT_python"
    ),
    param_others(
      id = "nCells",
      type = "integer",
      default = "ncol(ref_data)",
      description = "Total number of cells.",
      function_name = "PROSSTT_python"
    ),
    param_numeric(
      id = "alpha",
      default = 0.2,
      description = "Average alpha value.",
      function_name = "PROSSTT_python"
    ),
    param_numeric(
      id = "beta",
      default = 3,
      description = "Average beta value.",
      function_name = "PROSSTT_python"
    ),
    param_integer(
      id = "seed",
      default = 687680L,
      function_name = "PROSSTT_python"
    ))

  PROSSTT_method <- method_definition(
    method = "PROSSTT",
    programming = "Python",
    url = "http://wwwuser.gwdg.de/~compbiol/prosstt/doc/",
    authors = authors_definition(
      first = "Nikolaos",
      last = "Papadopoulos",
      email = NULL,
      github = "https://github.com/soedinglab/prosstt/",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "PROSSTT: probabilistic simulation of single-cell RNA-seq data for complex differentiation processes",
      doi = "10.1093/bioinformatics/btz078",
      journal = "Bioinformatics",
      date = "2019",
      peer_review = TRUE
    ),
    description = "PROSSTT is a package with code for the simulation of scRNAseq data for dynamic processes such as cell differentiation.")

  list(PROSSTT_method = PROSSTT_method,
       PROSSTT_parameters = PROSSTT_parameters)
}
