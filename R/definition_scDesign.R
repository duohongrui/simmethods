#' Get Information of scDesign
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' scDesign_method_definition <- scDesign_method_definition()
#'
scDesign_method_definition <- function(...){

  scDesign_parameters <- parameter_sets(
    param_reference(
      id = "realcount",
      type = "matrix",
      default = NULL,
      process = "simulation",
      force = TRUE,
      description = "A numeric matrix with rows representing genes and columns representing cells. Gene names are given as row names.",
      function_name = "design_data"
    ),
    param_integer(
      id = "S",
      default = 1e+08L,
      process = "simulation",
      description = "A number specifying the total number of RNA-seq reads. Default to 1e8. When ngroup > 1, S should be a vector specifying the read number in each cell state.",
      function_name = "design_data"
    ),
    param_others(
      id = "ncell",
      type = c("integer", "vector"),
      force = TRUE,
      process = "simulation",
      description = "An integer specifying the number of cells. When ngroup > 1, ncell is vector specifying the number of cells in each cell state.",
      function_name = "design_data"
    ),
    param_integer(
      id = "ngroup",
      default = 1L,
      process = "simulation",
      description = "An integer giving the number of cell states to simulate. Defaults to 1.",
      function_name = "design_data"
    ),
    param_numeric(
      id = "pUp",
      process = "simulation",
      default = 0.05,
      lower = 0,
      upper = 1,
      description = "A value between 0 and 1 specifying the proportion of up regulated genes between two adjacent cell states. Defaults to 0.05 and only used when ngroup > 1.",
      function_name = "design_data"
    ),
    param_numeric(
      id = "pDown",
      process = "simulation",
      default = 0.05,
      lower = 0,
      upper = 1,
      description = "A value between 0 and 1 specifying the proportion of down regulated genes between two adjacent cell states. Defaults to 0.05 and only used when ngroup > 1.",
      function_name = "design_data"
    ),
    param_numeric(
      id = "fU",
      default = 5,
      process = "simulation",
      description = "A value specifying the upper bound of fold changes of differentially expressed genes. Deaults to 5.",
      function_name = "design_data"
    ),
    param_numeric(
      id = "fL",
      default = 1.5,
      process = "simulation",
      description = "A value specifying the upper bound of fold changes of differentially expressed genes. Deaults to 5.",
      function_name = "design_data"
    ),
    param_integer(
      id = "ncores",
      default = 1L,
      process = "simulation",
      description = "An integer specifying the number of cores used for parallel computation. Defaults to 1.",
      function_name = "design_data"
    ),
    param_vector(
      id = "exprmean",
      default = NULL,
      description = "A named vector of user-specified gene mean expression parameters. The names of exprmean should match the rownames of realcount in the same order. Supplied gene expression mean should be on the $log_10$ scale.",
      process = "simulation",
      function_name = "design_data"
    )
  )

  scDesign_method <- method_definition(
    method = "scDesign",
    programming = "R",
    url = "https://github.com/Vivianstats/scDesign",
    authors = authors_definition(
      first = "Wei Li",
      last = "Li",
      email = "vivian.li@rutgers.edu",
      github = "https://github.com/Vivianstats/scDesign",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "A statistical simulator scDesign for rational scRNA-seq experimental design",
      doi = "10.1093/bioinformatics/btz321",
      journal = "Bioinformatics",
      date = "2019",
      peer_review = TRUE
    ),
    description = "A statistical simulator for rational scRNA-seq experimental design.",
    vignette = "http://47.254.148.113/software/Simsite/references/methods/11-scdesign/")

  list(scDesign_method = scDesign_method,
       scDesign_parameters = scDesign_parameters)
}
