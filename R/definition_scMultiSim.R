#' Get Information of scMultiSim
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' scMultiSim_method_definition <- scMultiSim_method_definition()
#'
scMultiSim_method_definition <- function(...){

  scMultiSim_parameters <- parameter_sets(
    param_reference(
      id = "input_refernece",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      description = "Note that the input matrix for estimation is nessesary unless you can provide newick tree format string.",
      function_name = "simutils::make_trees"
    ),
    param_others(
      id = "num.cells",
      type = "integer",
      default = "ncol(data)",
      process = "simulation",
      description = "Specify the number of cells.",
      function_name = "scMultiSim"
    ),
    param_others(
      id = "num.genes",
      type = "integer",
      default = "nrow(data)",
      process = "simulation",
      description = "Specify the number of genes.",
      function_name = "scMultiSim"
    ),
    param_Boolean(
      id = "discrete.cif",
      default = FALSE,
      process = "simulation",
      description = "Whether generating discrete or continuous cell population.",
      function_name = "scMultiSim"
    ),
    param_integer(
      id = "rand.seed",
      default = 0L,
      process = "simulation",
      description = "scMultiSim should produce the same result if all other parameters are the same.",
      function_name = "scMultiSim"
    ),
    param_integer(
      id = "threads",
      default = 1L,
      process = "simulation",
      description = "Use multithreading only when generating the CIF matrix. It will not speed up the simulation a lot, thus not recommended.",
      function_name = "scMultiSim"
    ),
    param_numeric(
      id = "unregulated.gene.ratio",
      process = "simulation",
      default = 0.1,
      description = "Ratio of unreulated to regulated genes. When a GRN is supplied with N genes, scMultiSim will simulate N * r_u extra (unregulated) genes.",
      function_name = "scMultiSim"
    ),
    param_others(
      id = "tree",
      type = "phylo",
      default = NULL,
      process = "simulation",
      force = TRUE,
      description = "The cell differential tree, which will be used to generate cell trajectories (if discrete.cif = T) or clusters (if discrete.cif = F). In discrete population mode, only the tree tips will be used. Three demo trees, Phyla5(), Phyla3() and Phyla1(), are provided.",
      function_name = "scMultiSim"
    ),
    param_integer(
      id = "num.cifs",
      default = 50L,
      process = "simulation",
      description = "Total number of differential and non-differential CIFs, which can be viewed as latent representation of cells.",
      function_name = "scMultiSim"
    ),
    param_numeric(
      id = "diff.cif.fraction",
      default = 0.9,
      process = "simulation",
      description = "Fraction of differential CIFs. Differential CIFs encode the cell type information, while non-differential CIFs are randomly sampled for each cell.",
      function_name = "scMultiSim"
    ),
    param_numeric(
      id = "cif.sigma",
      default = 0.1,
      process = "simulation",
      description = "The distribution used to sample CIF values.",
      function_name = "scMultiSim"
    ),
    param_numeric(
      id = "cif.center",
      default = 1,
      process = "simulation",
      description = "The distribution used to sample CIF values.",
      function_name = "scMultiSim"
    ),
    param_Boolean(
      id = "use.impulse",
      default = FALSE,
      description = "In continuous population mode, when sampling CIFs along the tree, use the impulse model rather than the default gaussian random walk.",
      process = "simulation",
      function_name = "scMultiSim"
    )
  )

  scMultiSim_method <- method_definition(
    method = "scMultiSim",
    programming = "R",
    url = "https://github.com/ZhangLabGT/scMultiSim",
    authors = authors_definition(
      first = "Hechen",
      last = "Li",
      email = NULL,
      github = "https://github.com/ZhangLabGT/scMultiSim",
      orcid = "0000-0002-1455-0041"
    ),
    manuscript = manuscript_definition(
      title = "scMultiSim: simulation of multi-modality single cell data guided by cell-cell interactions and gene regulatory networks",
      doi = "10.1101/2022.10.15.512320",
      journal = "bioRxiv",
      date = "2022",
      peer_review = FALSE
    ),
    description = "scMultiSim is an in silico simulator that generates multi-modality data of single-cells, including gene expression, chromatin accessibility, RNA velocity, and spatial location of cells. It takes a cell differential tree and a gene regulatory network (GRN) as input, and simulates spliced and unspliced counts while accounting for the relationships between modalities.",
    vignette = "http://47.254.148.113/software/Simsite/references/methods/31-scmultisim/")

  list(scMultiSim_method = scMultiSim_method,
       scMultiSim_parameters = scMultiSim_parameters)
}
