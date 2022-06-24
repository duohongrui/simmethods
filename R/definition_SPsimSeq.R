#' Get Information of SPsimSeq
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' SPsimSeq_method_definition <- SPsimSeq_method_definition()
#'
SPsimSeq_method_definition <- function(...){

  SPsimSeq_parameters <- parameter_sets(
    param_reference(
      id = "s.data",
      type = "matrix",
      default = NULL,
      force = TRUE,
      description = "A real source dataset (a SingleCellExperiment class or a matrix/data.frame of counts with genes in rows and samples in columns)",
      function_name = "SPsimSeq"
    ),
    param_integer(
      id = "n.sim",
      default = 1L,
      description = "An integer for the number of simulations to be generated",
      function_name = "SPsimSeq"
    ),
    param_vector(
      id = "batch",
      default = "rep(1, ncol(s.data))",
      description = "NULL or a vector containing batch indicators for each sample/cell in the source data.",
      function_name = "SPsimSeq"
    ),
    param_vector(
      id = "group",
      default = "rep(1, ncol(s.data))",
      description = "NULL or a vector containing group indicators for each sample/cell in the source data.",
      function_name = "SPsimSeq"
    ),
    param_integer(
      id = "n.genes",
      default = 1000L,
      description = "A numeric value for the total number of genes to be simulated",
      function_name = "SPsimSeq"
    ),
    param_vector(
      id = "batch.config",
      default = 1,
      description = "A numerical vector containing fractions for the composition of samples/cells per batch. The fractions must sum to 1. The number of batches to be simulated is equal to the size of the vector. (Example, batch.config=c(0.6, 0.4) means simulate 2 batches with 60% of the simulated samples/cells in batch 1 and the rest 40% in the second batch. Another example, batch.config=c(0.3, 0.35, 0.25) means simulate 3 batches with the first, second and third batches contain 30%, 35% and 25% of the samples/cells, respectively).",
      function_name = "SPsimSeq"
    ),
    param_vector(
      id = "group.config",
      default = 1,
      description = "A numerical vector containing fractions for the composition of samples/cells per group. Similar definition to 'batch.config'. The number of groups to be simulated is equal to the size of the vector. The fractions must sum to 1.",
      function_name = "SPsimSeq"
    ),
    param_numeric(
      id = "pDE",
      default = 0.1,
      description = "A numeric value between 0 and 1 indicating the desired fraction of DE genes in the simulated data.",
      lower = 0,
      upper = 1,
      function_name = "SPsimSeq"
    ),
    param_others(
      id = "cand.DE.genes",
      type = "list",
      default = NULL,
      description = "A list object contatining canidiate null and non-null (DE/predictor) genes. If NULL (the default), an internal function determines candidate genes based on log-fold-change and other statistics. The user can also pass a list of canidate null and non-null genes (they must be disjoint). In particular, the list should contain two character vectors (for the name of the features/genes in the source data) with names 'null.genes' and 'nonnull.genes'. For example, cand.DE.genes=list(null.genes=c('A', 'B'), nonnull.genes=c('C', 'D')).",
      function_name = "SPsimSeq"
    ),
    param_numeric(
      id = "lfc.thrld",
      default = 0.5,
      lower = 0,
      border = FALSE,
      description = "A positive numeric value for the minimum absolute log-fold-change for selecting candidate DE genes in the source data (when group is not NULL, pDE>0 and cand.DE.genes is NULL).",
      function_name = "SPsimSeq"
    ),
    param_numeric(
      id = "t.thrld",
      default = 2.5,
      lower = 0,
      border = FALSE,
      description = "A positive numeric value for the minimum absolute t-test statistic for the log-fold-changes of genes for selecting candidate DE genes in the source data (when group is not NULL, pDE>0 and cand.DE.genes is NULL).",
      function_name = "SPsimSeq"
    ),
    param_numeric(
      id = "llStat.thrld",
      default = 5,
      lower = 0,
      border = FALSE,
      description = "A positive numeric value for the minimum squared test statistics from the log-linear model to select candidate DE genes in the source data (when group is not NULL, pDE>0 and cand.DE.genes is NULL).",
      function_name = "SPsimSeq"
    ),
    param_others(
      id = "tot.samples",
      type = "integer",
      default = "ncol(s.data)",
      description = "A numerical value for the total number of samples/cells to be simulated.",
      function_name = "SPsimSeq"
    ),
    param_Boolean(
      id = "model.zero.prob",
      default = FALSE,
      description = "A logical value whether to model the zero expression probability separately (suitable for simulating of single-cell RNA-seq data or zero-inflated data).",
      function_name = "SPsimSeq"
    ),
    param_Boolean(
      id = "genewiseCor",
      default = TRUE,
      description = "A logical value, if TRUE (default) the simulation will retain the gene-to-gene correlation structure of the source data using Gausian-copulas . Note that if it is TRUE, the program will be slow or it may fail for a limited memory size.",
      function_name = "SPsimSeq"
    ),
    param_Boolean(
      id = "log.CPM.transform",
      default = TRUE,
      description = "A logical value. If TRUE, the source data will be transformed into log-(CPM+const) before estimating the probability distributions.",
      function_name = "SPsimSeq"
    ),
    param_vector(
      id = "lib.size.params",
      description = "NULL or a named numerical vector containing parameters for simulating library sizes from log-normal distribution. If lib.size.params =NULL (default), then the package will fit a log-normal distribution for the library sizes in the source data to simulate new library sizes. If the user would like to specify the parameters of the log-normal distribution for the desired library sizes, then the log-mean and log-SD params of rlnorm() functions can be passed using this argument. Example, lib.size.params = c(meanlog=10, sdlog=0.2).",
      function_name = "SPsimSeq"
    ),
    param_Boolean(
      id = "variable.lib.size",
      default = FALSE,
      description = "A logical value. If FALSE (default), then the expected library sizes are simulated once and remains the same for every replication (if n.sim>1).",
      function_name = "SPsimSeq"
    ),
    param_others(
      id = "w",
      type = "unknown",
      description = "see? hist",
      function_name = "SPsimSeq"
    ),
    param_character(
      id = "result.format",
      default = "SCE",
      alternatives = c("SCE", "list"),
      description = "A character value for the type of format for the output. Choice can be 'SCE' for SingleCellExperiment class or 'list' for a list object that contains the simulated count, column information and row information.",
      function_name = "SPsimSeq"
    ),
    param_Boolean(
      id = "return.details",
      default = FALSE,
      description = "A logical value. If TRUE, detailed results including estimated parameters and densities will be returned.",
      function_name = "SPsimSeq"
    ),
    param_Boolean(
      id = "verbose",
      default = TRUE,
      description = "A logical value, if TRUE a message about the status of the simulation will be printed on the console.",
      function_name = "SPsimSeq"
    ),
    param_numeric(
      id = "prior.count",
      default = 1,
      lower = 0,
      border = FALSE,
      description = "A positive constant to be added to the CPM before log transformation, to avoid log(0). The default is 1.",
      function_name = "SPsimSeq"
    ),
    param_numeric(
      id = "const.mult",
      default = 1e+06,
      lower = 0,
      border = FALSE,
      description = "A constant by which the count are scaled. Usually 1e6 to get CPM.",
      function_name = "SPsimSeq"
    ),
    param_numeric(
      id = "n.mean.class",
      default = 0.2,
      lower = 0,
      border = FALSE,
      description = "A fraction of the number of genes for the number of groups to be created for the mean log CPM of genes.",
      function_name = "SPsimSeq"
    ),
    param_numeric(
      id = "minFracZeroes",
      default = 0.25,
      lower = 0,
      border = FALSE,
      description = "Minimum fraction of zeroes before a zero inflation model is fitted.",
      function_name = "SPsimSeq"
    )
  )

  SPsimSeq_method <- method_definition(
    method = "SPsimSeq",
    programming = "R",
    url = "https://www.bioconductor.org/packages/release/bioc/html/SPsimSeq.html",
    authors = authors_definition(
      first = "Alemu Takele",
      last = "Assefa",
      email = "alemutakele.assefa@ugent.be",
      github = "https://github.com/CenterForStatistics-UGent/SPsimSeq",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "SPsimSeq: semi-parametric simulation of bulk and single-cell RNA-sequencing data",
      doi = "10.1093/bioinformatics/btaa105",
      journal = "Bioinformatics",
      date = "2020",
      peer_review = TRUE
    ),
    description = "SPsimSeq uses a specially designed exponential family for density estimation to constructs the distribution of gene expression levels from a given real RNA sequencing data.")

  list(SPsimSeq_method = SPsimSeq_method,
       SPsimSeq_parameters = SPsimSeq_parameters)
}
