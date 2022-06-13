#' Get Information of powsimR
#'
#' @param ... ...
#'
#' @return A list contains the information of method and default parameters
#' @import simutils
#' @export
#'
#' @examples
#' powsimR_method_definition <- powsimR_method_definition()
#'
powsimR_method_definition <- function(...){

  powsimR_parameters <- parameter_sets(
    param_reference(
      id = "countData",
      type = c("matrix"),
      default = NULL,
      process = "estimation",
      description = "countData is a (UMI) count matrix of gene expression. Rows correspond to genes, columns to samples. The gene names should be given as rownames without '_' in the names. The samples names should be given as colnames. The count matrix should only contain the expression of one group, e.g. wildtype / untreated control / one cell type population.",
      function_name = "estimateParam"
    ),
    param_others(
      id = "readData",
      type = "matrix",
      process = "estimation",
      description = "readData is a the matching read count matrix of gene expression if countData is UMI and the same formatting should be applied. Default is NULL and users only need to supply the read count matrix if they plan to apply downsampling of UMI counts for simulations, see Setup.",
      function_name = "estimateParam"
    ),
    param_dataframe(
      id = "batchData",
      type = "data.frame",
      process = "estimation",
      description = "batchData is a data.frame for batch annotation. Rows correspond to samples. The first column should contain the batches, e.g. 'a', 'b', 'c', etc.",
      function_name = c("estimateParam", "estimateSpike")
    ),
    param_others(
      id = "spikeData",
      type = "matrix",
      default = NULL,
      process = "estimation",
      description = "spikeData is a count matrix. Rows correspond to spike-ins, columns to samples. The order of columns should be the same as in the countData. This is only needed for spike-in aware normalisation methods ('MR', 'Linnorm', 'scran', 'SCnorm', 'bayNorm', 'Census').",
      function_name = "estimateParam"
    ),
    param_others(
      id = "spikeInfo",
      type = "matrix",
      default = NULL,
      process = "estimation",
      description = "spikeInfo is a molecule count matrix of spike-ins. Rows correspond to spike-ins. The order of rows should be the same as in the spikeData. The column names should be 'SpikeID' and 'SpikeInput' for molecule counts of spike-ins. This is only needed for spike-in aware normalisation methods ().",
      function_name = c("estimateParam", "estimateSpike")
    ),
    param_vector(
      id = "Lengths",
      process = "estimation",
      description = "is a numeric vector of transcript lengths with the same length and order as the rows in countData. This variable is only needed for internal gene length corrections (TPM).",
      function_name = "estimateParam"
    ),
    param_vector(
      id = "MeanFragLengths",
      process = "estimation",
      description = "MeanFragLengths is a numeric vector of mean fragment lengths with the same length as columns in countData. This variable is only needed for internal gene length corrections (TPM).",
      function_name = c("estimateParam", "estimateSpike")
    ),
    param_character(
      id = "RNAseq",
      default = NULL,
      alternatives = c('bulk', 'singlecell'),
      force = TRUE,
      process = "estimation",
      description = "RNAseq is a character value: 'bulk' or 'singlecell'.",
      function_name = c("estimateParam", "estimateSpike")
    ),
    param_character(
      id = "Protocol",
      default = NULL,
      alternatives = c('UMI', 'Read'),
      force = TRUE,
      process = "estimation",
      description = "Protocol is a character value defining the type of counts given in countData. Options are 'UMI' (e.g. 10X Genomics, CEL-seq2) or 'Read' (e.g. Smart-seq2).",
      function_name = c("estimateParam", "estimateSpike")
    ),
    param_character(
      id = "Distribution",
      default = NULL,
      alternatives = c('NB', 'ZINB'),
      force = TRUE,
      process = "estimation",
      description = "Distribution is a character value: 'NB' for negative binomial or 'ZINB' for zero-inflated negative binomial distribution fitting.",
      function_name = "estimateParam"
    ),
    param_character(
      id = "Normalisation",
      default = NULL,
      alternatives = c("TMM", "MR", "PosCounts", "UQ",
                       "scran", "Linnorm", "sctransform",
                       "SCnorm", "Census", "depth", "none"),
      force = TRUE,
      process = "estimation",
      description = "Normalisation is a character value: 'TMM', 'MR', 'PosCounts', 'UQ', 'scran', 'Linnorm', 'SCnorm', 'Census', 'depth', 'none'. For more information, please consult the Details section.",
      function_name = "estimateParam"
    ),
    param_numeric(
      id = "GeneFilter",
      default = 0.05,
      lower = 0,
      upper = 1,
      border = FALSE,
      process = "estimation",
      description = "GeneFilter is a numeric vector indicating the minimal proportion of nonzero expression values for a gene across all samples to be considered expressed and used for normalisation and parameter estimation. The default is 0.05, i.e. at least 5% of the expression values per gene need to be nonzero.",
      function_name = "estimateParam"
    ),
    param_numeric(
      id = "SampleFilter",
      default = 5,
      lower = 0,
      border = FALSE,
      process = "estimation",
      description = "SampleFilter is a numeric vector indicating the minimal number of MADs (median absolute deviation) away from the median number of features detected as well as sequencing depth across all samples so that outlying samples are removed prior to normalisation and parameter estimation. The default is 5, i.e. at least 5 MADs away from the median. Choose higher values if you want to filter out less samples. This parameter is particularly important for single cells to ensure reliable parameter estimation.",
      function_name = "estimateParam"
    ),
    param_numeric(
      id = "sigma",
      default = 1.96,
      lower = 0,
      border = FALSE,
      process = "estimation",
      description = "The variability band width for mean-dispersion loess fit defining the prediction interval for read count simulation. Default is 1.96, i.e. 95% interval.",
      function_name = "estimateParam"
    ),
    param_integer(
      id = "NCores",
      default = 1L,
      lower = 1L,
      border = TRUE,
      process = "estimation",
      description = "The number of cores for normalisation method SCnorm and Census. The default NULL means 1 core.",
      function_name = c("estimateParam", "simulateDE")
    ),
    param_character(
      id = "Normalisation",
      default = NULL,
      alternatives = c('depth','none'),
      force = TRUE,
      process = "estimation",
      description = "Normalisation is a character value: 'depth' or 'none'. For more information, please consult the details section. This is for estimateSpike function in powsimR",
      function_name = "estimateSpike"
    ),
    param_numeric(
      id = "SampleFilter",
      default = 3,
      lower = 0,
      border = FALSE,
      process = "estimation",
      description = "SampleFilter is a numeric vector indicating the minimal number of MADs (median absolute deviation) away from the median number of features detected as well as sequencing depth across all samples so that outlying samples are removed prior to normalisation and parameter estimation. The default is 5, i.e. at least 5 MADs away from the median. Choose higher values if you want to filter out less samples. This parameter is particularly important for single cells to ensure reliable parameter estimation. This is for estimateSpike function in powsimR",
      function_name = "estimateSpike"
    ),
    param_integer(
      id = "nsims",
      default = 25L,
      lower = 1L,
      border = TRUE,
      description = "Number of simulations to run. Default is 25.",
      function_name = "Setup"
    ),
    param_others(
      id = "ngenes",
      type = "integer",
      default = NULL,
      description = "ngenes is a numeric vector specifying the number of genes to simulate. Default is NULL, i.e. the number of genes that were deemed expressed after filtering.",
      function_name = "Setup"
    ),
    param_numeric(
      id = "p.DE",
      default = 0.1,
      lower = 0,
      upper = 1,
      border = TRUE,
      description = "Numeric vector between 0 and 1 representing the percentage of genes being differentially expressed due to phenotype, i.e. biological signal. Default is 0.1.",
      function_name = "Setup"
    ),
    param_others(
      id = "pLFC",
      type = c("numeric, vector, function"),
      default = 1,
      description = "The log phenotypic fold change for DE genes. This can be: (1) a constant, e.g. 2; (2) a vector of values with length being number of DE genes. If the input is a vector and the length is not the number of DE genes, it will be sampled with replacement to generate log fold changes; (3) a function that takes an integer n, and generates a vector of length n, e.g. function(x) rnorm(x, mean=0, sd=1.5). Default is 1. Please note that the fold change should be on log2 scale!",
      function_name = "Setup"
    ),
    param_vector(
      id = "p.G",
      default = 1,
      description = "The log phenotypic fold change for DE genes. This can be: (1) a constant, e.g. 2; (2) a vector of values with length being number of DE genes. If the input is a vector and the length is not the number of DE genes, it will be sampled with replacement to generate log fold changes; (3) a function that takes an integer n, and generates a vector of length n, e.g. function(x) rnorm(x, mean=0, sd=1.5). Default is 1. Please note that the fold change should be on log2 scale!",
      function_name = "Setup"
    ),
    param_vector(
      id = "p.B",
      description = "Numeric vector between 0 and 1 representing the percentage of genes being differentially expressed between batches. Default is NULL, i.e. no batch effect.",
      function_name = "Setup"
    ),
    param_others(
      id = "bLFC",
      type = c("numeric, vector, function"),
      default = NULL,
      description = "The log batch fold change for all genes. This can be: (1) a constant, e.g. 2; (2) a vector of values with length being number of all genes. If the input is a vector and the length is not the number of total genes, it will be sampled with replacement to generate log fold changes; (3) a function that takes an integer n, and generates a vector of length n, e.g. function(x) rnorm(x, mean=0, sd=1.5). Note that the simulations of only two batches is implemented. Default is NULL, i.e. no batch effect. Please note that the fold change should be on log2 scale!",
      function_name = "Setup"
    ),
    param_character(
      id = "bPattern",
      default = "uncorrelated",
      alternatives = c("uncorrelated", "orthogonal", "correlated"),
      description = "Character vector for batch effect pattern if p.B is non-null. Possible options include: 'uncorrelated', 'orthogonal' and 'correlated'. Default is 'uncorrelated'.",
      function_name = "Setup"
    ),
    param_vector(
      id = "n1",
      default = c(20,50,100),
      description = "Integer vectors specifying the number of biological replicates in each group. Default values are n1=c(20,50,100). The vectors need to have the same number of entries, i.e. length.",
      function_name = "Setup"
    ),
    param_vector(
      id = "n2",
      default = c(30,60,120),
      description = "Integer vectors specifying the number of biological replicates in each group. Default values are n2=c(30,60,120). The vectors need to have the same number of entries, i.e. length.",
      function_name = "Setup"
    ),
    param_vector(
      id = "Thinning",
      default = NULL,
      description = "Numeric vector specifying the downsampling. It has to be the same length and order as the vector n1. This is an implementation of thinCounts function. The vector entries should be a proportion between 0 and 1, e.g. 0.75 means that each sample in that group will have on average 75% sequencing depth compared to the original sequencing depth as listed in estimateParam function. Note that no upsampling is possible, i.e. defining a proportion greater than 1. The default is NULL, meaning no downsampling.",
      function_name = "Setup"
    ),
    param_character(
      id = "LibSize",
      default = "equal",
      alternatives = c("equal", "given"),
      description = "Size factors representing sample-specific differences/biases in expected mean values of counts: 'equal' or 'given'. The default is 'equal', i.e. equal size factor of 1. If the user defines it as 'given', the size factors are sampled from the estimated size factors of estimateParam function.",
      function_name = "Setup"
    ),
    param_others(
      id = "estParamRes",
      force = TRUE,
      type = "estParamRes",
      description = "The estimated simulation parameters for genes. This can be: (1) The output of estimateParam function. (2) A string specifying the name of precalculated estimates.",
      function_name = "Setup"
    ),
    param_others(
      id = "estSpikeRes",
      type = "estSpikeRes",
      description = "The spike-in simulation parameters generated by estimateSpike. These are needed for applying spike-in-dependent normalisation methods. Default is NULL, i.e. no spike-in count simulations.",
      function_name = "Setup"
    ),
    param_others(
      id = "DropGenes",
      type = "Boolean",
      description = "By default, the estimated parameters and fit based on genes that were defined as expressed are used for simulations. By setting this parameter to TRUE, a fraction of genes will be dropouts. The dropout genes are defined in estimateParam function using the GeneFilter option and can be plotted with plotParam function. Default is FALSE, i.e. no gene expression dropouts.",
      function_name = "Setup"
    ),
    param_numeric(
      id = "setup.seed",
      default = 111,
      lower = 0,
      description = "Setup seed.",
      function_name = "Setup"
    ),
    param_others(
      id = "SetupRes",
      type = "SetupRes",
      description = "This object specifies the simulation setup. This must be the return object from Setup function.",
      function_name = "simulateDE"
    ),
    param_character(
      id = "Prefilter",
      default = NULL,
      alternatives = c("CountFilter", "FreqFilter"),
      description = "A character vector specifying the gene expression filtering method to be used prior to normalisation (and possibly imputation). Default is NULL, i.e. no filtering.",
      function_name = "simulateDE"
    ),
    param_character(
      id = "Imputation",
      default = NULL,
      alternatives = c("scImpute", "DrImpute", "SAVER", "scone", "MAGIC"),
      description = "A character vector specifying the gene expression imputation method to be used prior to normalisation. Default is NULL, i.e. no imputation.",
      function_name = "simulateDE"
    ),
    param_character(
      id = "Normalisation",
      force = TRUE,
      alternatives = c("TMM", "UQ", "MR", "PosCounts", "scran", "SCnorm", "Linnorm", "sctransform", "Census", "depth"),
      description = "Normalisation method to use.",
      function_name = "simulateDE"
    ),
    param_vector(
      id = "Label",
      default = "none",
      description = "A character vector to define whether information about group labels should be used for normalisation. This is only implemented for scran and SCnorm. Possible options include the default 'none' which means that no sample group information is considered for normalisation; 'known' means that the simulated group labels are used and 'clustering' which applies an unsupervised hierarchical clustering to determine the group labels.",
      function_name = "simulateDE"
    ),
    param_character(
      id = "DEmethod",
      force = TRUE,
      alternatives = c("T-Test", "limma-trend", "limma-voom", "edgeR-LRT", "edgeR-QL", "DESeq2", "ROTS", "sctransform", "baySeq", "NOISeq", "EBSeq", "MAST", "BPSC", "scDD", "DECENT", "edgeR-zingeR", "DESeq2-zingeR", "edgeR-ZINB-WaVE", "DESeq2-ZINB-WaVE"),
      description = "A character vector specifying the DE detection method to be used. Please consult the Details section for available options.",
      function_name = "simulateDE"
    ),
    param_Boolean(
      id = "DEFilter",
      default = FALSE,
      description = "A logical vector indicating whether to run DE testing on filtered and/or imputed count data. Default is FALSE.",
      function_name = "simulateDE"
    ),
    param_Boolean(
      id = "Counts",
      default = FALSE,
      description = "A logical vector indicating whether the simulated count matrix is also provided as output. Default is FALSE since the output can be quite large. Note that if DEFilter is TRUE, then the returned count matrix will countain the filtered and/or imputed count data.",
      function_name = "simulateDE"
    )
  )

  powsimR_method <- method_definition(
    method = "powsimR",
    programming = "R",
    url = "https://github.com/bvieth/powsimR",
    authors = authors_definition(
      first = "Beate",
      last = "Vieth",
      email = NULL,
      github = "https://github.com/bvieth/powsimR",
      orcid = NULL
    ),
    manuscript = manuscript_definition(
      title = "powsimR: power analysis for bulk and single cell RNA-seq experiments",
      doi = "10.1093/bioinformatics/btx435",
      journal = "Bioinformatics",
      date = "2017",
      peer_review = TRUE
    ),
    description = "Power analysis for bulk and single cell RNA-seq experiments")

  list(powsimR_method = powsimR_method,
       powsimR_parameters = powsimR_parameters)
}
