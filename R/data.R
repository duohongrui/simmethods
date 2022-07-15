#' External ERCC information of 92 molecules.
#'
#' A dataset containing the information of 92 molecules ERCC.
#'
#' @format A data frame with 92 rows and 6 variables:
#' \describe{
#'   \item{\code{ERCC_id}}{The name of the ERCC molecules}
#'   \item{\code{subgroup}}{Which group does the molecule belongs to}
#'   \item{\code{con_Mix1_attomoles_ul}}{The concentration of the molecules in Mix1 liquid}
#'   \item{\code{con_Mix2_attomoles_ul}}{The concentration of the molecules in Mix2 liqui}
#'   \item{\code{expected_fc}}{The expected fold change between Mix1 and Mix2}
#'   \item{\code{log2_fc}}{Log2 transformation of the \code{expected_fc}}
#' }
#' @source \url{https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt}
"ERCC_info"


#' Single-cell RNA-seq Data from GSE54695
#'
#' @description A matrix contains the single-cell gene expression data from
#' GSE54695. Total 160 cells are cultured in two types of medium, two-inhibitor (2i)
#' and serum. In addtion, the count matrix contains ERCC spike-in information and
#' it can be very important when uses are intended to simulate datasets by BEARscc
#' and other methods.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54695}
"data"


#' Group information of 160 Cells in GSE54695
#'
#' @description Group information of 160 clls in GSE54695. Total 160 cells are
#' cultured in two types of medium, 0 for two-inhibitor (2i) and 1 for serum.
#'
"group_condition"





