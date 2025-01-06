# TODO: should we specify the source?
# TODO: document example data files
# TODO: meaning of net0 and net1

#' SERGIO DS4 gene expression matrix
#'
#' The DS4 gene expression matrix generated with SERGIO
#'
#' @format ## `SERGIO_DS4_net0_gene_expression_matrix`
#' A gzipped data frame with 100 rows and 900 columns:
#' \describe{
#'   \item{Rows}{Unique gene names}
#'   \item{Columns}{Unique sample names}
#'   \item{Values}{Expression levels of each gene in each sample}
#' }
#' @source <https://github.com/hlab1/scRNAseqWithCICT/blob/main/inputs/SERGIO_DS4/net0/ExpressionData.csv>
"SERGIO_DS4_net0_gene_expression_matrix"

#' SERGIO DS4 ground truth
#'
#' A ground truth gene regulatory network for the SERGIO DS4 generated gene expression data
#'
#' @format ## `SERGIO_DS4_net0_ground_truth`
#' A data frame with 137 rows and 2 columns:
#' \describe{
#'   \item{src}{Source gene}
#'   \item{trgt}{Target gene}
#' }
#' @source <https://github.com/hlab1/scRNAseqWithCICT/blob/main/inputs/SERGIO_DS4/net0/refNetwork.csv>
"SERGIO_DS4_net0_ground_truth"
