% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateRawEdges.R
\name{calculateRawEdges}
\alias{calculateRawEdges}
\title{calculateRawEdges}
\usage{
calculateRawEdges(
  gene_expression_matrix = NULL,
  cict_raw_edge_col = "Spearman",
  n.workers = 5,
  ...
)
}
\arguments{
\item{gene_expression_matrix}{Data.frame of gene expression where rows are
genes and columns are samples}

\item{cict_raw_edge_col}{gene-gene association metric to calculate}

\item{n.workers}{Integer. If running parallel, number of threads on which to
run calculations}

\item{...}{Options to be passed to subfunctions in calculateRawEdges}
}
\value{
A list in the CICT data object format. Contains
\code{gene_expression_matrix} and \code{raw_edges}.
}
\description{
Returns a gene-gene association matrix from a gene expression matrix
}
\examples{
# Download data from the external data folder of the CICTv2 GitHub repo
download.file("https://raw.githubusercontent.com/hlab1/CICTv2/refs/heads/main/inst/extdata/SERGIO_DS4_net0_gene_expression_matrix.csv", "SERGIO_DS4_net0_gene_expression_matrix.csv")
gene_expression_matrix <- read.csv("SERGIO_DS4_net0_gene_expression_matrix.csv", header = TRUE, row.names = 1)

calculateRawEdges(gene_expression_matrix = gene_expression_matrix)

# Reset workspace
unlink("SERGIO_DS4_net0_gene_expression_matrix.csv")
rm(gene_expression_matrix)
}
