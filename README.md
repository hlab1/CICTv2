# CICTv2
A package for gene regulatory network discovery from single-cell gene expression data using the CICT (Causal Inference Using Composition of Transactions) method.

# Installation
With the `remotes` package installed:
```
remotes::install_github("hlab1/CICTv2")
```

# Getting started with CICT

## 1. Loading test data
Test data includes a gene expression dataset and ground-truth network from the SERGIO single-cell simulator benchmarking set (https://doi.org/10.1016/j.cels.2020.08.003). This data can be found as rdata files under `./data` or as csv files in `./inst/extdata`. We have also included an example config file in `./inst/extdata`.

## 2. Running CICT
### Running CICT with separate R inputs
```
gene_expression_matrix <- read.csv("SERGIO_DS4_net0_gene_expression_matrix.csv", header = TRUE, row.names = 1)
ground_truth <- read.table("SERGIO_DS4_net0_ground_truth.csv",  header=TRUE, sep = ",")
out <- runCICT(gene_expression_matrix = gene_expression_matrix, ground_truth = ground_truth)
```

### Running CICT with config file
```
runCICT(config_path = "SERGIO_DS4_net0_config.yaml", in_format = "config_file")
```