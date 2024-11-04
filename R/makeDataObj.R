# chooses how to put inputs together, verifies data
makeDataObj <- function(gene_expression_matrix = NULL,
                        ground_truth = NULL,
                        gene_association_matrix = NULL,
                        rf_features = NULL,
                        rf_outputs = NULL,
                        gene_regulatory_network = NULL,
                        in_data_obj = NULL,
                        in_format = "separate",
                        ...) {
  # TODO: see if passing in ... to the data obj makes sense
  if(in_format == "separate") {
    out_data_obj <- list(gene_expression_matrix = gene_expression_matrix,
                         ground_truth = ground_truth,
                         gene_association_matrix = gene_association_matrix,
                         rf_features = rf_features,
                         rf_outputs = rf_outputs,
                         gene_regulatory_network = gene_regulatory_network)
  }
  else if(in_format == "data_obj") {
    # warn in the docs that user makes their own cict_data_obj at their own risk
    out_data_obj <- in_data_obj
  }
  else if(in_format == config) {
    # TODO: pass in using config
    out_data_obj <- list(gene_expression_matrix = gene_expression_matrix,
                         ground_truth = ground_truth,
                         gene_association_matrix = gene_association_matrix,
                         rf_features = rf_features,
                         rf_outputs = rf_outputs,
                         gene_regulatory_network = gene_regulatory_network)
  }

  names = c("gene_expression_matrix", "ground_truth","gene_association_matrix", "rf_features", "rf_outputs", "gene_regulatory_network")

  tryCatch({
    for(name in names) {
      if(!switch(name, "gene_expression_matrix" = validateGEM(out_data_obj$gene_expression_matrix),
                          "ground_truth" = validateGT(out_data_obj$gene_expression_matrix),
                          "gene_association_matrix" = validateGAM(out_data_obj$gene_association_matrix),
                          "rf_features" = validateRFFeatures(out_data_obj$rf_features),
                          "rf_outputs" = validateRFOut(out_data_obj$rf_outputs),
                          "gene_regulatory_network" = validateGRN(out_data_obj$gene_regulatory_network))) {
        stop(name, " is invalid")
      }
    }
    return(out_data_obj)
  }, error = function(e) {
    message(e$message)
  })
  return(NULL)
}

validateGEM <- function(gem) {
  tryCatch({
    # TODO: check if this is necessary or if dims are all that matter
    if(!is.data.frame(gem) & !is.matrix(gem)) {
      stop("Gene expression matrix must be a matrix or DataFrame")
    }

    # TODO: use this instead if any data type with 2 dims can be used
    # if(length(dim(gem)) != 2) {
    #   stop("Gene expression matrix must have 2 dimensions")
    # }

    if(is.data.frame(gem)) {
      if(!all(sapply(gem, is.numeric))) {
        stop("Gene expression matrix must have numeric data")
      }
    }
    else {
      if(!is.numeric(gem)) {
        stop("Gene expression matrix must have numeric data")
      }
    }

    return(TRUE)
  }, error = function(e) {
    message(e$message)
  })
  return(FALSE)
}

validateGT <- function(gt) {
  return(TRUE)
}

validateGAM <- function(gam) {
  print("Gene association matrix is valid")
  return(TRUE)
}

validateRFFeatures <- function(rf_features) {
  return(TRUE)
}

validateRFOut <- function(rf_out) {
  return(TRUE)
}

validateGRN <- function(grn) {
  # only useful for overall validation
  return(TRUE)
}

askUserProceed <- function() {
  print("Do you want to proceed? (y/n)")
  print("Proceed")
}
