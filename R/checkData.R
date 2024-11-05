# chooses how to put inputs together, verifies data
checkData <- function(gene_expression_matrix = NULL,
                        ground_truth = NULL,
                        gene_association_matrix = NULL,
                        rf_features = NULL,
                        rf_outputs = NULL,
                        gene_regulatory_network = NULL,
                        in_data_obj = NULL,
                        config_path = NULL,
                        in_format = "separate",
                        suppress_warnings = FALSE,
                        ...) {
  # TODO: warnings are not actually suppressed
  tryCatch({
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
    else if(in_format == "config_file") {
      # TODO: pass in using config
      out_data_obj <- list(gene_expression_matrix = gene_expression_matrix,
                           ground_truth = ground_truth,
                           gene_association_matrix = gene_association_matrix,
                           rf_features = rf_features,
                           rf_outputs = rf_outputs,
                           gene_regulatory_network = gene_regulatory_network)
    }
    else {
      stop("Invalid input format specified")
    }

    names = c("gene_expression_matrix", "ground_truth","gene_association_matrix", "rf_features", "rf_outputs", "gene_regulatory_network")

    for(name in names) {
      if(!switch(name, "gene_expression_matrix" = validateGEM(out_data_obj$gene_expression_matrix, suppress_warnings = suppress_warnings),
                          "ground_truth" = validateGT(out_data_obj$ground_truth, suppress_warnings = suppress_warnings),
                          "gene_association_matrix" = validateGAM(out_data_obj$gene_association_matrix, suppress_warnings = suppress_warnings),
                          "rf_features" = validateRFFeatures(out_data_obj$rf_features, suppress_warnings = suppress_warnings),
                          "rf_outputs" = validateRFOut(out_data_obj$rf_outputs, suppress_warnings = suppress_warnings),
                          "gene_regulatory_network" = validateGRN(out_data_obj$gene_regulatory_network, suppress_warnings = suppress_warnings))) {
        stop(name, " is invalid")
      }
    }
    return(out_data_obj)
  }, error = function(e) {
    message("Error: ", e$message)
    return(NULL)
  })
}

validateGEM <- function(gem, suppress_warnings = FALSE) {
  tryCatch({
    if(!is.data.frame(gem) & !is.matrix(gem)) {
      stop("Gene expression matrix must be a matrix or DataFrame")
    }

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
    message("Error: ", e$message)
    return(FALSE)
  })
}

validateGT <- function(gt, suppress_warnings = FALSE) {
  # TODO: migrate traintestreport ground truth-gene expression matrix gene name checking here
  # TODO: can have extra columns but check that there exist cols labeled with the names specified in the gc
  tryCatch({
    if(is.null(gt)) {
      if(!suppress_warnings) {
        warning("Ground truth table was not given. trainTestReport will not work.")
      }
    }

    return(TRUE)
  }, error = function(e) {
    message("Error: ", e$message)
    return(FALSE)
  }, warning = function(w) {
    message("Warning: ", w$message)
    if(askUserProceed()) {
      return(TRUE)
    }
    return(FALSE)
  })
}

# TODO: validation of other inputs so checkData can be used for general data checking

validateGAM <- function(gam, suppress_warnings = FALSE) {
  return(TRUE)
}

validateRFFeatures <- function(rf_features, suppress_warnings = FALSE) {
  return(TRUE)
}

validateRFOut <- function(rf_out, suppress_warnings = FALSE) {
  return(TRUE)
}

validateGRN <- function(grn, suppress_warnings = FALSE) {
  return(TRUE)
}

askUserProceed <- function() {
  input <- ""
  while(input != "y" && input != "n") {
    input = readline(prompt = "Do you want to proceed? (y/n)")
  }
  if(input == "y") {
    return(TRUE)
  }
  return(FALSE)
}
