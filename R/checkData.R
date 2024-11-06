# verifies data
# data that is valid is returned in the object
# data that is not valid is null in the object and the reason it is not valid is printed
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
  # TODO: see if passing in ... to the data obj makes sense
  tryCatch({
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

    names = c("gene_expression_matrix", "ground_truth")

    for(name in names) {
      check_result <- switch(name, "gene_expression_matrix" = checkGEM(out_data_obj$gene_expression_matrix),
                             "ground_truth" = checkGT(out_data_obj$ground_truth),
                             "gene_association_matrix" = checkGAM(out_data_obj$gene_association_matrix),
                             "rf_features" = checkRFFeatures(out_data_obj$rf_features),
                             "rf_outputs" = checkRFOut(out_data_obj$rf_outputs),
                             "gene_regulatory_network" = checkGRN(out_data_obj$gene_regulatory_network))
      if(!check_result$valid) {
        print(paste0("Invalid input: ", check_result$message))
        out_data_obj[name] <- NULL
      }
      if(!suppress_warnings & check_result$warning) {
        print(paste0(check_result$message))
      }
    }
    return(out_data_obj)
  }, error = function(e) {
    message(e$message)
    return(NULL)
  })
}

checkGEM <- function(gem) {
  valid <- TRUE
  warning <- FALSE
  message <- "Gene expression matrix is valid"

  if(is.null(gem)) {
    valid <- FALSE
    message <- "Gene expression matrix was not given"
    return(list(valid = valid, warning = warning, message = message))
  }


  if(!is.data.frame(gem) & !is.matrix(gem)) {
    valid <- FALSE
    message <- "Gene expression matrix must be a matrix or DataFrame"
    return(list(valid = valid, warning = warning, message = message))
  }

  if(is.data.frame(gem)) {
    if(!all(sapply(gem, is.numeric))) {
      valid <- FALSE
      message <- "Gene expression matrix must have numeric data"
      return(list(valid = valid, warning = warning, message = message))
    }
  }
  else {
    if(!is.numeric(gem)) {
      valid <- FALSE
      message <- "Gene expression matrix must have numeric data"
      return(list(valid = valid, warning = warning, message = message))
    }
  }

  return(list(valid = valid, warning = warning, message = message))
}

checkGT <- function(gt) {
  # TODO: migrate traintestreport ground truth-gene expression matrix gene name checking here
  # TODO: can have extra columns but check that there exist cols labeled with the names specified in the gc

  valid <- TRUE
  warning <- FALSE
  message <- "Ground truth table is valid"

  if(is.null(gt)) {
    warning <- TRUE
    message <- "Ground truth table was not given. trainTestReport will not work."
    return(list(valid = valid, warning = warning, message = message))
  }

  if(!("src" %in% colnames(gt) & "trgt" %in% colnames(gt))) {
    valid <- FALSE
    message <- "Ground truth table must contain 'src' and 'trgt' columns"
    return(list(valid = valid, warning = warning, message = message))
  }

  return(list(valid = valid, warning = warning, message = message))
}
