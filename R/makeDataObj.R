makeDataObj <- function(gene_expression_matrix = NULL,
                        ground_truth = NULL,
                        gene_association_matrix = NULL,
                        rf_features = NULL,
                        rf_outputs = NULL,
                        gene_regulatory_network = NULL,
                        in_data_obj = NULL,
                        suppress_warnings = FALSE,
                        ...)
{
  # warn in the docs that user makes their own cict_data_obj at their own risk
  # see if passing in ... to the data obj makes sense
  out_data_obj <- list(gene_expression_matrix = NULL,
                       ground_truth = NULL,
                       gene_association_matrix = NULL,
                       rf_features = NULL,
                       rf_outputs = NULL,
                       gene_regulatory_network = NULL)

  names = c("gene_expression_matrix", "ground_truth","gene_association_matrix", "rf_features", "rf_outputs", "gene_regulatory_network")
  for (name in names) {
    # TODO: prevent nulls from getting taken out of the list in this step
    out_data_obj[[name]] <- getInput(passed_in = get(name), obj_in = in_data_obj[[name]], name = name)
  }

  return(out_data_obj)
}

getInput <- function(passed_in, obj_in, name) {
  # possibly print some indicators of which is being used
  if(!is.null(passed_in)) {
    input <- passed_in
  }
  else {
    input <- obj_in
  }

  if(!is.null(obj_in) & !is.null(passed_in)) {
    print(paste("A", name, "was passed in as input, but the input cict_data_obj already contains a", name))
    askUserProceed()
  }

  if(validateInput(input, name)) {
    print(paste("The", name, "is valid"))
  }
  else {
    print(paste("The", name, "is invalid. Set anyway?"))
    askUserProceed()
  }

  return(input)
}

validateInput <- function(input, name) {
  if(name == "gene_expression_matrix") {
    return(validateGEM(input))
  }
  if(name == "ground_truth") {
    return(validateGT(input))
  }
  if(name == "gene_association_matrix") {
    return(validateGAM(input))
  }
  if(name == "rf_outputs") {
    return(validateRFOut(input))
  }
  if(name == "rf_features") {
    return(validateRFFeatures(input))
  }
  if(name == "gene_regulatory_network") {
    return(validateGRN(input))
  }
  return(FALSE)
}

validateGEM <- function(gem) {
  return(TRUE)
}

validateGAM <- function(gem) {
  return(TRUE)
}

validateGT <- function(gem) {
  return(TRUE)
}

validateRFFeatures <- function(gem) {
  return(TRUE)
}

validateRFOut <- function(gem) {
  return(TRUE)
}

validateGRN <- function(gem) {
  # only useful if we want to make the data creation function an overall validation function
  return(TRUE)
}

askUserProceed <- function() {
  print("Do you want to proceed? (y/n)")
  print("Proceed")
}
