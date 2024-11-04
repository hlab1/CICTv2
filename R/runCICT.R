runCICT <- function(gene_expression_matrix = NULL, ground_truth = NULL, in_data_obj = NULL, suppress_warnings = FALSE, ...) {
  cict_data_obj <- makeDataObj(gene_expression_matrix = gene_expression_matrix,
                               ground_truth = ground_truth,
                               in_data_obj = in_data_obj,
                               suppress_warnings,
                               ...)
  # above function takes inputs, verifies them, and puts them in an object
  # resulting object is fed through below functions

  tryCatch({
    if(is.null(cict_data_obj)) {
      stop("Failed to create data object")
    }

    print("continued driver")
    # calculateRawEdges
    # prepareEdgeFeatures
    # trainTestReport

    return(cict_data_obj)
  }, error = function(e) {
    message(e$message)
  })
  return(NULL)
}
