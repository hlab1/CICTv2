runCICT <- function(gene_expression_matrix = NULL,
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
  cict_data_obj <- checkData(gene_expression_matrix = gene_expression_matrix,
                               ground_truth = ground_truth,
                               gene_association_matrix = gene_association_matrix,
                               rf_features = rf_features,
                               rf_outputs = rf_outputs,
                               gene_regulatory_network = gene_regulatory_network,
                               in_data_obj = in_data_obj,
                               config_path = "",
                               in_format = in_format,
                               suppress_warnings = suppress_warnings,
                               ...)
  # above function takes inputs, verifies them, and puts valid inputs in an object

  tryCatch({
    # gene expression matrix and ground truth are required for driver to complete correctly
    if(is.null(cict_data_obj) | is.null(cict_data_obj$gene_expression_matrix) | is.null(cict_data_obj$ground_truth)) {
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
