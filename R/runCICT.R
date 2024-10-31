runCICT <- function(gene_expression_matrix = NULL, ground_truth = NULL, in_data_obj = NULL, suppress_warnings = FALSE, ...) {
  cict_data_obj <- makeDataObj(gene_expression_matrix = gene_expression_matrix,
                                   ground_truth = ground_truth,
                                   in_data_obj = in_data_obj,
                                   suppress_warnings = suppress_warnings,
                                   ...)
  # validated object is fed through below functions

  # calculateRawEdges
  # prepareEdgeFeatures
  # trainTestReport
}
