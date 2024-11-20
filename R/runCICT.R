

# TODO: Write description of CICT data list format and link to it from function documentation
# TODO: add example data
# TODO: use example data to write usage examples with recommended practices (from wickham r packages book)
# TODO: more in-depth explanation for other options parameters

#' Driver for the CICT pipeline
#'
#' Takes a gene expression matrix and a ground truth table, calculates raw
#' edges, prepares edge features, and uses a random forest model to predict a
#' gene regulatory network. Does not currently work in config mode.
#'
#' @param gene_expression_matrix The gene expression matrix, where each row
#'   represents a gene and each column represents a sample. A matrix or
#'   DataFrame with numeric data.
#' @param ground_truth A DataFrame representing the ground-truth gene regulatory
#'   network for model training and evaluation. Each row represents a
#'   source-target relationship, with the source gene in the column labeled
#'   `"src"` and the target gene in the column labeled `"trgt"`
#' @param in_data_obj A list in the CICT data object format. Produced by a CICT
#'   function. Must contain `gene_expression_matrix` and `ground_truth`.
#' @param config_path Path to the YAML config file. Config must contain paths to
#'   the gene expression matrix and ground truth.
#' @param in_format String indicating expected input format. `"separate"` if
#'   passing inputs through `gene_expression_matrix` and `ground_truth`,
#'   `"data_obj"` if passing inputs through `in_data_obj`, `"config"` if passing
#'   inputs through `config_path`.
#' @param ... Other options
#'
#' @return A list in the CICT data object format, with the data from each step
#'   in the CICT pipeline. Contains `gene_expression_matrix`, `ground_truth`,
#'   `gene_association_matrix`, `edge_features`, `model`, `model_assessment`,
#'   `predicted_edges`, and potentially other data.
#' @export
#'
#' @examples print("TODO")
runCICT <- function(gene_expression_matrix = NULL,
                    ground_truth = NULL,
                    in_data_obj = NULL,
                    config_path = NULL,
                    in_format = "separate",
                    ...) {
  cict_data_obj <-
    checkData(
      gene_expression_matrix = gene_expression_matrix,
      ground_truth = ground_truth,
      in_data_obj = in_data_obj,
      config_path = config_path,
      in_format = in_format,
      ...
    )
  # above function takes inputs, verifies them, and puts valid inputs in an object

  tryCatch({
    # gene expression matrix and ground truth are required for driver to complete correctly
    if (is.null(cict_data_obj) |
        is.null(cict_data_obj$gene_expression_matrix) |
        is.null(cict_data_obj$ground_truth)) {
      stop("Failed to create data object")
    }

    print("continued driver")
    # calculateRawEdges
    # prepareEdgeFeatures
    # pass in raw_edges from cre
    pef_out <- prepareEdgeFeatures(
      gene_expression_matrix = cict_data_obj$gene_expression_matrix,
      ...
    )
    cict_data_obj$edge_features <- pef_out$edge_features

    pe_out <-
      predictEdges(
        gene_expression_matrix = cict_data_obj$gene_expression_matrix,
        edge_features = cict_data_obj$edge_features,
        ground_truth = cict_data_obj$ground_truth,
        ...
      )
    cict_data_obj$model <- pe_out$model
    cict_data_obj$model_assessment <- pe_out$model_assessment
    cict_data_obj$predicted_edges <- pe_out$predicted_edges

    return(cict_data_obj)
  }, error = function(e) {
    message(e$message)
  })
  return(NULL)
}
