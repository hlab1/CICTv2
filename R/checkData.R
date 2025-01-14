# TODO: may be better to change from valid-warning-message system now that error handling has been updated
# TODO: make sure raw_edges documentation is correct

# TODO: confirm type of grn and add to docs
# TODO: see where vignettes are more appropriate than function docs
# TODO: more in-depth explanation for other options parameters

#' Check validity of data
#'
#' Checks if inputted data is either valid for use in CICT functions or is a
#' valid output of a CICT function. Valid data is returned in the CICT data
#' list. Invalid data has a null field in the list. Does not currently verify
#' `raw_edges`, `edge_features`, `model`, `model_assessment`, or
#' `predicted_edges`.
#'
#' @param gene_expression_matrix The gene expression matrix, where each row
#'   represents a gene and each column represents a sample. A matrix or
#'   DataFrame with numeric data.
#' @param ground_truth A DataFrame representing the ground-truth gene regulatory
#'   network for model training and evaluation. Each row represents a
#'   source-target relationship, with the source gene in the column labeled
#'   `"src"` and the target gene in the column labeled `"trgt"`
#' @param raw_edges Edges for the gene association matrix calculated from
#'   expressioon data. Each row represents a source-target relationship, with
#'   the source gene in the column labeled `"src"`, the target gene in the
#'   column labeled `"trgt"`, and the correlation coefficient in a row named
#'   with `Pearson` or `Spearman` depending on the type of coefficient used.
#' @param edge_features Features calculated from raw edge data to be used in model
#'   prediction
#' @param model The model used to predict edges
#' @param model_assessment Evaluation of the model used to predict edges
#' @param predicted_edges A predicted gene regulatory network.
#' @param in_data_obj A list in the CICT data object format. Produced by a CICT
#'   function.
#' @param config_path Path to the YAML config file.
#' @param in_format String indicating expected input format. `"separate"` if
#'   passing inputs as separate arguments, `"data_obj"` if passing inputs
#'   through `in_data_obj`, `"config"` if passing inputs through `config_path`.
#' @param ... Other options
#'
#' @return A list in the CICT data object format, where fields for which valid
#'   data were provided contain the input data and fields for which data
#'   provided was invalid are null.
#' @noRd
#'
#' @examples
#' # separate inputs
#' checkData(gene_expression_matrix = SERGIO_DS4_gene_expression_matrix, ground_truth = SERGIO_DS4_ground_truth)
#'
#' # input a data list
#' checkData(in_data_obj = checkData(gene_expression_matrix = SERGIO_DS4_gene_expression_matrix,
#'                         ground_truth = SERGIO_DS4_ground_truth),
#'           in_format = "data_obj")
#'
#' # inputs from a config file
#' # create a config YAML file with absolute paths to data
#' yaml::write_yaml(list(gene_expression_matrix = system.file("extdata", "SERGIO_DS4_gene_expression_matrix.csv", package = "CICTv2", mustWork = TRUE),
#'                 ground_truth = system.file("extdata", "SERGIO_DS4_ground_truth.csv", package = "CICTv2", mustWork = TRUE),
#'                 raw_edges = NULL,
#'                 edge_features = NULL,
#'                 model = NULL,
#'                 model_assessment = NULL,
#'                 predicted_edges = NULL),
#'            "/home/syz248/CICTv2/inst/extdata/SERGIO_DS4_config.yaml")
#' checkData(config_path = system.file("extdata", "SERGIO_DS4_config.yaml", package = "CICTv2", mustWork = TRUE), in_format = "config_file")
checkData <- function(gene_expression_matrix = NULL,
                      ground_truth = NULL,
                      raw_edges = NULL,
                      edge_features = NULL,
                      model = NULL,
                      model_assessment = NULL,
                      predicted_edges = NULL,
                      in_data_obj = NULL,
                      config_path = NULL,
                      in_format = "separate",
                      suppress_warnings = FALSE,
                      ...) {
  # TODO: decide how to pass in ... to the data obj
  tryCatch({
    if (in_format == "separate") {
      out_data_obj <-
        list(
          gene_expression_matrix = gene_expression_matrix,
          ground_truth = ground_truth,
          raw_edges = raw_edges,
          edge_features = edge_features,
          model = model,
          model_assessment = model_assessment,
          predicted_edges = predicted_edges
        )
    }
    else if (in_format == "data_obj") {
      out_data_obj <- in_data_obj
    }
    else if (in_format == "config_file") {
      # TODO: check if config path is a string
      config <- yaml::yaml.load_file(config_path)
      # TODO: allow passing in of read.csv arguments for more flexibility
      out_data_obj <-
        list(
          gene_expression_matrix = (
            if (!is.null(config$gene_expression_matrix))
              read.csv(config$gene_expression_matrix, row.names = 1)
            else
              NULL
          ),
          ground_truth = (if (!is.null(config$ground_truth))
            read.csv(config$ground_truth)
            else
              NULL),
          raw_edges = (if (!is.null(config$raw_edges))
            read.csv(config$raw_edges)
            else
              NULL),
          edge_features = (if (!is.null(config$edge_features))
            read.csv(config$edge_features)
            else
              NULL),
          model = (if (!is.null(config$model))
            read.csv(config$model)
            else
              NULL),
          model_assessment = (if (!is.null(config$model_assessment))
            read.csv(config$model_assessment)
            else
              NULL),
          predicted_edges = (if (!is.null(config$predicted_edges))
            read.csv(config$predicted_edges)
            else
              NULL)
        )
    }
    else {
      stop("Invalid input format specified")
    }

    names = c("gene_expression_matrix", "ground_truth")

    for (name in names) {
      check_result <-
        switch(
          name,
          "gene_expression_matrix" = checkGEM(out_data_obj$gene_expression_matrix),
          "ground_truth" = checkGT(out_data_obj$ground_truth)
        )
      if (!check_result$valid) {
        print(paste0("Invalid input: ", check_result$message))
        out_data_obj[name] <- NULL
      }
      if (!suppress_warnings & check_result$warning) {
        print(paste0(check_result$message))
      }
    }
    return(out_data_obj)
  }, error = function(e) {
    message(e$message)
    return(NULL)
  })
}

#' Check validity of gene expression matrix
#'
#' Checks if the gene expression matrix exists, is a matrix or DataFrame, and
#' has numeric data.
#'
#' @param gem The object that will be checked to determine whether it is a valid
#'   gene expression matrix.
#'
#' @return A list with fields `valid`, `warning`, and `error`. `valid` is a
#'   boolean that is `TRUE` if the input exists, is a matrix or DataFrame, and
#'   has numeric data; it is otherwise `FALSE`. `warning` is a boolean with
#'   value `FALSE`. `message` is a string that contains an explanation of why
#'   the input triggered a warning, was valid, or was invalid.
#' @noRd
#'
#' @examples
#' checkGEM(SERGIO_DS4_gene_expression_matrix)
checkGEM <- function(gem) {
  # TODO: check for gene and sample names that exist and are unique
  valid <- TRUE
  warning <- FALSE
  message <- "Gene expression matrix is valid"

  if (is.null(gem)) {
    valid <- FALSE
    message <- "Gene expression matrix was not given"
    return(list(
      valid = valid,
      warning = warning,
      message = message
    ))
  }


  if (!is.data.frame(gem) & !is.matrix(gem)) {
    valid <- FALSE
    message <-
      "Gene expression matrix must be a matrix or DataFrame"
    return(list(
      valid = valid,
      warning = warning,
      message = message
    ))
  }

  if (is.data.frame(gem)) {
    if (!all(sapply(gem, is.numeric))) {
      valid <- FALSE
      message <- "Gene expression matrix must have numeric data"
      return(list(
        valid = valid,
        warning = warning,
        message = message
      ))
    }
  }
  else {
    if (!is.numeric(gem)) {
      valid <- FALSE
      message <- "Gene expression matrix must have numeric data"
      return(list(
        valid = valid,
        warning = warning,
        message = message
      ))
    }
  }

  gene_names = rownames(gem)
  if (is.null(gene_names) |
      !identical(unique(gene_names), gene_names)) {
    valid <- FALSE
    message <- "Gene expression matrix must have unique gene names"
    return(list(
      valid = valid,
      warning = warning,
      message = message
    ))
  }

  sample_names = colnames(gem)
  if (is.null(sample_names) |
      !identical(unique(sample_names), sample_names)) {
    valid <- FALSE
    message <-
      "Gene expression matrix must have unique sample names"
    return(list(
      valid = valid,
      warning = warning,
      message = message
    ))
  }

  return(list(
    valid = valid,
    warning = warning,
    message = message
  ))
}

#' Check validity of ground truth table
#'
#' Checks if the ground truth table exists and contains the columns needed for
#' model training and evaluation.
#'
#' @param gt The object that will be checked to determine whether it is a valid
#'   ground truth table.
#'
#' @return A list with fields `valid`, `warning`, and `error`. `valid` is a
#'   boolean that is `TRUE` if the input exists, and contains the columns "src"
#'   and "trgt"; it is otherwise `FALSE`. `warning` is a boolean that is `TRUE`
#'   if the ground truth table does not exists, and otherwise `FALSE`. `message`
#'   is a string that contains an explanation of why the input triggered a
#'   warning, was valid, or was invalid.
#' @noRd
#'
#' @examples
#' checkGT(SERGIO_DS4_ground_truth)
checkGT <- function(gt) {
  # TODO: migrate ground truth-gene expression matrix gene name checking from traintestreport to here
  # TODO: do we want to enforce tibble?
  valid <- TRUE
  warning <- FALSE
  message <- "Ground truth table is valid"

  if (is.null(gt)) {
    warning <- TRUE
    message <-
      "Ground truth table was not given. trainTestReport will not work."
    return(list(
      valid = valid,
      warning = warning,
      message = message
    ))
  }

  if (!("src" %in% colnames(gt) & "trgt" %in% colnames(gt))) {
    valid <- FALSE
    message <-
      "Ground truth table must contain 'src' and 'trgt' columns"
    return(list(
      valid = valid,
      warning = warning,
      message = message
    ))
  }

  return(list(
    valid = valid,
    warning = warning,
    message = message
  ))
}
