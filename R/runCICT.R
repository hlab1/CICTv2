#' Driver for the CICT pipeline
#'
#' Takes a gene expression matrix and a ground truth table, calculates raw
#' edges, prepares edge features, and uses a random forest model to predict a
#' gene regulatory network. If in_format is specified as config_path, inputs and
#' arguments will come from a YAML config file.
#'
#' @param gene_expression_matrix The gene expression matrix, where each row
#'   represents a gene and each column represents a sample. A matrix or Data
#'   Frame with numeric data.
#' @param ground_truth A Data Frame representing the ground-truth gene
#'   regulatory network for model training and evaluation. Each row represents a
#'   source-target relationship, with the source gene in the column labeled
#'   `"src"` and the target gene in the column labeled `"trgt"`
#' @param config_path Path to the YAML config file, if one is used.
#' @param in_format String indicating expected input format. `"separate"` if
#'   passing inputs through `gene_expression_matrix` and `ground_truth`,
#'   `"config_file"` if passing inputs through `config_path`. If using the config
#'   format, all other arguments passed to the function are ignored except
#'   `config_path`, and the resulting predicted edges and CICT data object  will
#'   be stored as a CSV and .Rds file respectively. Use "results_dir" argument
#'   in config file to specify where these files should be saved; otherwise,
#'   they will be saved to the working directory.
#' @param ... Options to be passed to calculateRawEdges, prepareEdgeFeatures,
#'   and/or predictEdges
#'
#' @return A list in the CICT data object format, with the data from each step
#'   in the CICT pipeline. Contains `gene_expression_matrix`, `ground_truth`,
#'   `raw_edges`, `edge_features`, `model`, `model_assessment`,
#'   `predicted_edges`, and potentially other data.
#' @export
#'
#' @examples
#' # From the data folder of the CICTv2 GitHub repo, download and load
#' # SERGIO_DS4_net0_gene_expression_matrix.rda and
#' # SERGIO_DS4_net0_ground_truth.rda
#' runCICT(gene_expression_matrix = SERGIO_DS4_gene_expression_matrix,
#'         ground_truth = SERGIO_DS4_ground_truth)
runCICT <- function(gene_expression_matrix = NULL,
                    ground_truth = NULL,
                    config_path = NULL,
                    in_format = "separate",
                    ...) {
  cict_data_obj <-
    checkData(
      gene_expression_matrix = gene_expression_matrix,
      ground_truth = ground_truth,
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

    # set unnamed args based on input format
    # TODO: once other functions can take both config and input formats, should
    # delete this and change to passing all arguments directly
    unnamed_args <- list(...)
    if (in_format == "config_file") {
      # args from ellipses are ignored when using config file
      config <- yaml::yaml.load_file(config_path)
      config_names <- names(config)
      unnamed_args <-
        config[!(config_names %in% names(cict_data_obj))]
      unnamed_args$in_format <- "separate"

      # set results dir to work dir if no otherwise specified
      if(is.null(unnamed_args$results_dir)) {
        unnamed_args$results_dir = "."
      }

      # if there is an existing log file, clear it
      log_file <- file.path(unnamed_args$results_dir, "log")
      file.create(log_file)
      # redirect printed notes to log file
      sink(log_file)
    }
    args <- list(c(cict_data_obj, unnamed_args))

    # run pipeline
    print(paste("[",Sys.time(),"]","Began calculating raw edge weights",sep=" "));
    cict_data_obj$raw_edges <-
      do.call("calculateRawEdges", c(cict_data_obj, unnamed_args))$raw_edges
    print(paste("[",Sys.time(),"]","Finished calculating raw edge weights",sep=" "));

    print(paste("[",Sys.time(),"]","Began calculating edge features",sep=" "));
    cict_data_obj$edge_features <-
      do.call("prepareEdgeFeatures", c(cict_data_obj, unnamed_args))$edge_features
    print(paste("[",Sys.time(),"]","Finished calculating edge features",sep=" "));

    print(paste("[",Sys.time(),"]","Began predicting edges",sep=" "));
    pe_out <-
      do.call("predictEdges", c(cict_data_obj, unnamed_args))
    cict_data_obj$model <- pe_out$model
    cict_data_obj$model_assessment <- pe_out$model_assessment
    cict_data_obj$predicted_edges <- pe_out$predicted_edges
    print(paste("[",Sys.time(),"]","Finished predicting edges",sep=" "));

    if(in_format == "config_file") {
      write.csv(cict_data_obj$predicted_edges, file = file.path(unnamed_args$results_dir, "predicted_edges.csv"), row.names = FALSE)
      saveRDS(cict_data_obj, file = file.path(unnamed_args$results_dir, "cict_data.Rds"))

      # close connection to log file
      sink()
    }

    return(cict_data_obj)
  }, error = function(e) {
    message(e$message)
  })
  return(NULL)
}
