#NEW PEF


#prepareEdgeFeatures(raw_edges = trial_dt_edge, gene_expression_matrix = trial_dt_geneexp, cict_raw_edge_col = 'Spearman')
# trial_dt_edge <- read.csv("/Users/vvtch/Desktop/sigmafolder/inputs/rawEdges.csv")
# trial_dt_geneexp <- read.csv("/Users/vvtch/Desktop/sigmafolder/inputs/SERGIO_DS4/net1/ExpressionData.csv")
#' Prepare Edge Features
#'
#' This function prepares edge features from raw edges and a gene expression matrix.
#'
#' @param raw_edges Data frame. Raw edges data. Default is NULL.
#' @param gene_expression_matrix Data frame. Gene expression matrix. Default is NULL.
#' @param cict_raw_edge_col Character. Column name for raw edge calculation. Default is 'Spearman'.
#' @param in_format Character. Format of the input data. Can be "separate" or "data_obj". Default is "separate".
#' @param ... Additional arguments to be passed to other functions.
#' @return A list containing the prepared edge features and other relevant data.
#' @details This function processes the raw edges and gene expression matrix to prepare edge features. It supports two input formats: "separate" and "data_obj".
#' @examples
#' \dontrun{
#' # Example usage:
#' raw_edges <- read.csv("path/to/raw_edges.csv")
#' gene_expression_matrix <- read.csv("path/to/gene_expression_matrix.csv")
#' edge_features <- prepareEdgeFeatures(
#'   raw_edges = raw_edges,
#'   gene_expression_matrix = gene_expression_matrix,
#'   cict_raw_edge_col = 'Spearman'
#' )
#' }
#' @export
prepareEdgeFeatures <-
  function(raw_edges = NULL,
           gene_expression_matrix = NULL,
           cict_raw_edge_col = 'Spearman',
           in_format = "separate",
           prior = NULL, ...) {
    # TODO: allow config and throw error if in_format is not valid
    if (in_format == "separate") {
      dt_edge <- raw_edges
      dt_geneexp <- gene_expression_matrix
    }
    else if (in_format == "data_obj") {
      dt_edge <- in_data_obj$raw_edges
      dt_geneexp <- in_data_obj$gene_expression_matrix
    }
    backup_raw_edges <- dt_edge
    # define the hardcoded variables
    earlyThresholdForGraphAnalysis <- 0
    #third, make sure the data is read
    #colnames(dt_geneexp)[1] <-'gene'

    #finally, functions
    results2 <-
      calculate_f0(
        prepare_table_pef(
          dt_edge,
          dt_geneexp,
          cict_raw_edge_col,
          earlyThresholdForGraphAnalysis
        )
      )
    results3 <- calculate_f1(results2)
    results4 <- calculate_f2(results3)

    results5 <- results4
    results5$raw_edges <- tibble::as_tibble(backup_raw_edges)

    # Convert prior to a data frame if it is not already
    prior <- as.data.frame(prior)

    # Perform the left join
    results6 <- results5
    if (!is.null(prior)) {
    # Convert prior to a data frame if it is not already
    prior <- as.data.frame(prior)
    # Add a new column "prior" to edge_features
    results6$edge_features$prior <- ifelse(
    paste(results6$edge_features$src, results6$edge_features$trgt) %in%
      paste(prior$src, prior$trgt), 1, 0
  )
  } else {
    results6 <- results5
    results6$edge_features <- results5$edge_features
  }
  return(results6)
  }

prepare_table_pef <-
  function(dt_edge,
           dt_geneexp,
           cict_raw_edge_col = 'Spearman',
           earlyThresholdForGraphAnalysis) {
    #dt_edge_back00 = dt_edge # create a backup of the originаl data table
    dt_edge$edgeTyp = "" #creating new empty column edgeTyp

    dt_edge <-
      dt_edge %>% dplyr::mutate(Weight = .data[[cict_raw_edge_col]]) %>% dplyr::select(src, trgt, everything()) #make src and trg the first 2 columns
    # defining a new data table with vertices
    dt_vertices = data.table::data.table()  #change n.itm.v -> dt.vertices

    #pipeline to format dt_vertices properly

    # This code performs the following operations on the data frame `dt_geneexp`:
    # 1. Converts `dt_geneexp` to a data frame using `setDF`.
    # 2. Adds a new column `ocr` which is the sum of all columns except 'gene'.
    # 3. Replaces values in `ocr` that are 0 or NA with `PSUDEO_ZERO`.
    # 4. Selects only the 'gene' and 'ocr' columns.
    # 5. Renames the 'gene' column to 'vID'.

    dt_vertices <- data.frame(gene = rownames(dt_geneexp),
                              rowSums(dt_geneexp)) %>%
      dplyr::mutate(ocr = rowSums(.[which(!colnames(.) %in% c('gene'))])) %>%
      dplyr::mutate(ocr = ifelse(ocr == 0 |
                                   is.na(ocr), PSUDEO_ZERO, ocr)) %>%
      dplyr::select(gene, ocr)
    dt_vertices <-
      data.table::setDT(dt_vertices)[, `:=`(OcrInp = ocr, OcrOut = ocr)][, ocr := NULL]

    #convert our data about gene-gene assotiation from matrix to non directed graph
    #rm(ig)
    ig = igraph::graph_from_data_frame(dt_edge[dt_edge$Weight > earlyThresholdForGraphAnalysis,], directed =
                                 FALSE)

    return(list(dt_edge, dt_geneexp, dt_vertices, ig))
  }

calculate_f0 <- function (results) {
  pseudo_zero <- 0
  dt_edge <-
    results[[1]]  # modified dt_edge raw edges with weight values of selected type added
  dt_geneexp <-
    results[[2]] # unmodified dt_geneexp, expression data
  dt_vertices <- results[[3]] # a new dataframe with features
  ig <- results[[4]] # the graph

  # Calculate out-degree and in-degree for each vertex
  outd <- igraph::degree(ig, v = igraph::V(ig), mode = "out")
  outd <- data.frame(Outdegree = outd, subcat = names(outd))

  ind <- igraph::degree(ig, v = igraph::V(ig), mode = "in")
  ind <- data.frame(Indegree = ind, subcat = names(ind))

  # Merge in-degree and out-degree with dt_vertices
  dt_vertices <- dt_vertices %>%
    dplyr::left_join(ind, by = c("gene" = "subcat")) %>%
    dplyr::left_join(outd, by = c("gene" = "subcat")) %>%
    dplyr::mutate(
      Indegree = ifelse(is.na(Indegree), pseudo_zero, Indegree),
      Outdegree = ifelse(is.na(Outdegree), pseudo_zero, Outdegree)
    )

  # Calculate probabilities
  vertexOcrSumIn <- sum(dt_vertices$OcrInp, na.rm = TRUE)
  vertexOcrSumOut <- sum(dt_vertices$OcrOut, na.rm = TRUE)
  dt_vertices <- dt_vertices %>%
    dplyr::mutate(probInp = OcrInp / vertexOcrSumIn,
           probOut = OcrOut / vertexOcrSumOut)

  vinfcols <-
    c("gene",
      "OcrInp",
      "OcrOut",
      "probInp",
      "probOut",
      "Indegree",
      "Outdegree")

  # Merge dt_edge with dt_vertices
  dt_edge <- dt_edge %>%
    dplyr::left_join(dt_vertices[, ..vinfcols], by = c("src" = "gene")) %>%
    dplyr::left_join(dt_vertices[, ..vinfcols], by = c("trgt" = "gene"))

  # Summarize OcrInp and OcrOut
  dt_vertices_trgtsum <- dt_vertices %>%
    dplyr::left_join(dt_edge, by = c("gene" = "src")) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      SumOcrInp.y = sum(OcrInp.y, na.rm = TRUE),
      SumOcrOut.y = sum(OcrOut.y, na.rm = TRUE)
    )

  dt_vertices_srcsum <- dt_vertices %>%
    dplyr::left_join(dt_edge, by = c("gene" = "trgt")) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      SumOcrInp.x = sum(OcrInp.x, na.rm = TRUE),
      SumOcrOut.x = sum(OcrOut.x, na.rm = TRUE)
    )

  # Merge summarized data back to dt_vertices
  dt_vertices <- dt_vertices %>%
    dplyr::left_join(dt_vertices_trgtsum, by = "gene") %>%
    dplyr::inner_join(dt_vertices_srcsum, by = "gene")

  cols <-
    c("SumOcrInp.x",
      "SumOcrOut.x",
      "SumOcrInp.y",
      "SumOcrOut.y",
      "OcrOut",
      "OcrInp")

  # Convert factors to numeric
  factorToNumeric <- function(x)
    as.numeric(as.character(x))
  dt_vertices_num <-
    data.table::setDT(dt_vertices)[, lapply(.SD, factorToNumeric), .SDcols = cols]
  dt_vertices <-
    dt_vertices[, !names(dt_vertices) %in% cols, with = FALSE]
  dt_vertices <- cbind(dt_vertices_num, dt_vertices)
  rm(dt_vertices_num)
  gc()

  # Replace NA values
  data.table::setDT(dt_vertices)
  dt_vertices[is.na(SumOcrOut.x), SumOcrOut.x := 1]
  dt_vertices[is.na(SumOcrInp.y), SumOcrInp.y := 1]

  # Merge dt_edge with dt_vertices again
  vinfcols <-
    c("gene",
      "SumOcrInp.x",
      "SumOcrOut.x",
      "SumOcrInp.y",
      "SumOcrOut.y")
  dt_edge <- dt_edge %>%
    dplyr::inner_join(dt_vertices[, ..vinfcols], by = c("src" = "gene")) %>%
    dplyr::inner_join(dt_vertices[, ..vinfcols], by = c("trgt" = "gene"))

  # Replace NA values in dt_edge
  dt_edge <-
    replace_na(
      dt_edge,
      replace = list(
        SumOcrOut.x.x = 0,
        SumOcrOut.x.y = 0,
        SumOcrInp.x.x = 0,
        SumOcrInp.x.y = 0,
        SumOcrInp.y.x = 0,
        SumOcrInp.y.y = 0,
        SumOcrOut.y.x = 0,
        SumOcrOut.y.y = 0
      )
    )

  # Calculate confidence and contribution
  dt_edge <- dt_edge %>%
    dplyr::mutate(conf = Weight / OcrOut.x, contrib = Weight / OcrInp.y) %>%
    dplyr::mutate(
      conf = ifelse(is.na(conf) | is.infinite(conf), 0, conf),
      contrib = ifelse(is.na(contrib) |
                         is.infinite(contrib), 0, contrib)
    )

  # Discretize confidence and contribution
  nbins <- as.integer(sqrt(ncol(dt_geneexp)) * 1)
  dt_edge <- dt_edge %>%
    dplyr::mutate(
      confdisc = unlist(infotheo::discretize(conf, "equalfreq", nbins)),
      contribdisc = unlist(infotheo::discretize(contrib, "equalfreq", nbins))
    )
  e_cnfcnt <- dt_edge %>% dplyr::select(src, trgt, confdisc, contribdisc)

  # Remove unnecessary columns
  data.table::setDT(e_cnfcnt)
  data.table::setkeyv(e_cnfcnt, c('src', 'trgt'))

  return(list(dt_edge, dt_geneexp, dt_vertices, ig))
}

myMean <- function (x,
                    default = 0,
                    na.rm = TRUE,
                    ...) {
  result = mean(x, na.rm = na.rm, ...)
  if (is.na(result))
    default
  else
    result
}

myMedian <- function (x,
                      default = 0,
                      na.rm = TRUE,
                      ...) {
  result = median(x, na.rm = na.rm, ...)
  if (is.na(result))
    default
  else
    result
}

mySD <- function (x, default = 0, na.rm = TRUE) {
  result = sd(x, na.rm)
  if (is.na(result))
    default
  else
    result
}

mySkewness <- function (x,
                        default = 0,
                        na.rm = TRUE,
                        ...) {
  result = moments::skewness(x, na.rm, ...)
  if (is.na(result))
    default
  else
    result
}

myKurtosis <- function (x,
                        default = 3,
                        na.rm = TRUE,
                        ...) {
  result = PerformanceAnalytics::kurtosis(x, na.rm, ...)
  if (is.na(result))
    default
  else
    result
}

myMedianAbsoluteDeviation <-
  function (x,
            center = median(x),
            constant = 1.4826,
            na.rm = FALSE,
            low = FALSE,
            high = FALSE,
            default = 0) {
    result = mad(x, center , constant , na.rm, low , high)
    if (is.na(result))
      default
    else
      result
  }

myMin = function(x, lowwerLimit = -Inf) {
  minim = Inf

  t = lapply(x, function(x)
    if (is.finite(x) & x < minim)
      x
    else
      minim)
  mn = min(unlist(t), na.rm = TRUE)
  if (lowwerLimit == -Inf)
    mn
  else
    max(mn, lowwerLimit)
}

#' Replace NA values in a DataFrame based on a replacement list
#'
#' This function replaces NA values in a DataFrame with specified values based on a pattern matching of column names.
#'
#' @param df A DataFrame in which NA values need to be replaced.
#' @param rplist A named list where the names are patterns to match column names in the DataFrame, and the values are the replacement values for NA.
#'
#' @return A DataFrame with NA values replaced according to the specified patterns and replacement values.
#'
#' @examples
#' df <- data.frame(a1 = c(NA, 2, 3), a2 = c(4, NA, 6), b1 = c(7, 8, NA))
#' rplist <- list(a = 0, b = 1)
#' my_replace_na(df, rplist)
#'
#' @importFrom tidyr replace_na
#' @noRd
my_replace_na <- function(df, rplist)
{
  rplptrns = names(rplist)
  rplvals = unlist(rplist)

  dfcols = colnames(df)
  rplcols = lapply(rplptrns, function(x)
    dfcols[stringr::str_detect(dfcols, paste0(".*", x))]) #pattern checking

  rplist2 = list()
  for (i in 1:length(rplist)) {
    l = unlist(rplcols[i])
    v = rplvals[i]
    l1 = rep(v, length(l))
    names(l1) <- l
    rplist2 = append(rplist2, l1)
  }

  df = replace_na(df, replace = rplist2)
  df
}


calculate_moments <- function(data, group_col, value_col, prefix, tn) {
  lmoments <-
    data %>% dplyr::group_by(!!dplyr::sym(group_col)) %>% dplyr::do(extractLmoments(.[[paste0(value_col, "N")]]))
  summary_stats <-
    data %>% dplyr::group_by(!!dplyr::sym(group_col)) %>% dplyr::summarise(
      Total = sum(!!dplyr::sym(value_col), na.rm = TRUE),
      Mean = myMean(!!dplyr::sym(value_col), 0, na.rm = TRUE),
      Median = myMedian(!!dplyr::sym(value_col), 0, na.rm = TRUE),
      SD = mySD(!!dplyr::sym(value_col), 0, na.rm = TRUE),
      Skew = mySkewness(!!dplyr::sym(value_col), 0, na.rm = TRUE),
      Kurt = myKurtosis(!!dplyr::sym(value_col), 3, na.rm = TRUE),
      MADconst = ifelse(is.na(sn::qsc(.75, Mean, SD, Skew)), 1.488, sn::qsc(.75, Mean, SD, Skew)),
      MAD = myMedianAbsoluteDeviation(
        !!dplyr::sym(value_col),
        Median,
        MADconst,
        na.rm = TRUE,
        default = 0
      ),
      NTotal = sum(!!dplyr::sym(paste0(value_col, "N")), na.rm = TRUE),
      NMean = myMean(!!dplyr::sym(paste0(value_col, "N")), 0, na.rm = TRUE),
      NMedian = myMedian(!!dplyr::sym(paste0(value_col, "N")), 0, na.rm = TRUE),
      NSD = mySD(!!dplyr::sym(paste0(value_col, "N")), na.rm = TRUE, 0),
      NSkew = mySkewness(!!dplyr::sym(paste0(value_col, "N")), 0, na.rm = TRUE),
      NKurt = myKurtosis(!!dplyr::sym(paste0(value_col, "N")), 3, na.rm = TRUE),
      NMADconst = ifelse(is.na(sn::qsc(
        .75, NMean, NSD, NSkew
      )), 1.488, sn::qsc(.75, NMean, NSD, NSkew)),
      NMAD = myMedianAbsoluteDeviation(
        !!dplyr::sym(paste0(value_col, "N")),
        NMedian,
        NMADconst,
        na.rm = TRUE,
        default = 0
      ),
      tnTotal=sum(tn,na.rm= TRUE),
       tnMean=myMean(tn, 0,na.rm= TRUE),
       tnMedian=myMedian(tn, 0,na.rm= TRUE),
       tnSD=mySD(tn, na.rm = TRUE,0),
       tnSkew=mySkewness(tn, 0,na.rm =TRUE),
       tnKurt=myKurtosis(tn, 3,na.rm = TRUE),
       tnMADconst=ifelse(is.na(sn::qsc(.75,tnMean,tnSD,tnSkew)),1.488,sn::qsc(.75,tnMean,tnSD,tnSkew)),
       tnMAD= myMedianAbsoluteDeviation(tn,tnMedian,tnMADconst,na.rm=TRUE,default=0)
    )
  full_stats <-
    summary_stats %>% dplyr::inner_join(lmoments, by = group_col)
  names(full_stats) <- paste0(prefix, names(full_stats))
  return(full_stats)
}

calculate_f1 <- function(results2) {
  #same results unpacking process as in calculate_f0

  dt.edge <-
    results2[[1]]  # modified dt.edge raw edges with weight values of selected type added
  dt.geneexp <-
    results2[[2]] # unmodified dt_geneexp, expression data
  dt.vertices <- results2[[3]] # a new dataframe with features
  ig <- results2[[4]] # the graph

  pseudo_zero <- 0
  PSUDEO_ZERO <- 0
  emissionRate <- 1
  ngens = nrow(dt.geneexp)
  Nobservations = ncol(dt.geneexp) - 1
  nbins = as.integer(sqrt(Nobservations) * 1)
  breaks = 1:nbins
  breaks = NULL

  # This script processes vertices and their associated edges by applying discretization
  # and histogram generation. The main steps include:
  # 1. Initialization of variables.
  # 2. Iteration over a range of generations (ngens).
  # 3. For each generation, the script:
  #    a. Retrieves the vertex ID (gene) of the first vertex.
  #    b. Applies discretization to the edges' confidence and contribution values.
  #    c. Generates histograms for the discretized values.
  #    d. Updates the vertices data table with the generated histograms.
  # 4. Error handling is implemented to catch and count any errors that occur during processing.
  thisrow = 0
  j = 1
  i = 1
  errCnt = 0
  for (j in (1):ngens) {
    tryCatch({
      e.a = e.b = e.ab = NULL
      firstg = dt.vertices[j,]$gene

      if (!globalDistDiscretization) {
        #Discretization was applied on edges of each particular node
        a.confdisc = myhist(infotheo::discretize(e.cnfcnt[src == firstg,]$conf, "equalwidth", nbins)$X,
                            nbins,
                            breaks)
        a.contribdisc =  myhist(
          infotheo::discretize(e.cnfcnt[src == firstg,]$contrib, "equalwidth", nbins)$X,
          nbins,
          breaks
        )

      } else if (globalDistDiscretization) {
        #Discretization was applied on all edges
        suppressMessages({
          a.confhist = myhist(e.cnfcnt[src == firstg,]$confdisc,
                              nbins,
                              breaks,
                              plot = F,
                              prob = F)
          a.contribhist =  myhist(e.cnfcnt[src == firstg,]$contribdisc,
                                  nbins,
                                  breaks,
                                  plot = F,
                                  prob = F)
        })
      }
      dt.vertices[gene == firstg, `:=`(confhist = list(a.confhist))]
      dt.vertices[gene == firstg, `:=`(contribhist = list(a.contribhist))]
    }, error = function(e) {
      errCnt <<- errCnt + 1
      sprintf("Error in row: %s, first: %s, " , j, firstg)
    })
  }

  # This script prints a message indicating the number of nodes that didn't have edges for
  # confhist or contribhist calculations, and then removes several data tables from the
  # environment to free up memory.
  #
  # Variables:
  # - errCnt: Integer, the count of nodes without edges for confhist or contribhist calculations.
  #

  # This script performs several data transformations and calculations on the 'dt.edge' dataframe.
  #
  # The following operations are performed:
  # 1. Two new columns are created:
  #    - 'srctrgtSum': Sum of 'OcrOut.x' and 'OcrInp.y'.
  #    - 'srctrgtProduct': Product of 'OcrOut.x' and 'OcrInp.y'.
  # 2. The 'srctrgtSum' and 'srctrgtProduct' columns are updated to replace NA values or zeros with 'PSUDEO_ZERO'.
  # 3. Two additional columns are calculated:
  #    - 'HTR': Harmonized Transition Rate, calculated as (Weight * srctrgtSum) / srctrgtProduct.
  #    - 'EE': Some efficiency measure, calculated as (OcrOut.x^2 * OcrInp.y) / srctrgtSum.
  # 4. A subset of 'dt.edge' is created, containing rows where 'HTR' is NA.

  dt.edge = dt.edge %>% dplyr::mutate(srctrgtSum = OcrOut.x + OcrInp.y,
                                      srctrgtProduct = OcrOut.x * OcrInp.y) %>%
    dplyr::mutate(
      srctrgtSum = ifelse(is.na(srctrgtSum) |
                            srctrgtSum == 0, PSUDEO_ZERO, srctrgtSum),
      srctrgtProduct = ifelse(
        is.na(srctrgtProduct) |
          srctrgtProduct == 0 ,
        PSUDEO_ZERO,
        srctrgtProduct
      )
    ) %>%
    dplyr::mutate(
      HTR = (Weight * srctrgtSum) / srctrgtProduct,
      EE = OcrOut.x ^ 2 * OcrInp.y / srctrgtSum
    )

  b = dt.edge %>% dplyr::filter(is.na(HTR))

  # This code snippet performs the following operations:
  # 1. Groups the 'dt.edge' data frame by the 'src' column and calculates the sum of 'HTR' and 'EE' columns for each group, ignoring NA values.
  #    The results are stored in a new data frame 'self' with columns 'scfSumHTR' and 'scfSumEE'.
  # 2. Groups the 'dt.edge' data frame by the 'trgt' column and calculates the sum of 'HTR' and 'EE' columns for each group, ignoring NA values.
  #    The results are stored in a new data frame 'others' with columns 'ocbSumHTR' and 'ocbSumEE'.
  self = dt.edge %>% dplyr::group_by(src) %>% dplyr::summarise(scfSumHTR = sum(HTR, na.rm = TRUE),
                                                 scfSumEE = sum(EE, na.rm = TRUE))
  others = dt.edge %>% dplyr::group_by(trgt) %>%  dplyr::summarise(ocbSumHTR = sum(HTR, na.rm = TRUE),
                                                     ocbSumEE = sum(EE, na.rm = TRUE))


  # This script performs the following operations:

  # 1. Defines a vector `vinfcols` containing column names related to vertex information.
  # 2. Converts the data tables `dt.edge` and `dt.vertices` to data frames using `setDF`.
  # 3. Joins the `dt.edge` data frame with the `dt.vertices` data frame twice:
  #    a. First join is performed on the "src" column of `dt.edge` with the "gene" column of `dt.vertices`.
  #    b. Second join is performed on the "trgt" column of `dt.edge` with the "gene" column of `dt.vertices`.
  # 1. Creates a backup of the original 'dt.vertices' data frame.
  # 2. Joins 'dt.vertices' with the 'self' data frame on the 'gene' and 'src' columns.
  # 3. Joins 'dt.vertices' with the 'others' data frame on the 'gene' and 'trgt' columns.
  # 4. Calculates the number of NA values in each column of 'dt.vertices' and prints the columns with NA values.
  dt.vertices.back = dt.vertices
  dt.vertices = dt.vertices %>% dplyr::left_join(self, by = c("gene" = "src"))
  dt.vertices = dt.vertices %>% dplyr::left_join(others, by = c("gene" = "trgt"))

  # This script processes the 'dt.vertices' data table by performing the following steps:
  # 1. Identifies columns with missing values and displays the count of missing values for each column.
  # 2. Converts 'dt.vertices' to a data frame and replaces NA values in specified columns with default values:
  #    - 'scfSumHTR', 'ocbSumHTR', 'scfSumEE', 'ocbSumEE' are replaced with 0.
  #    - 'DiagnosesAb' is replaced with an empty string.
  #    - 'SumOcrInp.x' is replaced with 0.
  # 3. Creates a new data table 'dt.vertices1' containing a subset of columns from 'dt.vertices':
  #    - 'gene', 'scfSumHTR', 'ocbSumHTR', 'scfSumEE', 'ocbSumEE'.
  a = sapply(dt.vertices, function(x)
    sum(is.na(x)))
  a[a > 0]
  #View(dt.vertices[,c("gene","scfSumHTR","ocbSumHTR"),with=FALSE])
  data.table::setDF(dt.vertices)
  dt.vertices = tidyr::replace_na(
    dt.vertices,
    replace = list(
      scfSumHTR = 0,
      ocbSumHTR = 0,
      scfSumEE =
        0,
      ocbSumEE = 0,
      DiagnosesAb = "",
      SumOcrInp.x = 0
    )
  )

  dt.vertices1 = dt.vertices[, c("gene", "scfSumHTR", "ocbSumHTR", "scfSumEE", "ocbSumEE")]

  # This code converts the data tables 'dt.edge' and 'dt.vertices1' to data frames.
  # It then performs two left joins on 'dt.edge':
  # 1. Joins 'dt.edge' with 'dt.vertices1' on the 'src' column of 'dt.edge' and 'gene' column of 'dt.vertices1'.
  # 2. Joins 'dt.edge' with 'dt.vertices1' on the 'trgt' column of 'dt.edge' and 'gene' column of 'dt.vertices1'.
  data.table::setDF(dt.edge)
  data.table::setDF(dt.vertices1)
  dt.edge = dt.edge %>% dplyr::left_join(dt.vertices1, by = c("src" = "gene"))
  dt.edge = dt.edge %>% dplyr::left_join(dt.vertices1, by = c("trgt" = "gene"))

  #   Enhance edges conf contrib -----
  # This script processes a data frame `dt.edge` by performing various calculations and transformations.
  #
  # Preconditions:
  # - The constant `PSUDEO_ZERO` must be defined before running this script.
  #
  # The script performs the following operations:
  # 1. Converts `dt.edge` to a data frame.
  # 2. Mutates the data frame by adding several new columns:
  #    - `V`: Absolute difference between `OcrOut.x` and `OcrInp.y`.
  #    - `R`: Ratio of `V` to `Weight`.
  #    - `P`: Product of `V` and `Weight`.
  #    - `NHTRfromSource`: Normalized HTR from source.
  #    - `NHTRtoTarget`: Normalized HTR to target.
  #    - `tNHTR`: Sum of `NHTRfromSource` and `NHTRtoTarget`.
  #    - `NEEfromSource`: Normalized EE from source.
  #    - `NEEtoTarget`: Normalized EE to target.
  #    - `NEE_NHTR_Diff`: Difference between `NHTRfromSource` and `NEEfromSource`.
  #    - `NEE_NHTR_Ratio`: Ratio of `NHTRfromSource` to `NEEfromSource`, using `PSUDEO_ZERO` to avoid division by zero.
  #    - `confN`: Normalized confidence.
  #    - `contribN`: Normalized contribution.
  #    - `tn`: Sum of `confN` and `contribN`.
  #    - `t`: Sum of `conf` and `contrib`.
  #    - `Pout`: Weighted difference between `OcrOut.x` and `OcrInp.y` normalized by `OcrOut.x`.
  #    - `Pin`: Weighted difference between `OcrOut.x` and `OcrInp.y` normalized by `OcrInp.y`.
  #    - `slope`: Absolute difference between `OcrOut.x` and `OcrInp.y` digeneed by `Weight`.
  #    - `hypotenuse`: Hypotenuse calculated from `V` and `Weight`.
  #    - `EO`: Emission rate normalized by `SumOcrOut.x.x`.
  #    - `EI`: Emission rate normalized by `SumOcrOut.x.y`.
  #    - `OEER`: Weight digeneed by `EO`.
  #    - `OERR`: Weight digeneed by `EI`.
  #    - `toe`: Sum of `OEER` and `OERR`.
  #    - `AM`: Arithmetic mean of `OEER` and `OERR`.
  #    - `GM`: Geometric mean of `OEER` and `OERR`.
  #    - `HM`: Harmonic mean of `OEER` and `OERR`.
  #    - `UR`: Unexpected resistance calculated as `V` digeneed by the difference between `Weight` and `EO`.
  #    - `URR`: Unexpected resistance calculated as `V` digeneed by the difference between `Weight` and `EI`.
  if (!exists("PSUDEO_ZERO"))
    stop("define PSUDEO_ZERO")
  dt.edge = data.table::setDF(dt.edge) %>% #dplyr::mutate(Weight=mf.mi) %>%
    dplyr::mutate(
      V = abs(OcrOut.x - OcrInp.y),
      R = V / Weight,
      P = V * Weight,

      NHTRfromSource = HTR / scfSumHTR.x,
      NHTRtoTarget = HTR / ocbSumHTR.y,
      tNHTR = NHTRfromSource + NHTRtoTarget,

      NEEfromSource = EE / scfSumEE.x,
      NEEtoTarget = EE / ocbSumEE.y,

      NEE_NHTR_Diff = NHTRfromSource - NEEfromSource,
      NEE_NHTR_Ratio = NHTRfromSource / ifelse(NEEfromSource ==
                                                 0, PSUDEO_ZERO, NEEfromSource),

      confN = conf * (OcrInp.y / srctrgtSum),
      #NHTRfromSource
      contribN = contrib * (OcrOut.x / srctrgtSum),
      #NHTRtoTarget, #
      tn = confN + contribN,
      t = conf + contrib,
      #+1 is necessary to protect recursive edges from going infinity

      Pout = (OcrOut.x - OcrInp.y) * Weight / OcrOut.x,
      Pin = (OcrOut.x - OcrInp.y) * Weight / OcrInp.y,
      #/maxItmstOcr
      slope = abs(OcrOut.x - OcrInp.y) / Weight,
      #tng
      hypotenuse = sqrt(V ^ 2 + Weight ^ 2),
      EO = emissionRate * (srctrgtProduct) / (SumOcrOut.x.x +
                                                1),
      EI = emissionRate * (srctrgtProduct) / (SumOcrOut.x.y +
                                                1),
      OEER = Weight / EO,
      OERR = Weight / EI,
      toe = OEER + OERR,
      AM = (OEER + OERR) / 2,
      GM = (OEER * OERR) ^ .5,
      HM = 2 / ((1 / OEER) + (1 / OERR)),
      UR = V / (Weight - EO) ,
      URR = V / (Weight - EI)
    ) #unexpected resistance

  # This script performs data quality checks on the 'dt.edge' data frame.
  # Specifically, it calculates the following:
  # 1. The number of NA values in each column of 'dt.edge' and filters out columns with no NA values.
  # 2. The number of zero values (excluding NAs) in each column of 'dt.edge' and filters out columns with no zero values.
  # 3. The number of infinite values in each column of 'dt.edge' and filters out columns with no infinite values.
  #
  # Note: The script contains commented-out code for additional data processing steps, such as:
  # - Filtering rows with infinite values in the 'contribZ' column and viewing specific columns.
  # - Filtering rows where 'Weight' is zero and viewing specific columns.
  a = sapply(dt.edge, function(x)
    sum(is.na(x)))
  a[a > 0]
  b = sapply(dt.edge, function(x)
    sum(x == 0 & !is.na(x)))
  b[b > 0]
  c = sapply(dt.edge, function(x)
    sum(is.infinite(x)))
  c[c > 0]

  #   Add Power parameters to vertices and edges ----
  #all the power input from different sources to each trgt
  # This script calculates summary statistics for power parameters in a dataset.
  #
  # The dataset `dt.edge` is grouped by the target (`trgt`) and source (`src`) columns.
  # For each group, the following summary statistics are calculated:
  #
  # - `Pinsum`: Sum of `Pin` values, ignoring missing values.
  # - `PinAbsSum`: Sum of the absolute values of `Pin`, ignoring missing values.
  # - `PinSD`: Standard deviation of `Pin`, ignoring missing values.
  # - `PinMean`: Mean of `Pin`, ignoring missing values.
  #
  # Similarly, for the source (`src`) groups, the following summary statistics are calculated:
  #
  # - `Poutsum`: Sum of `Pout` values, ignoring missing values.
  # - `PoutAbsSum`: Sum of the absolute values of `Pout`, ignoring missing values.
  # - `PoutSD`: Standard deviation of `Pout`, ignoring missing values.
  # - `PoutMean`: Mean of `Pout`, ignoring missing values.
  powerParamsIn = dt.edge %>% dplyr::group_by(trgt) %>%
    dplyr::summarise(
      Pinsum = sum(Pin, na.rm = TRUE),
      PinAbsSum = sum(abs(Pin), na.rm = TRUE),
      PinSD = sd(Pin, na.rm = TRUE),
      PinMean = mean(Pin, na.rm = TRUE)
    )
  powerParamsOut = dt.edge %>% dplyr::group_by(src) %>%
    dplyr::summarise(
      Poutsum = sum(Pout, na.rm = TRUE),
      PoutAbsSum = sum(abs(Pout), na.rm = TRUE),
      PoutSD = sd(Pout, na.rm = TRUE),
      PoutMean = mean(Pout, na.rm = TRUE)
    )

  # This script checks for missing values (NA) in the elements of two lists: powerParamsOut and powerParamsIn.
  # It uses the sapply function to apply a function over each element of the lists.
  # The function calculates the sum of NA values in each element.
  # The results are stored in variable 'a' and only the elements with more than 0 NA values are displayed.
  a = sapply(powerParamsOut, function(x)
    sum(is.na(x)))
  a[a > 0]
  a = sapply(powerParamsIn, function(x)
    sum(is.na(x)))
  a[a > 0]

  # This code performs the following operations:

  # 1. Joins the 'dt.vertices' data frame with 'powerParamsIn' data frame on the 'gene' column from 'dt.vertices' and 'trgt' column from 'powerParamsIn'.
  # 2. Joins the resulting 'dt.vertices' data frame with 'powerParamsOut' data frame on the 'gene' column from 'dt.vertices' and 'src' column from 'powerParamsOut'.
  dt.vertices = dt.vertices %>% dplyr::left_join(powerParamsIn, by = c("gene" =
                                                                  "trgt"))
  dt.vertices = dt.vertices %>% dplyr::left_join(powerParamsOut, by = c("gene" =
                                                                   "src"))

  # This line of code replaces NA (missing) values in the 'dt.vertices' data frame.
  # Specifically, it sets missing values in the 'PinSD' and 'PoutSD' columns to 0.
  # The 'replace_na' function from the 'tidyr' package is used for this purpose.
  dt.vertices = replace_na(dt.vertices, replace = list(PinSD = 0, PoutSD =
                                                         0))

  vinfcols = c(
    "gene",
    "Pinsum",
    "PinAbsSum",
    "PinSD",
    "PinMean",
    "Poutsum",
    "PoutAbsSum",
    "PoutSD",
    "PoutMean"
  )
  data.table::setDF(dt.edge)
  data.table::setDF(dt.vertices)
  # This section of the code performs inner joins on the dt.edge dataframe with the dt.vertices dataframe.
  # It first joins dt.edge with dt.vertices based on the "src" column matching the "gene" column in dt.vertices.
  # Then, it joins dt.edge with dt.vertices again based on the "trgt" column matching the "gene" column in dt.vertices.
  # The vinfcols variable specifies the columns to be selected from dt.vertices for the join operations.
  dt.edge = dt.edge %>% dplyr::inner_join(dt.vertices[, vinfcols], by = c("src" =
                                                                     "gene"))
  dt.edge = dt.edge %>% dplyr::inner_join(dt.vertices[, vinfcols], by = c("trgt" =
                                                                     "gene"))


  # This script appears to be a work in progress for handling data tables related to edges and vertices.
  #
  # The following operations are performed:
  # - A backup of the 'dt.edge' data table is created and stored in 'dt.edge.back'.
  # - A backup of the 'dt.vertices' data table is created and stored in 'dt.vertices.back'.
  #
  # Note:
  # - The commented lines suggest that there might be an intention to restore the original data tables
  #   from their backups at some point.
  dt.edge.back = dt.edge #dt.edge=dt.edge.back
  dt.vertices.back = dt.vertices
  #dt.vertices=dt.vertices.back;dt.edge=dt.edge.back

  gc()


  # This script processes a data frame `dt.edge` by adding BA factors and combined measures.
  #
  # Steps:
  # 1. Convert `dt.edge` to a data frame.
  # 2. Create a backup of `dt.edge`.
  # 3. Rename columns in `dt.edge` to include BA suffixes and select specific columns.
  # 4. Join `dt.edge` with the renamed version on `src` and `trgt` columns.
  # 5. Identify and handle missing and zero values in the joined data frame.
  # 6. Replace NA values in specific columns with predefined constants.
  # 7. Calculate differences and ratios between original and BA columns.
  # 8. Identify and handle infinite values in the final data frame.
  #
  # Key Variables:
  # - `PSUDEO_ZERO_1`, `PSUDEO_ZERO_2`: Constants used for replacing zero values.
  # - `tBAreplacement`, `tnBAreplacement`, `toeBAreplacement`, `confNBAreplacement`, `contribNBAreplacement`: Minimum values used for replacement.
  #
  # Key Functions:
  # - `setDF()`: Converts the input to a data frame.
  # - `plyr::rename()`: Renames columns in the data frame.
  # - `dplyr::select()`: Selects specific columns from the data frame.
  # - `left_join()`: Joins two data frames.
  # - `replace_na()`: Replaces NA values in the data frame.
  # - `dplyr::mutate()`: Adds new columns or modifies existing ones in the data frame.
  #
  # Notes:
  # - The script includes commented-out sections for potential future use.
  # - The script uses `gc()` to trigger garbage collection and free up memory.
  #    Add BA factors and combined measures -----
  if (TRUE) {
    data.table::setDF(dt.edge)
    gc()
    dt.edge.back = dt.edge
    dt.edge1 = dt.edge
    dt.edge1 = plyr::rename(
      dt.edge,
      replace = c(
        'Weight' = 'IBA' ,
        'UR' = 'URBA' ,
        'AM' = 'AMBA' ,
        'OEER' = 'OEERBA' ,
        'OERR' = 'OERRBA' ,
        'conf' = 'confBA' ,
        'contrib' = 'contribBA',
        'confN' = 'confNBA' ,
        'contribN' = 'contribNBA' ,
        't' = 'tBA' ,
        'tn' = 'tnBA' ,
        'toe' = 'toeBA' ,
        'HTR' = 'HTRBA',
        'tNHTR' = 'tNHTRBA' ,
        'Pin' = 'PinBA',
        'Pout' = 'PoutBA',
        'NHTRfromSource' = 'NHTRfromSourceBA',
        'NHTRtoTarget' = 'NHTRtoTargetBA' ,
        'EE' = 'EEBA' ,
        'NEEfromSource' = 'NEEfromSourceBA' ,
        'NEEtoTarget' = 'NEEtoTargetBA'
      ),
      warn_missing = F
    )   %>%

      dplyr::select(any_of(
        c (
          'src',
          'trgt',
          'IBA',
          'URBA',
          'AMBA',
          'OEERBA',
          'OERRBA',
          'confBA',
          'contribBA',
          'confNBA',
          'contribNBA',
          'tBA',
          'tnBA',
          'toeBA',
          'HTRBA',
          'tNHTRBA',
          'PinBA',
          'PoutBA',
          'NHTRfromSourceBA',
          'NHTRtoTargetBA',
          'EEBA',
          'NEEfromSourceBA',
          'NEEtoTargetBA'
        )
      ))

    dt.edge2 = dt.edge %>% dplyr::left_join(dt.edge1, by = c("src" = "trgt", "trgt" =
                                                        "src"))
    rm(dt.edge1)
    a = sapply(dt.edge2, function(x)
      sum(is.na(x)))
    a[a > 0]
    b = sapply(dt.edge2, function(x)
      sum(x == 0 & !is.na(x)))
    b[b > 0]

    #Bounds finding a minimum to a specific boundary
    PSUDEO_ZERO_1 = 1e-5
    tBAreplacement = min(myMin(dt.edge2$tBA), PSUDEO_ZERO_1)
    if (tBAreplacement == 0)
      tBAreplacement = PSUDEO_ZERO_1
    tnBAreplacement = min(myMin(dt.edge2$tnBA), PSUDEO_ZERO_1)
    if (tnBAreplacement == 0)
      tnBAreplacement = PSUDEO_ZERO_1
    toeBAreplacement = min(myMin(dt.edge2$toeBA), PSUDEO_ZERO_1)
    if (toeBAreplacement == 0)
      toeBAreplacement = PSUDEO_ZERO_1
    confNBAreplacement = min(myMin(dt.edge2$confNBA), PSUDEO_ZERO_1)
    if (confNBAreplacement == 0)
      confNBAreplacement = PSUDEO_ZERO_1
    contribNBAreplacement = min(myMin(dt.edge2$contribNBA), PSUDEO_ZERO_1)
    if (contribNBAreplacement == 0)
      contribNBAreplacement = PSUDEO_ZERO_1
    tnBAreplacement = min(myMin(dt.edge2$tnBA), PSUDEO_ZERO_1)
    if (tnBAreplacement == 0)
      tnBAreplacement = PSUDEO_ZERO_1
    tBAreplacement = min(myMin(dt.edge2$tBA), PSUDEO_ZERO_1)
    if (tBAreplacement == 0)
      tBAreplacement = PSUDEO_ZERO_1
    toereplacement = min(myMin(dt.edge2$toeBA), PSUDEO_ZERO_1)
    if (toeBAreplacement == 0)
      toeBAreplacement = PSUDEO_ZERO_1

    PSUDEO_ZERO_2 = 0.0001

    dt.edge2 = replace_na(
      dt.edge2,
      replace =
        list(
          IBA = PSUDEO_ZERO_2,
          URBA = PSUDEO_ZERO_2,
          AMBA = 0,
          OEERBA = 0,
          OERRBA = 0,
          confBA = 0,
          contribBA = 0,
          confNBA = 0,
          contribNBA = 0,
          toeBA = toeBAreplacement,
          PinBA = 0,
          PoutBA = 0,
          HTRBA = 0,
          tBA = tBAreplacement,
          tnBA = tnBAreplacement,
          EEBA = 0,
          NEEfromSourceBA = 0,
          NEEtoTargetBA = 0,
          NHTRfromSourceBA = 0,
          NHTRtoTargetBA = 0,
          tNHTRBA = 0 ,
          HTRBA = 0,
          Pinsum.x = 0 ,
          PinAbsSum.x = 0,
          PinMean.x = 0,
          Poutsum.y = 0,
          PoutAbsSum.y = 0,
          PoutMean.y = 0
        )
    )

    a = sapply(dt.edge2, function(x)
      sum(is.na(x)))
    a[a > 0]
    rm(dt.edge)
    gc()
    dt.edge = dt.edge2 %>% dplyr::mutate(
      confDiff = conf - confBA,
      contribDiff = contrib - contribBA,
      confNDiff = confN - confNBA,
      contribNDiff = contribN - contribNBA,
      directionLogRatio = log((Weight /
                                 OcrOut.x) / (IBA / OcrInp.y)),
      Idiff = Weight - IBA,
      IRatio = Weight / IBA,
      URdiff = UR - URBA,
      AMdiff = AM - AMBA,
      tdiff = t - tBA,
      tRatio = t / tBA,
      tndiff = tn - tnBA,
      tnRatio = tn / tnBA,
      toediff = toe - toeBA,
      toeRatio = toe / toeBA,
      #observed to expected emission  minus observed to expected reception for each node
      deltaAttitude = OEER - OERRBA,
      Tendency = ifelse(is.na(UR) |
                          UR == 0, 0, 1 / UR) - ifelse(is.na(URBA) |
                                                         URBA == 0, 0, 1 / URBA)
    )

    a = sapply(dt.edge, function(x)
      sum(is.na(x)))
    a[a > 0]
    c = sapply(dt.edge, function(x)
      sum(is.infinite(x)))
    c[c > 0]
    b = sapply(dt.edge, function(x)
      sum(x == 0 & !is.na(x)))
    b[b > 0]
  }
  #remove the data tables dt.edge1, dt.edge2, and dt.vertices1 from the environment.
  rm(dt.edge2, dt.vertices1)

  # 1. Converts 'dt.edge' to a data frame.
  # 2. Uses 'anti_join' to find rows in 'dt.edge' where 'src' does not match any 'gene' in 'dt.vertices'.
  # 3. Prints a message listing the unique 'src' values that are not found in 'dt.vertices'.
  dt.edge = as.data.frame(dt.edge)
  n.notinVertices = dt.edge %>% dplyr::anti_join(dt.vertices, by = c("src" = "gene"))
  if(length(unique(n.notinVertices$src)) > 0) {
    print(paste0(
      "!!! nodes:" ,
      paste(unique(n.notinVertices$src), collapse = ","),
      "  do not exist in vertices"
    ))

  }

  # This script processes the 'dt.edge' data frame by adding a new column 'trnsparency'.
  # The 'trnsparency' column is calculated as the ratio of 'conf' to the maximum value of 'conf'.
  dt.edge = dt.edge %>% dplyr::mutate(trnsparency = conf / max(conf))

  # This code checks if the variable 'url.logfile' is not equal to "noLog".
  # If the condition is true, it writes the string 'Step5' to the file specified by 'url.logfile'.
  # The 'append = TRUE' argument ensures that 'Step5' is added to the end of the file without overwriting its existing contents.
  # This section enhances the vertices influx and outflux by removing duplicated columns
  # from the data tables `dt.edge` and `dt.vertices`. The `unique` function is used to
  # ensure that only unique column names are retained. The `setDF` function is then
  # applied to `dt.edge` to convert it to a data frame.
  # Enhance vertices influx outflux----
  dt.edge = dt.edge[, unique(colnames(dt.edge))]
  dt.vertices = dt.vertices[, unique(colnames(dt.vertices))] #removing duplicated columns
  data.table::setDF(dt.edge)
  # This script processes vertex and edge data to calculate the outflux for each vertex.
  #
  # Steps:
  # 1. Create a backup of the original vertex data.
  # 2. Perform an inner join between the vertex data and a subset of the edge data (selecting "src" and "Weight" columns).
  # 3. Group the joined data by vertex ID (gene).
  # 4. Summarize the data to calculate the total outflux (sum of weights) for each vertex.
  #
  # Variables:
  # - dt.vertices.back: Backup of the original vertex data.
  # - dt.edge.tmp: Temporary variable for edge data (currently not used).
  # - dt.vertices.outflux: Data frame containing the outflux for each vertex.
  dt.vertices.back = dt.vertices
  dt.edge.tmp =
    dt.vertices.outflux = dt.vertices %>% dplyr::inner_join(dt.edge[, c("src", "Weight")], by =
                                                       c("gene" = "src")) %>%
    dplyr::group_by(gene) %>% dplyr::summarise(outflux = sum(Weight))

  # This script performs the following operations on the data tables:
  # 1. Calculates the influx for each vertex by joining the 'dt.vertices' table with the 'dt.edge' table on the 'trgt' column,
  #    grouping by 'gene', and summarizing the total 'Weight' for each 'gene'.
  # 2. Merges the 'dt.vertices' table with 'dt.vertices.outflux' on the 'gene' column.
  # 3. Merges the 'dt.vertices' table with the previously calculated 'dt.vertices.influx' on the 'gene' column.
  dt.vertices.influx = dt.vertices %>% dplyr::inner_join(dt.edge[, c("trgt", "Weight")], by =
                                                    c("gene" = "trgt")) %>%
    dplyr::group_by(gene) %>% dplyr::summarise(influx = sum(Weight))

  dt.vertices = dt.vertices %>% dplyr::left_join(dt.vertices.outflux, by = c("gene" =
                                                                        "gene"))
  dt.vertices = dt.vertices %>% dplyr::left_join(dt.vertices.influx, by = c("gene" =
                                                                       "gene"))


  # This script processes a data table of vertices to identify and count missing values (NA).
  # It uses the sapply function to apply a custom function to each element of dt.vertices.
  # The custom function calculates the sum of NA values for each element.
  # The result is stored in variable 'a', which is then filtered to show only elements with missing values.
  # Note: Missing outflux indicates that some nodes do not have a target.
  a = sapply(dt.vertices, function(x)
    sum(is.na(x)))
  a[a > 0] # missing outflux means no target for some nodes and

  dt.vertices = as.data.frame(dt.vertices)
  dt.vertices = replace_na(
    dt.vertices,
    replace = list(
      Pinsum = 0,
      PinAbsSum = 0,
      PinSD = 0,
      PinMean = 0,
      Poutsum = 0,
      PoutAbsSum = 0,
      PoutSD = 0,
      PoutMean = 0,
      outflux = 0,
      influx = 0
    )
  )


  rm(
    dt.vertices.influx,
    dt.vertices.outflux,
    powerParamsIn,
    powerParamsOut,
    dt.edge.tmp
  )
  data.table::setDT(dt.edge)
  dt.edge.noselfedge = dt.edge [src != trgt, .(src, contrib, conf, trgt)]
  data.table::setDF(dt.edge)
  data.table::setDF(dt.vertices)

  selfContribNLMoments = dt.edge %>% dplyr::group_by(src) %>% dplyr::do(extractLmoments(.$contrib))

  # Calculate moments for self contributions and confidence

  dt.edge.back1 <- dt.edge
  dt.vertices.back1 <- dt.vertices

  # Calculate moments for others contributions and confidence
  othersConfsFull <-
    calculate_moments(dt.edge, "trgt", "conf", "ocf", tn)
  othersContribsFull <-
    calculate_moments(dt.edge, "trgt", "contrib", "ocb", tn)
  # Calculate moments for self contributions and confidence
  selfConfsFull <- calculate_moments(dt.edge, "src", "conf", "scf", tn)
  selfContribsFull <-
    calculate_moments(dt.edge, "src", "contrib", "scb", tn)


  dt.vertices.othersparams <-
    othersConfsFull %>% dplyr::inner_join(othersContribsFull, by = c("ocftrgt" = "ocbtrgt"))
  dt.vertices.selfparams <-
    selfConfsFull %>% dplyr::inner_join(selfContribsFull, by = c("scfsrc" = "scbsrc"))


  a <- sapply(dt.vertices.othersparams, function(x)
    sum(is.na(x)))
  a[a > 0]

  #   Integrating moments with dt.vertices -----
  # This script performs the following operations on the 'dt.vertices' data frame:
  # 1. Calculates the number of NA values in each column of 'dt.vertices' and stores the result in 'a'.
  # 2. Filters out columns with no NA values from 'a'.
  # 3. Creates a backup of 'dt.vertices' named 'dt.vertices.BeforeMoments'.
  # 4. Performs a left join of 'dt.vertices' with 'dt.vertices.othersparams' on the 'gene' and 'ocftrgt' columns.
  # 5. Performs a left join of 'dt.vertices' with 'dt.vertices.selfparams' on the 'gene' and 'scbsrc' columns.
  # 6. Recalculates the number of NA values in each column of 'dt.vertices' and stores the result in 'a'.
  # 7. Filters out columns with no NA values from 'a'.
  a = sapply(dt.vertices, function(x)
    sum(is.na(x)))
  a[a > 0]
  dt.vertices.BeforeMoments = dt.vertices

  dt.vertices = dt.vertices %>% dplyr::left_join(dt.vertices.othersparams, by =
                                            c("gene" = "ocftrgt"))
  dt.vertices = dt.vertices %>% dplyr::left_join(dt.vertices.selfparams, by = c("gene" =
                                                                           "scfsrc"))

  a = sapply(dt.vertices, function(x)
    sum(is.na(x)))
  a[a > 0]

  # This script processes the 'dt.vertices' data frame by replacing NA values with specified defaults.
  # The 'my_replace_na' function is used to replace NA values in 'dt.vertices' with the values progeneed in 'rplist'.
  # The replacement values are 0, except for 'MADcost', which is set to 1.488.
  # After processing, several temporary variables related to 'dt.vertices' are removed from the environment
  dt.vertices = my_replace_na(
    dt.vertices,
    rplist = list(
      Mean = 0,
      Median = 0,
      Skew = 0,
      Kurt = 3,
      SD = 0,
      Total = 0,
      NSD = 0,
      flux = 0,
      sum.x = 0,
      sum.y = 0,
      Pinsum = 0,
      PinBAsum = 0,
      Poutsum = 0,
      PoutBAsum = 0,
      MADcost = 1.488,
      L1 = 0,
      L2 = 0,
      tau3 = 0,
      tau4 = 0.1226
    )
  )

  rm(
    dt.vertices.othersparams,
    dt.vertices.selfparams,
    othersConfsFull,
    othersContribsFull
  )

  rm(selfConfsFull, selfContribsFull)

  return(list(dt.edge, dt.vertices, dt.geneexp))
}


extractLmoments = function(v)
{
  if(length(v)==1){
    # if just one value reutrns defualst for normal distirbution
    #Normal 	L1=mean, L2=	σ /√π ,tau3=l_skewness=	0 ,tau4=L_kurtosis=	0.1226
    data.frame(L1=v,L2=0,L3=NA,L4=NA,tau3=0,tau4=0.1226)
  } else
  {
    l =Lmoments::Lmoments(v,returnobject = TRUE)
    if(is.null(l$ratios)) l$ratios=c(NA,NA,0,0.1226)  #happens when length(v)==2
    p=data.frame(L1=l$lambdas[1],L2=l$lambdas[2],L3=l$lambdas[3],L4=l$lambdas[4],
                 tau3=l$ratios[3],tau4=l$ratios[4])
    p
  }
}

removeDups1 <- function(dt, excluded, samplesize = NA) {
  if (!is.na(samplesize))
    dt = dplyr::sample_n(dt, samplesize)
  dupslist = lapply(dt, digest::digest)
  hashs = data.frame(
    cls = names(dupslist),
    hash = unlist(dupslist),
    dups = duplicated(dupslist)
  )
  dups = duplicated(hashs$hash)
  table(dups)

  unique(hashs$cls[dups])
}

calculate_f2 <- function (results3) {
  #same results unpacking process as in calculate_f0
  dt.edge <-
    results3[[1]] # modified dt.edge raw edges with weight values of selected type added
  dt.vertices <-
    results3[[2]] # unmodified dt_geneexp, expression data
  dt.geneexp <- results3[[3]]
  data.table::setDF(dt.vertices)
  L2 = dt.vertices %>% dplyr::ungroup() %>% dplyr::summarise(
    MeanOcfSkew = myMean(ocfSkew, 0, na.rm = TRUE),
    MedianOcfSkew = myMedian(ocfSkew, 0, na.rm = TRUE),
    SDOcfSkew = mySD(ocfSkew, 0, na.rm = TRUE),
    # SD of less than 2 inputs =
    SkewOcfSkew = mySkewness(ocfSkew, 0, na.rm = TRUE),
    KurtOcfSkew = myKurtosis(ocfSkew, 3, na.rm = TRUE),
    MADconstOcfSkew = sn::qsc(.75, MeanOcfSkew, SDOcfSkew, SkewOcfSkew),
    #myQSnorm
    MADOcfSkew = myMedianAbsoluteDeviation(
      ocfSkew,
      MedianOcfSkew,
      MADconstOcfSkew,
      na.rm = TRUE,
      default = 0
    ),

    MeanScbNSkew = myMean(scbNSkew, 0, na.rm = TRUE),
    MedianScbNSkew = myMedian(scbNSkew, 0, na.rm = TRUE),
    SDScbNSkew = mySD(scbNSkew, 0, na.rm = TRUE),
    # SD of less than 2 inputs =
    SkewScbNSkew = mySkewness(scbNSkew, 0, na.rm = TRUE),
    KurtScbNSkew = myKurtosis(scbNSkew, 3, na.rm = TRUE),
    MADconstScbNSkew = sn::qsc(.75, MeanScbNSkew, SDScbNSkew, SkewScbNSkew),
    #myQSnorm
    MADScbNSkew = myMedianAbsoluteDeviation(
      scbNSkew,
      MedianScbNSkew,
      MADconstScbNSkew,
      na.rm = TRUE,
      default = 0
    ),

    MeanOcfNMedian = myMean(ocfNMedian, 0, na.rm = TRUE),
    MedianOcfNMedian = myMedian(ocfNMedian, 0, na.rm = TRUE),
    SDOcfNMedian = mySD(ocfNMedian, 0, na.rm = TRUE),
    # SD of less than 2 inputs =
    SkewOcfNMedian = mySkewness(ocfNMedian, 0, na.rm = TRUE),
    KurtOcfNMedian = myKurtosis(ocfNMedian, 3, na.rm = TRUE),
    MADconstOcfNMedian = sn::qsc(.75, MeanOcfNMedian, SDOcfNMedian, SkewOcfNMedian),
    #myQSnorm
    MADOcfNMedian = myMedianAbsoluteDeviation(
      ocfNMedian,
      MedianOcfNMedian,
      MADconstOcfNMedian,
      na.rm = TRUE,
      default = 0
    ),

    MeanPoutsum = myMean(Poutsum, 0, na.rm = TRUE),
    MedianPoutsum = myMedian(Poutsum, 0, na.rm = TRUE),
    SDPoutsum = mySD(Poutsum, 0, na.rm = TRUE),
    # SD of less than 2 inputs =
    SkewPoutsum = mySkewness(Poutsum, 0, na.rm = TRUE),
    KurtPoutsum = myKurtosis(Poutsum, 3, na.rm = TRUE),
    MADconstPoutsum = sn::qsc(.75, MeanPoutsum, SDPoutsum, SkewPoutsum),
    #myQSnorm
    MADPoutsum = myMedianAbsoluteDeviation(
      Poutsum,
      MedianPoutsum,
      MADconstPoutsum,
      na.rm = TRUE,
      default = 0
    ),


    MeanScbL3 = myMean(scbL3, 0, na.rm = TRUE),
    MedianScbL3 = myMedian(scbL3, 0, na.rm = TRUE),
    SDScbL3 = mySD(scbL3, 0, na.rm = TRUE),
    # SD of less than 2 inputs =
    SkewScbL3 = mySkewness(scbL3, 0, na.rm = TRUE),
    KurtScbL3 = myKurtosis(scbL3, 3, na.rm = TRUE),
    MADconstScbL3 = sn::qsc(.75, MeanScbL3, SDScbL3, SkewScbL3),
    #myQSnorm
    MADScbL3 = myMedianAbsoluteDeviation(
      scbL3,
      MedianScbL3,
      MADconstScbL3,
      na.rm = TRUE,
      default = 0
    )
  )

  a = sapply(L2, function(x)
    sum(x == 0))
  a[a > 0]
  #TODO: add more L2 calculations to dt.vertices
  dt.vertices = dt.vertices %>% dplyr::mutate(
    ZOcfSkew = (ocfSkew - L2$MedianOcfSkew) / L2$MADOcfSkew,
    ZOcfNMedian = (Poutsum - L2$MedianOcfNMedian) / L2$MADOcfNMedian,
    ZScbNSkew = (Poutsum - L2$MedianScbNSkew) / L2$MADScbNSkew,
    ZScbScbL3 = (scbL3 - L2$MedianScbL3) / L2$MADScbL3,
    ZPoutsum = (Poutsum - L2$MedianPoutsum) / L2$MADPoutsum
  )

  a = sapply(dt.vertices, function(x)
    sum(is.na(x)))
  a[a > 0]

  dt.vertices = my_replace_na(
    dt.vertices,
    rplist = list(
      Mean = 0,
      Median = 0,
      Skew = 0,
      Kurt = 3,
      MAD = 0,
      SD = 0,
      Total = 0,
      NSD = 0,
      flux = 0,
      sum.x = 0,
      sum.y = 0,
      MADcost = 1.488,
      L1 = 0,
      L2 = 0,
      L3 = 0,
      L4 = 0,
      tau3 = 0,
      tau4 = 0
    )
  )

  #   Add source and target charectristics to edges -----
  dt.edge.BeforeAugmenting = dt.edge

  collist = c(
    'MAD',
    'Median',
    'MADconst',
    'Mean',
    'SD',
    'Skew',
    'Kurt',
    'Total',
    'NMAD',
    'NMedian',
    'NMADconst',
    'NMean',
    'NSD',
    'NSkew',
    'NKurt',
    'NTotal',
    'tnMedian',
    'tnMADconst',
    'tnMAD',
    'L1',
    'L2',
    'L3',
    'L4',
    'tau3',
    'tau4'
  )
  collist = c(
    paste0("scb", collist),
    paste0("scf", collist),
    paste0("ocb", collist),
    paste0("ocf", collist)
  )
  colList = c(
    'gene',
    'DiagnosesAb',
    'evenTyp',
    'Pinsum',
    'PinBAsum',
    'Poutsum',
    'PoutBAsum',
    'ZOcfSkew',
    'ZPoutsum',
    'ZOcfNMedian',
    'ZScbNSkew',
    'ZScbScbL3',
    collist
  )

  finalCoolList = colList[colList %in% colnames(dt.vertices)]
  dt.vertices1 = dt.vertices[, finalCoolList]


  a = sapply(dt.vertices1, function(x)
    sum(is.na(x)))
  a[a > 0]

  data.table::setDF(dt.edge)
  data.table::setDF(dt.vertices1)
  #!!!NOTE that parameters with .x point to source parameters and those with .y point to target parameters

  dt.edge.back2 = dt.edge
  dt.edge = as.data.frame(dt.edge)
  dt.edge = dt.edge %>% dplyr::left_join(dt.vertices1, by = c("src" = "gene"))
  dt.edge = dt.edge %>%   dplyr::left_join(dt.vertices1, by = c("trgt" = "gene"))

  dupcols = removeDups1(
    dt.edge,
    excluded = c('OcrInp.x', 'Indegree.x'),
    samplesize = nrow(dt.edge) * 0.02
  )
  dt.edge = dt.edge %>% dplyr::select(!any_of(dupcols), any_of(c('Weight', "scftnMedian.x", "scftnMADconst.x", "scftnMAD.x", "ocfNMAD.x", "scftnMedian.y", "scftnMADconst.y", "scftnMAD.y", "ocfNMAD.y")))

  dt.edge = dt.edge %>% dplyr::select(c(-"scbtnMedian.x", -"scbtnMADconst.x", -"scbtnMAD.x", -"scbtnMedian.y", -"scbtnMADconst.y", -"scbtnMAD.y"))

  a = sapply(dt.edge, function(x)
    sum(is.na(x)))
  a[a > 0]

  dt.edge = my_replace_na(dt.edge, rplist =
                            list(
                              L3.x = 0,
                              L3.y = 0,
                              L4.x = 0,
                              L4.y = 0
                            ))

  rm(dt.vertices1)
  #   Adding confz and contribz ----
  if (T) {
    dt.edge = dt.edge %>%
      dplyr::mutate(
        sconfZ.x = (conf - scfMean.x) / ifelse(scfSD.x == 0, PSUDEO_ZERO_2, scfSD.x),
        #Zscore based on Mean and MAD
        scontribZ.x = (contrib - scbMean.x) / ifelse(scbSD.x == 0, PSUDEO_ZERO_2 , scbSD.x) ,

        sconfZ.y = (conf - scfMean.y) / ifelse(scfSD.y == 0, PSUDEO_ZERO_2, scfSD.y),
        #Zscore based on Mean and MAD
        scontribZ.y = (contrib - scbMean.y) / ifelse(scbSD.y == 0, PSUDEO_ZERO_2 , scbSD.y)
      )
  }
  data.table::setDT(dt.edge)

  #   Finalize calculations on edges -----

  dt.edge = dt.edge %>% dplyr::mutate(DiagnosesAb.x = '', DiagnosesAb.y = '')
  rm(dt.edge.BeforeAugmenting)

  #   Correct field names x denotes source and y denotes target ----
  names(dt.edge) <- gsub("_x", ".x", names(dt.edge), fixed = TRUE)
  names(dt.edge) <- gsub("_y", ".y", names(dt.edge), fixed = TRUE)
  return (
    list(
      gene_expression_matrix = dt.geneexp,
      ground_truth = NULL,
      raw_edges = dt.edge,
      edge_features = dt.edge,
      model = NULL,
      model_assessment = NULL,
      predicted_edges = NULL
    )
  )
}
