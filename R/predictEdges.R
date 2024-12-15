#########################################################################################@
# © 2024 Graeme Vissers
# predict_edges.R takes edge_features and sets up the following
# - The train set
# - The test set
# runRF.R then trains a random forest model using parameters specified by the user
# with the caret R package for machine learning.
#########################################################################################@

#Method could be CICT, and RF or XGB which will only use the basic non CICT features

#' predictEdges
#'
#' Implements CICT supervised learning and prediction of regulatory edges. Currently heavily depends on global variables
#'
#' @param dt.edge CICT edges produced by prepareEdgeFeatures
#' @param edge_features NEW NAME FOR dt.edge
#' @param ground_truth NEW NAME for tbl.goldStandard
#' @param gene_expression_matrix gene expression matrix not required for this function
#' @param in_data_obj cict object replacing rcrd that has all of the results stored as a list
#' @param rcrd A list object that accumulates intermediary objects, results and performance measures and will be stored as an rds file for later use
#' @param method Default value is: "CICT", which uses CICT features in a supervised learning. Also can take "RF" (random forest) or "XGB" (xgboost) which these supervised methods on relevance measures
#' @param evaluateUnseenSample Defualt is: TRUE, calculates regulatory edges on a subsample of up to 50000 edges.
#' @param evaluateAllEdges Defualt is: FALSE. If TRUE, CICT tries predicting all potential edges in the given network which might take longer time.
#' @param Debug If TRUE, function enters debugging mode in critical stops along the exectuion allowing verification of variables and features
#' @param preset.train Defualt is: NA. If provided a path to proper CSV, uses that for training. Useful for sensitivity analysis as well as comparision with other methods on similar set of edges/features
#' @param preset.test Defualt is: NA. If provided a path to proper CSV, uses that for training. Useful for sensitivity analysis as well as comparision with other methods on similar set of edges/features
#' @return Returns a list consisted of three objects
#' rcrd: is a list object of intermediary objects
#' edges: a dataframe of edge objects and CICT features for edges
#' Vertices: a dataframe of vertices objects and CICT features for vertices
#' @examples
#' c(rcrd,edges,vertices) %<-z% prepareEdgeFeatures(Debug=Debug)
#' @export
#'
predictEdges <- function(edge_features = NULL,
                         gene_expression_matrix = NULL,
                         ground_truth = NULL,
                         in_data_obj = NULL,
                         in_format = 'separate',
                         method = 'rf',
                         evaluateUnseenSample = T,
                         evaluateAllEdges = F,
                         Debug = F,
                         url.preset.train = NA,
                         url.preset.test = NA,
                         minGroundTruth.ratio.learning = 0.3,
                         maxGroundTruth = 500,
                         randomEdgesFoldCausal = 5,
                         RF_ntrees = 50,
                         RF_max_depth = 20,
                         exportTrainAndTest = T,
                         returnDat = T,
                         include.negative = 'random',
                         remove.tfs = T,
                         split.val.tfs = F,
                         learning.params = NA,
                         runOnAllEdges = T,
                         trainingTarget = 'class2',
                         tstPercent = 0.3,
                         url.outputFolder='./cict_output/',
                         ...) {
  library(PRROC)
  # PARSE DATA
  {
    if (in_format == "data_obj") {
      edge_features <- in_data_obj$edge_features
      gene_expression_matrix <- in_data_obj$gene_expression_matrix
      ground_truth <- in_data_obj$ground_truth
    }
    out_data_obj <- in_data_obj
  }

  # LABELS CAUSAL, REVERSE CAUSAL, IRRELEVANTA, AND NEGATIVE EXAMPLES
  {
    # t1.c is 'CAUSAL' edges and is the intersect of ground truth and the edges
    # from expression data
    t1.c = ground_truth %>% dplyr::select(src, trgt) %>%
      dplyr::inner_join(edge_features, by = c("src" = "src", "trgt" = "trgt"))
    if (nrow(t1.c) < nrow(ground_truth) / 2)
      warning("Problem in ground truth. More than half of ground truth was not found in edges")
    if (nrow(t1.c) <= 0)
      print("!!! No causal edge in the groundtruth? check gene names")
    t1.c$predicate = "CAUSES"

    # t1.rc is 'REVERSE CAUSAL' edges, and is the same as t1.c, but with src and trgt
    # swapped
    t1.rc = edge_features %>% dplyr::inner_join(t1.c %>% dplyr::select(src, trgt),
                                         by = c("src" = "trgt", "trgt" = "src"))
    t1.rc$predicate = "REV_CAUSES"

    if (include.negative == 'random') {
      # t1.n is 'NEGATIVE' edges to serve as true negative examples in the learning set
      t1.n = edge_features %>%
        dplyr::anti_join(ground_truth, by = c("src" = "src", "trgt" = "trgt")) %>%
        dplyr::filter(src %in% ground_truth$src)
      t1.n$predicate = "NEGATIVE"
    }

    # Joining t1.c and t1.rc
    t1.causal_reversecausal = rbind(t1.c, t1.rc)

    # Adding non-causal edges
    # t1.rnd are edges that are not found in either t1.c or t1.rc (anti-joining by
    # just t1.c should produce the same result)
    t1.rnd = edge_features %>% #dplyr::mutate(edgetyptruth = NA) %>%
      dplyr::anti_join(t1.c, by = c("src" = "src", "trgt" = "trgt")) %>%
      dplyr::anti_join(t1.rc, by = c("src" = "trgt", "trgt" = "src"))

    if (include.negative == 'random') {
      t1.rnd = t1.rnd %>% dplyr::anti_join(t1.n, by = c("src" = "trgt", "trgt" = "src"))
    }
    t1.rnd$predicate = "IRRELEVANTA"

    # t1 consists of:
    # 1. The overlap of causual edges in gt and edges in rawEdges
    # 2. The reverse of (1)
    # 3. All other edges that are in rawEdges and NOT in the gt or its reverse
    t1 = rbind(t1.c, t1.rc, t1.rnd)

    if (include.negative == 'random') {
      t1 = rbind(t1, t1.n)
    }

    # Adds 'class1', 'class2', and 'class3', and src-trgt as 'shared_name'
    # Class1 will either be 'c', 'rc', or 'u'
    # Class2 will be T if Class1 is 'c'
    # Class3 will be T if Class1 is 'c' or 'rc'
    t1$class1 = hutils::Switch(
      t1$predicate,
      CAUSES = 'c',
      REV_CAUSES = 'rc',
      IRRELEVANTA = 'ir',
      IF_NA = 'u',
      DEFAULT = 'u',
      NEGATIVE = 'n'
    )

    t1 = t1 %>%
      dplyr::mutate(
        SUID = dplyr::row_number(),
        class2 = ifelse(class1 %in% c('c'), TRUE, FALSE),
        class3 = ifelse(class1 %in% c('c', 'rc'), TRUE, FALSE),
        shared_name = paste0(src, "-", trgt)
      )

    # Number of causal edges in the ground truth will be the minimum of 'maxGroundTruth'
    # and the minGroundTruth.ratio.learning multiplied by the number of 'causal' edges in t1.c. If this
    # value is less than 300, NcausalEdges will be 2x the initial value.

    # Note that this may bring the value of NcausalEdges HIGHER than maxGroundTruth!!

    NcausalEdges = min(maxGroundTruth,
                       nrow(t1.c) * minGroundTruth.ratio.learning) %>% as.integer()
    # 300 minimum seems arbitrary...
    if (NcausalEdges < 300)
      NcausalEdges = floor(2 * nrow(t1.c) * minGroundTruth.ratio.learning)
    NcausalEdges = min(maxGroundTruth, NcausalEdges) %>% as.integer()
    NrandomEdges = NcausalEdges * randomEdgesFoldCausal

    # TODO: Add preset functionality after non-preset conditions are set
    # If preset is provided, preset.train and preset.test sets are loaded
    # Otherwise, causal edges and reverse-causal edges are randomly selected with NcausalEdges,
    # and random edges are selected at
    if (!(is.na(url.preset.train) | is.na(url.preset.test))) {
      # Load presets
      preset.train <- read.csv(url.preset.train)
      preset.test <- read.csv(url.preset.test)
      names(preset.test)[1:2] <- c('src', 'trgt')
      names(preset.train)[1:2] <- c('src', 'trgt')
      t2 = rbind(
        preset.train %>% dplyr::select(src, trgt) %>% dplyr::inner_join(t1, by = c(
          "src" = "src", "trgt" = "trgt"
        )),
        preset.test %>% dplyr::select(src, trgt) %>% dplyr::inner_join(t1, by =
                                                                  c(
                                                                    "src" = "src", "trgt" = "trgt"
                                                                  ))
      )
    } else{
      t2 = rbind(
        t1 %>% dplyr::filter(class1 == 'c') %>% dplyr::sample_n(size = NcausalEdges),
        t1 %>% dplyr::filter(class1 == 'rc') %>% dplyr::sample_n(size = NcausalEdges),
        t1 %>% dplyr::filter(class1 == 'ir') %>% dplyr::sample_n(size = NrandomEdges)
      )

      # If negative, add negative class
      if (include.negative == 'random') {
        t2 = rbind(t2,
                   t1 %>% dplyr::filter(class1 == 'n') %>% dplyr::sample_n(size = NcausalEdges))
      }
    }

    # t2.complement is everything that is not included in the learning set
    t2.complement = t1  %>%  dplyr::anti_join(t2, by = c("src" = "src","trgt" = "trgt"))

    # Record everything to a log file and in_data_obj object
    msg = sprintf(
      ' Network density= %s, learning set density= %s, \n total edges= %s,
                learning set=%s, test fraction = %s , \n learning set random= %s, learning set causal edges= %s',
      round(nrow(t1.c) / nrow(edge_features), 4),
      round(sum(t2$predicate == 'CAUSES') / nrow(t2), 4),
      prettyNum(nrow(edge_features), ','),
      prettyNum(nrow(t2), ','),
      tstPercent,
      round(NrandomEdges, 2),
      round(NcausalEdges, 2)
    )
    in_data_obj = c(
      in_data_obj,
      c(
        density.net = round(nrow(t1.c) / nrow(edge_features), 4),
        density.learningset = round(sum(t2$predicate == 'CAUSES') /
                                      nrow(t2), 4),
        edges.total = prettyNum(nrow(edge_features), ','),
        set.learning = prettyNum(nrow(t2), ','),
        set.learning.testPcnt = tstPercent,
        set.learning.rndm = round(NrandomEdges, 2),
        set.learning.causal = round(NcausalEdges, 2)
      )
    )
  }

  # SETS TRAINING PARAMETERS
  {
    # Set the parameters for learning and evaluation
    if (!is.na(learning.params)) {
      mdlColNames = learning.params
    } else {
      mdlColNames = colnames(t2)
      mdlColNames = unique(mdlColNames)
    }
    mdlColNames = mdlColNames[mdlColNames %in% colnames(t2)]
    mdlColNames = mdlColNames[!mdlColNames %in% c(
      'src',
      'trgt',
      'class1',
      'SUID',
      'class2',
      'class3',
      'shared_name',
      'predicate',
      'DiagnosesAb.x',
      'DiagnosesAb.y',
      'edgeTyp'
    )]

    objectiveColNames = c("class1", "class2", "class3")
    objectiveColNames = objectiveColNames[objectiveColNames %in% colnames(t2)]

    evaluationColNames = c("predicate", "src", "trgt", "SUID")
    evaluationColNames = evaluationColNames[evaluationColNames %in% colnames(t2)]

    colnamesToExport = unique(c(
      mdlColNames,
      objectiveColNames,
      c('srctrgtSum', 'srctrgtProduct')
    ))
    mdlChkColNames = c(evaluationColNames, objectiveColNames, mdlColNames)

    mdlChkColNames = mdlChkColNames[mdlChkColNames %in% colnames(t2)]
    tst1.totalset = t2[, intersect(colnames(t2), unique(c(mdlChkColNames)))]

    # This removes NAs from the totalset. In this dataframe, only 'R', 'directionLogRatio', and
    # 'iRatio' had NA values. Unclear if this should be kept. Should discuss with Carol
    a = sapply(tst1.totalset, function(x)
      sum(is.na(x)))
    a[a > 0]
    rp = rep(0, length(a[a > 0]))
    names(rp) <- names(a[a > 0])
    as.list(rp)
    tst1.totalset = replace_na(tst1.totalset, as.list(rp))

    # Checking for nas in the test set
    try({
      tst1.totalset = na.omit(tst1.totalset)
    })
  }

  # PARTITIONS DATA FOR LEARNING SET
  {
    ##########################################
    # Creates learning set for random forest training
    ntrgtClass =  nrow(tst1.totalset[trainingTarget == TRUE,])

    # If a preset train and test is provided, use that
    # Otherwise, learning sets are randomly partitioned from the entire tst1.tst set
    if (!(is.na(url.preset.train) | is.na(url.preset.test))) {
      print('Using preset test and preset train')
      tst <-
        preset.train %>% dplyr::inner_join(preset.test, by = c('src', 'trgt'))
      tst1.tst =  preset.test %>% dplyr::select(src, trgt) %>%
        dplyr::inner_join(tst1.totalset, by = c("src" = "src", "trgt" = "trgt"))
      tst1.train = tst1.totalset %>% dplyr::anti_join(tst1.tst, by = c("src" =
                                                                  "src", "trgt" = "trgt"))

    } else{
      while (TRUE) {
        set.seed(as.integer(runif(1, 1, 10000)))
        spltIdx = as.vector(caret::createDataPartition(
          1:nrow(tst1.totalset),
          p = (1 - tstPercent),
          list = FALSE,
          times = 1
        ))
        tst1.train = tst1.totalset[spltIdx,]
        tst1.tst = tst1.totalset[-spltIdx,]
        if (nrow(tst1.tst[trainingTarget == TRUE,]) >= (tstPercent - 0.01) * ntrgtClass)
          break
      }
    }

    tst1.totalset.tfs <- tst1.totalset[tst1.totalset$class2 == T,]
    tst1.totalset.tfs <- as.character(unique(tst1.totalset.tfs$src))

    # Save train and test sets
    # TODO: think about output folder name and if we should create it or throw an error
    if (exportTrainAndTest) {
      # try({
      #   if (!dir.exists(url.outputFolder))
      #     dir.create(url.outputFolder)
      #   tst1.train %>% dplyr::select(src , trgt, class2, class3) %>%
      #     rename(
      #       Gene1 = src,
      #       Gene2 = trgt,
      #       Type = class2,
      #       Association = class3
      #     ) %>%
      #     fwrite(
      #       file = paste0(url.outputFolder, 'train.csv'),
      #       row.names = F,
      #       sep = '\t'
      #     )
      #   tst1.tst %>% dplyr::select(src , trgt, class2, class3) %>%
      #     rename(
      #       Gene1 = src,
      #       Gene2 = trgt,
      #       Type = class2,
      #       Association = class3
      #     ) %>%
      #     fwrite(
      #       file = paste0(url.outputFolder, 'test.csv'),
      #       row.names = F,
      #       sep = '\t'
      #     )
      # })
    }
  }

  # TRAINS RANDOM FOREST MODEL
  {
    # TRAINING CARET  For Prediction and feature selection ----
    # MAX_MEM_SIZE not defined
    # H2OCnn = h2o.init(nthreads = parallel::detectCores(), enable_assertions = TRUE,
    #                   max_mem_size = '100G',strict_version_check=FALSE, port = port)

    if (split.val.tfs) {
      src.groups <- groupKFold(group = tst1.totalset$src, k = 10)
      fit_control <- caret::trainControl(method = "cv",
                                  index = src.groups)
    } else {
      fit_control <- caret::trainControl(method = "cv",
                                  number = 10)
    }

    # Set up df for caret
    caret.y = as.factor(as.data.frame(tst1.totalset)[, trainingTarget])
    caret.x = as.data.frame(tst1.totalset)[, mdlColNames]
    rownames(caret.x) <- tst1.totalset$shared_name

    # Run model on caret
    caret.model <- caret::train(caret.x, caret.y, method = method, trControl = fit_control)

    rownames(tst1.tst) <- tst1.tst$shared_name
    tstset.preds <- caret::extractProb(list(caret.model), tst1.tst)
    tstset.preds$src <- tst1.totalset$src
    tstset.preds$trgt <- tst1.totalset$trgt

    print(tst1.totalset[,trainingTarget])

    out_data_obj$predicted_edges <- tstset.preds
    out_data_obj$variable_importance <- caret::varImp(caret.model)

    # Assigns caret model to model slot
    out_data_obj$model <- caret.model

    # Runs assessment using ground truth and caret functionality


    out_data_obj$model_assessment <- NULL
    print('Data produced successfuly ==================================')
    return(out_data_obj)
  }

  # PREDICT ON ALL EDGES
  {

  }
}