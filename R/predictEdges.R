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
#' @param learning_ratio percent of ground truth to be used for learning
#' @param in_data_obj cict object replacing rcrd that has all of the results stored as a list
#' @param rcrd A list object that accumulates intermediary objects, results and performance measures and
#' will be stored as an rds file for later use.
#' @param preset.train Defualt is: NA. If provided a path to proper CSV, uses that for training.
#' Useful for sensitivity analysis as well as comparision with other methods on similar set of edges/features
#' @param preset.test Defualt is: NA. If provided a path to proper CSV, uses that for training.
#' Useful for sensitivity analysis as well as comparision with other methods on similar set of edges/features
#' @param predict_on Default is: 'none'. User can provide strings 'none', 'all', or a vector
#' of genes on which the user wants to predict.
#' @param sample_tfs Default is: 'random'. Used only when param split_ground_truth_by is set to 'tfs'.
#' User can provide strings 'quartiles' or 'random'. Option 'quartiles' will  randomly sample from each
#' quartile of TFs ranked by number of targets. 'random' will randomly sample from the entire list in an unbiased fashion.
#' @param randomEdgesFoldCausal Numeric representing the number of irrelevant edges to be included
#' in the learning set relative to the number of causal edges. Default is 5.
#' @param negativeEdgesFoldCausal Numeric representing the number of negative edges to be included
#' in the learning set relative to the number of causal edges. Default is 1.
#' @param remove_learning Removes the learning set from the final predictions for performance
#' evaluation. Prevents bias when learning set TFs will always have inflated predictive scores.
#' Default is TRUE.
#' @return Returns a list consisting of three objects
#' edges: a dataframe of edge objects and CICT features for edges
#' Vertices: a dataframe of vertices objects and CICT features for vertices
#' @examples
#' out <- predictEdges(edge_features = SERGIO_DS4_net0_edge_features,
#'                     ground_truth = SERGIO_DS4_net0_ground_truth,
#'                     in_format = 'separate',
#'                     randomEdgesFoldCausal = 5,
#'                     learning_ratio = 0.7)
#' @export
#'
predictEdges <- function(edge_features = NULL,
                         ground_truth = NULL,
                         in_data_obj = NULL,
                         in_format = 'separate',
                         url.preset.train = NA,
                         url.preset.test = NA,
                         learning_ratio = 0.8,
                         maxGroundTruth = 500,
                         randomEdgesFoldCausal = 5,
                         negativeEdgesFoldCausal = 1,
                         exportTrainAndTest = T,
                         returnDat = T,
                         include.negative = 'random',
                         remove.tfs = T,
                         split.val.tfs = F,
                         learning_params = NA,
                         trainingTarget = 'class2',
                         tstPercent = 0.3,
                         predict_on='none',
                         split_ground_truth_by='nodes',
                         sample_tfs='random',
                         remove_learning=TRUE,
                         ...) {
  # PARSE DATA
  {
    if (in_format == "data_obj") {
      edge_features <- in_data_obj$edge_features
      ground_truth <- in_data_obj$ground_truth
    }
  }

  # SUBSETS GROUND TRUTH FOR LEARNING AND EVALUATION
  {
    ground_truth_tfs <- unique(ground_truth$src)
    n_ground_truth_tfs <- length(ground_truth_tfs)
    if (split_ground_truth_by == "nodes") {
      n_tfs_learn <- floor(learning_ratio * n_ground_truth_tfs)
      learning_tfs <- sample(ground_truth_tfs, size = n_tfs_learn)
      eval_tfs <- setdiff(ground_truth_tfs, learning_tfs)
      learning_edges <- ground_truth[ground_truth$src %in% learning_tfs,]
      eval_edges <- ground_truth[ground_truth$src %in% eval_tfs,]
    } else {
      n_edges_learn <- floor(learning_ratio * n_ground_truth_tfs)
      learning_rows <- sample(seq_len(nrow(ground_truth)), size = n_edges_learn)
      eval_rows <- setdiff(seq_len(nrow(ground_truth)), learning_edges)
      learning_edges <- ground_truth[learning_rows,]
      eval_edges <- ground_truth[eval_rows,]
    }
  }

  # LABELS CAUSAL, REVERSE CAUSAL, IRRELEVANTA, AND NEGATIVE EXAMPLES
  {
    # t1.c is 'CAUSAL' edges and is the intersect of the ground truth and the edges
    # from expression data
    t1.c = edge_features %>%
      dplyr::inner_join(ground_truth, by = c("src" = "src", "trgt" = "trgt"))
    if (nrow(t1.c) < nrow(ground_truth) / 2)
      warning("More than half of ground truth was not found in edges")
    if (nrow(t1.c) <= 0)
      print("No causal edge in the ground truth. Check gene names.")
    t1.c$predicate = "CAUSES"
    
    # t1.n is 'NEGATIVE' edges to serve as true negative examples
    # in the learning set
    if (include.negative == 'random') {
      t1.n = edge_features %>%
        dplyr::anti_join(ground_truth, by = c("src" = "src", "trgt" = "trgt")) %>%
        dplyr::filter(src %in% ground_truth$src)
      t1.n$predicate = "NEGATIVE"
    }

    # Adding non-causal edges
    # t1.rnd are edges that are not found in either t1.c
    t1.rnd = edge_features %>%
      dplyr::anti_join(t1.c, by = c("src" = "src", "trgt" = "trgt"))

    if (include.negative == 'random') {
      t1.rnd = t1.rnd %>% dplyr::anti_join(t1.n, by = c("src" = "src", "trgt" = "trgt"))
    }
    t1.rnd$predicate = "IRRELEVANTA"

    # t1 consists of:
    # 1. The overlap of causual edges in gt and edges in rawEdges
    # 2. The reverse of (1)
    # 3. All other edges that are in rawEdges and NOT in the gt or its reverse
    t1 = rbind(t1.c, t1.rnd)

    if (include.negative == 'random') {
      t1 = rbind(t1, t1.n)
    }

    # Adds 'class1', 'class2', and src-trgt as 'shared_name'
    # Class1 will either be 'c', or 'u'
    # Class2 will be T if Class1 is 'c'
    t1$class1 = hutils::Switch(
      t1$predicate,
      CAUSES = 'c',
      IRRELEVANTA = 'ir',
      IF_NA = 'u',
      DEFAULT = 'u',
      NEGATIVE = 'n'
    )

    t1 = t1 %>%
      dplyr::mutate(
        SUID = dplyr::row_number(),
        class2 = ifelse(class1 %in% c('c'), TRUE, FALSE),
        shared_name = paste0(src, "-", trgt)
      )

    # Number of causal edges in the ground truth will be the minimum of 'maxGroundTruth'
    # and the learning_ratio multiplied by the number of 'causal' edges in t1.c. If this
    # value is less than 300, nCausalEdges will be 2x the initial value.

    # Note that this may bring the value of nCausalEdges HIGHER than maxGroundTruth!!

    nCausalEdges <- nrow(learning_edges)
    nRandomEdges <- nCausalEdges * randomEdgesFoldCausal
    if (include.negative == "random") {
      nNegativeEdges <- nCausalEdges * negativeEdgesFoldCausal
    } else {
      nNegativeEdges <- 0
    }

    # TODO: Add preset functionality after non-preset conditions are set
    # If preset is provided, preset.train and preset.test sets are loaded
    # Otherwise, causal edges and reverse-causal edges are randomly selected with nCausalEdges,
    # and random edges are selected at
    if (!(is.na(url.preset.train) | is.na(url.preset.test))) {
      # Load presets
      preset.train <- read.csv(url.preset.train)
      preset.test <- read.csv(url.preset.test)
      names(preset.test)[1:2] <- c('src', 'trgt')
      names(preset.train)[1:2] <- c('src', 'trgt')
      t2 = rbind(
        preset.train %>% dplyr::select(src, trgt) %>%
          dplyr::inner_join(t1, by = c("src" = "src", "trgt" = "trgt")),
        preset.test %>% dplyr::select(src, trgt) %>%
          dplyr::inner_join(t1, by = c("src" = "src", "trgt" = "trgt"))
      )
    } else{
      t2 = rbind(
        t1 %>% dplyr::inner_join(y = learning_edges, by = c("src" = "src", "trgt" = "trgt")),
        t1 %>% dplyr::filter(class1 == 'ir') %>% dplyr::sample_n(size = nRandomEdges)
      )
      # If negative, add negative class
      if (include.negative == 'random') {
        t2 = rbind(t2,
                   t1 %>% dplyr::filter(class1 == 'n') %>% dplyr::sample_n(size = nCausalEdges))
      }
    }

    # t2.complement is everything that is not included in the learning set
    t2.complement = t1  %>%  dplyr::anti_join(t2, by = c("src" = "src","trgt" = "trgt"))

    # Record everything to a log file. This should be coordinated with Stella's log file
    msg = sprintf(
      ' Network density= %s, learning set density= %s, \n total edges= %s,
                learning set=%s, test fraction = %s , \n learning set random= %s, learning set causal edges= %s',
      round(nrow(t1.c) / nrow(edge_features), 4),
      round(sum(t2$predicate == 'CAUSES') / nrow(t2), 4),
      prettyNum(nrow(edge_features), ','),
      prettyNum(nrow(t2), ','),
      tstPercent,
      round(nRandomEdges, 2),
      round(nRandomEdges, 2)
    )
    model_info <- list(
      "Model training info" = list(
        "Percent causal" = round(nrow(t1.c) / nrow(edge_features), 4),
        "Percent causal learning set" = round(sum(t2$predicate == 'CAUSES') /
                                      nrow(t2), 4),
        "Total edges" = prettyNum(nrow(edge_features), ','),
        "Learning set" = prettyNum(nrow(t2), ','),
        "Learning set random" = round(nRandomEdges, 2),
        "Learning set causal" = round(nCausalEdges, 2),
        "Learning set negative" = round(nNegativeEdges, 2)
      )
    )
  }

  # SETS TRAINING PARAMETERS
  {
    # Set the parameters for learning and evaluation
    if (!is.na(learning_params)) {
      mdlColNames <- learning_params
    } else {
      mdlColNames <- colnames(t2)
      mdlColNames <- unique(mdlColNames)
    }
    mdlColNames <- mdlColNames[mdlColNames %in% colnames(t2)]
    mdlColNames <- mdlColNames[!mdlColNames %in% c(
      'src',
      'trgt',
      'class1',
      'SUID',
      'class2',
      'shared_name',
      'predicate',
      'DiagnosesAb.x',
      'DiagnosesAb.y',
      'edgeTyp'
    )]

    objectiveColNames <- c("class1", "class2")
    objectiveColNames <- objectiveColNames[objectiveColNames %in% colnames(t2)]

    evaluationColNames <- c("predicate", "src", "trgt", "SUID")
    evaluationColNames <- evaluationColNames[evaluationColNames %in% colnames(t2)]

    colnamesToExport = unique(c(
      mdlColNames,
      objectiveColNames,
      c('srctrgtSum', 'srctrgtProduct')
    ))
    mdlChkColNames = c(evaluationColNames, objectiveColNames, mdlColNames)
    mdlChkColNames <- mdlChkColNames[mdlChkColNames %in% colnames(t2)]

    tst1.totalset <- as.data.frame(t2)[, mdlChkColNames]
    tst1.complement <- as.data.frame(t2.complement)[, mdlChkColNames]

    # This removes NAs from the totalset. In this dataframe, only 'R', 'directionLogRatio', and
    # 'iRatio' had NA values. Unclear if this should be kept. Should discuss with Carol
    a = sapply(tst1.totalset, function(x)
      sum(is.na(x)))
    rp <- rep(0, length(a[a > 0]))
    names(rp) <- names(a[a > 0])
    rp <- as.list(rp)
    tst1.totalset = tidyr::replace_na(tst1.totalset, rp)

    a = sapply(tst1.complement, function(x)
      sum(is.na(x)))
    a[a > 0]
    rp <- rep(0, length(a[a > 0]))
    names(rp) <- names(a[a > 0])
    rp <- as.list(rp)
    tst1.complement = tidyr::replace_na(tst1.complement, rp)

    # Checking for nas in the test set
    try({
      tst1.totalset = na.omit(tst1.totalset)
      tst1.complement = na.omit(tst1.complement)
    })
  }

  # PARTITIONS LEARNING SET FOR TRAINING AND TESTING
  {
    ##########################################
    # Creates learning set for random forest training
    ntrgtClass <-  nrow(tst1.totalset[trainingTarget == TRUE,])

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
      #   tst1.train %>% dplyr::select(src , trgt, class2) %>%
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
    if (split.val.tfs) {
      src.groups <- caret::groupKFold(group = tst1.totalset$src, k = 5)
      fit_control <- caret::trainControl(method = "cv",
                                  index = src.groups)
    } else {
      fit_control <- caret::trainControl(method = "cv",
                                  number = 5)
    }

    # Set up df for training with caret
    caret.y = as.factor(as.data.frame(tst1.totalset)[, trainingTarget])
    caret.x = as.data.frame(tst1.totalset)[, c('src', 'trgt', mdlColNames)]

    # Set up unknown df for pedictions
    if (length(predict_on) == 1) {
      if (predict_on == 'all') {
        pred.dat <- as.data.frame(rbind(
          tst1.totalset,tst1.complement))[, c('src', 'trgt', mdlColNames)]
      } else if (predict_on == 'none') {
        pred.dat <- caret.x
      } else {
        pred.dat <- as.data.frame(rbind(
          tst1.totalset,tst1.complement))[, c('src', 'trgt', mdlColNames)]
        pred.dat <- pred.dat[pred.dat$src %in% predict_on,]
      }
    } else {
      pred.dat <- as.data.frame(rbind(
          tst1.totalset,tst1.complement))[, c('src', 'trgt', mdlColNames)]
      pred.dat <- pred.dat[pred.dat$src %in% predict_on,]
    }

    # Train model on caret
    caret.model <- caret::train(caret.x, caret.y, method = "rf",
                                trControl = fit_control)

    # Predict on specified edges
    preds <- predict(caret.model, newdata=pred.dat, type='prob')
    preds <- cbind(pred.dat$src, preds)
    preds <- cbind(pred.dat$trgt, preds)
    colnames(preds) <- c("src", "trgt", "FALSE", "Weight")
    preds <- preds[, c("src", "trgt", "Weight")]
    preds <- preds[order(preds$Weight, decreasing = TRUE),]

    # Assigns caret model to model slot
    model_info[[length(model_info) + 1]] <- caret.model
    out_data_obj <- list()
    out_data_obj$model <- model_info

    out_data_obj$model_assessment <- NULL

    out_data_obj$predicted_edges <- preds
  }
  
  # Calculate AUPRC, pAUPRC, AUROC and pAUROC
  {
    if(remove_learning) {
      preds_assess <- preds %>% dplyr::anti_join(y = learning_edges,
                                                 by = c("src" = "src", "trgt" = "trgt"))
    }
    pos_class <- preds_assess %>% dplyr::inner_join(y = eval_edges,
                                             by = c("src" = "src", "trgt" = "trgt"))
    if (nrow(pos_class) == 0) {
      warning("No rows in the evaluation set were found in your predictions. Check gene names.")
      return(out_data_obj)
    }
    neg_class <- preds_assess %>% dplyr::anti_join(y = pos_class,
                                            by = c("src" = "src", "trgt" = "trgt"))
    
    roc <- PRROC::roc.curve(scores.class0 = pos_class$Weight,
                            scores.class1 = neg_class$Weight, curve = TRUE)
    pr <- PRROC::pr.curve(scores.class0 = pos_class$Weight,
                          scores.class1 = neg_class$Weight, curve = TRUE)
    assessments <- list(pr, roc)
    out_data_obj$model_assessment <- assessments
  }
  
  print('========== Data produced successfuly ==========')
  return(out_data_obj)
}