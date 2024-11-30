#########################################################################################@
# © 2024 Graeme Vissers
# predict_edges.R takes edge_features and sets up the following
# - The train set
# - The test set
# runRF.R then trains a random forest model using parameters specified by the user
# with the caret R package for machine learning.
#########################################################################################@

#Method could be CICT, and RF or XGB which will only use the basic non CICT features

#' PredictEdges
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
#' # rcrd: is a list object of intermediary objects
#' # edges: a dataframe of edge objects and CICT features for edges
#' # Vertices: a dataframe of vertices objects and CICT features for vertices
#' @examples
#' # Example usage of the function
#' c(rcrd,edges,vertices) %<-z% prepareEdgeFeatures(Debug=Debug)
#' @export
#' 
PredictEdges <-function(edge_features=NULL,
                        gene_expression_matrix=NULL,
                        ground_truth=NULL,
                        in_data_obj=NULL,
                        in_format='separate',
                        method='rf', 
                        evaluateUnseenSample = T, 
                        evaluateAllEdges = F,
                        Debug = F,
                        url.preset.train=NA, 
                        url.preset.test=NA,
                        minGroundTruth.ratio.learning = 0.3,
                        maxGroundTruth=500,
                        randomEdgesFoldCausal=5,
                        RF_ntrees=50,
                        RF_max_depth=20,
                        exportTrainAndTest=T,
                        returnDat=T,
                        include.negative='random',
                        rep = '',
                        remove.tfs = T,
                        split.val.tfs = F,
                        learning.params=NA,
                        runOnAllEdges = T,
                        trainingTarget='class2',
                        tstPercent = 0.3,
                        url.logfile = 'noLog'){
  
  library(dplyr)
  library(hutils)
  library(caret)
  
  # PARSE DATA
  {
    if (in_format == "data_obj"){
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
    t1.c = ground_truth %>% dplyr::select(src,trgt) %>% 
      inner_join(edge_features,by=c("src"="src","trgt"="trgt"))
    if(nrow(t1.c) < nrow(ground_truth)/2) warning("Problem in ground truth. More than half of ground truth was not found in edges")
    if(nrow(t1.c) <=0) print("!!! No causal edge in the groundtruth? check gene names")
    t1.c$predicate = "CAUSES"
    
    # t1.rc is 'REVERSE CAUSAL' edges, and is the same as t1.c, but with src and trgt
    # swapped
    t1.rc = edge_features %>% inner_join(t1.c %>% dplyr::select(src,trgt),by=c("src"="trgt","trgt"="src"))
    t1.rc$predicate = "REV_CAUSES"
    
    if (include.negative == 'random') {
      # t1.n is 'NEGATIVE' edges to serve as true negative examples in the learning set
      t1.n = edge_features %>%
        anti_join(ground_truth,by=c("src"="src","trgt"="trgt")) %>% 
        filter(src %in% ground_truth$src)
      t1.n$predicate = "NEGATIVE"
    }
    
    # This is to generate new port numbers for running H2O. Will remove when caret is used
    if (!is.na(rep)) {
      rep.int <- as.integer(substr(rep, nchar(rep), nchar(rep)))
      port = 54321 + 3*rep.int
    } else {
      port = 54321
    }
    
    # Joining t1.c and t1.rc
    t1.causal_reversecausal=rbind(t1.c,t1.rc)
    
    # Adding non-causal edges
    # t1.rnd are edges that are not found in either t1.c or t1.rc (anti-joining by
    # just t1.c should produce the same result)
    t1.rnd = edge_features %>% #dplyr::mutate(edgetyptruth = NA) %>%
      anti_join(t1.c,by=c("src"="src","trgt"="trgt")) %>%
      anti_join(t1.rc,by=c("src"="trgt","trgt"="src"))
    
    if (include.negative == 'random') {
      t1.rnd = t1.rnd %>% anti_join(t1.n,by=c("src"="trgt","trgt"="src"))
    }
    t1.rnd$predicate = "IRRELEVANTA"
    
    # t1 consists of:
    # 1. The overlap of causual edges in gt and edges in rawEdges
    # 2. The reverse of (1)
    # 3. All other edges that are in rawEdges and NOT in the gt or its reverse
    t1 = rbind(t1.c,t1.rc,t1.rnd)
    
    if (include.negative == 'random') {
      t1 = rbind(t1, t1.n)
    }
    
    # Adds 'class1', 'class2', and 'class3', and src-trgt as 'shared_name'
    # Class1 will either be 'c', 'rc', or 'u'
    # Class2 will be T if Class1 is 'c'
    # Class3 will be T if Class1 is 'c' or 'rc'
    t1$class1=hutils::Switch(t1$predicate, CAUSES = 'c', REV_CAUSES= 'rc',
                             IRRELEVANTA= 'ir', IF_NA= 'u', DEFAULT = 'u',
                             NEGATIVE = 'n')
    
    t1 = t1 %>% 
      dplyr::mutate(SUID=row_number(),
                    class2=ifelse(class1 %in% c('c'),TRUE,FALSE),
                    class3=ifelse(class1 %in% c('c','rc'),TRUE,FALSE),
                    shared_name=paste0(src,"-",trgt))
    
    # Number of causal edges in the ground truth will be the minimum of 'maxGroundTruth'
    # and the minGroundTruth.ratio.learning multiplied by the number of 'causal' edges in t1.c. If this
    # value is less than 300, NcausalEdges will be 2x the initial value.
    
    # Note that this may bring the value of NcausalEdges HIGHER than maxGroundTruth!!
    
    NcausalEdges = min(maxGroundTruth, nrow(t1.c)*minGroundTruth.ratio.learning) %>% as.integer()
    # 300 minimum seems arbitrary...
    if(NcausalEdges<300) NcausalEdges=floor(2*nrow(t1.c)*minGroundTruth.ratio.learning)
    NcausalEdges = min(maxGroundTruth, NcausalEdges) %>% as.integer()
    NrandomEdges = NcausalEdges * randomEdgesFoldCausal
    
    # TODO: Add preset functionality after non-preset conditions are set
    # If preset is provided, preset.train and preset.test sets are loaded
    # Otherwise, causal edges and reverse-causal edges are randomly selected with NcausalEdges,
    # and random edges are selected at 
    if(!(is.na(url.preset.train) | is.na(url.preset.test))) {
      # Load presets
      preset.train <- read.csv(url.preset.train)
      preset.test <- read.csv(url.preset.test)
      names(preset.test)[1:2] <- c('src', 'trgt')
      names(preset.train)[1:2] <- c('src','trgt')
      t2 = rbind( preset.train %>% dplyr::select(src,trgt) %>% inner_join(t1,by=c("src"="src","trgt"="trgt")),
                  preset.test %>% dplyr::select(src,trgt) %>% inner_join(t1,by=c("src"="src","trgt"="trgt")))
    }else{
      t2=rbind(t1 %>% filter(class1=='c') %>% sample_n(size = NcausalEdges),
               t1 %>% filter(class1=='rc') %>% sample_n(size = NcausalEdges),
               t1 %>% filter(class1=='ir')%>% sample_n(size=NrandomEdges))
      
      # If negative, add negative class
      if (include.negative == 'random') {
        t2 = rbind(t2,
                   t1 %>% filter(class1=='n') %>% sample_n(size = NcausalEdges))
      }
    }
    
    # t2.complement is everything that is not included in the learning set
    t2.complement = t1  %>%  anti_join(t2,by=c("src"="src","trgt"="trgt"))
    
    # Record everything to a log file and in_data_obj object
    msg=sprintf(' Network density= %s, learning set density= %s, \n total edges= %s, 
                learning set=%s, test fraction = %s , \n learning set random= %s, learning set causal edges= %s',  
                round(nrow(t1.c)/nrow(edge_features),4),round(sum(t2$predicate=='CAUSES')/nrow(t2),4),
                prettyNum(nrow(edge_features),','),prettyNum(nrow(t2),','),tstPercent,
                round(NrandomEdges,2),round(NcausalEdges,2))
    
    print(url.logfile)
    cat(msg)
    if(url.logfile!="noLog") write(msg, file = url.logfile, append = TRUE)
    
    in_data_obj=c(in_data_obj, c(density.net=round(nrow(t1.c)/nrow(edge_features),4),
                   density.learningset=round(sum(t2$predicate=='CAUSES')/nrow(t2),4),
                   edges.total=prettyNum(nrow(edge_features),','),
                   set.learning=prettyNum(nrow(t2),','),
                   set.learning.testPcnt=tstPercent,
                   set.learning.rndm=round(NrandomEdges,2),
                   set.learning.causal =round(NcausalEdges,2)))
    
    write(table(t2$predicate) %>% knitr::kable(),  file = url.logfile, append = TRUE)
  }
  
  # SETS TRAINING PARAMETERS
  {
    # Set the parameters for learning and evaluation
    if (!is.na(learning.params)) {
      mdlColNames=learning.params
    } else {
      mdlColNames=colnames(t2)
      mdlColNames=unique(mdlColNames)
    }
    mdlColNames=mdlColNames[mdlColNames %in% colnames(t2)]
    mdlColNames=mdlColNames[!mdlColNames %in% c('src','trgt', 'class1', 'SUID', 'class2', 'class3', 'shared_name', 'predicate',
                                                'DiagnosesAb.x', 'DiagnosesAb.y', 'edgeTyp')]
    
    objectiveColNames=c("class1","class2","class3")
    objectiveColNames=objectiveColNames[objectiveColNames %in% colnames(t2)]
    
    evaluationColNames=c("predicate","src","trgt","SUID")
    evaluationColNames=evaluationColNames[evaluationColNames %in% colnames(t2)]
    
    colnamesToExport= unique(c(mdlColNames,objectiveColNames, c('srctrgtSum','srctrgtProduct')))
    mdlChkColNames = c(evaluationColNames,objectiveColNames,mdlColNames)
    
    mdlChkColNames=mdlChkColNames[mdlChkColNames %in% colnames(t2)]
    tst1.totalset = t2[,intersect(colnames(t2),unique(c(mdlChkColNames)))]
    
    # This removes NAs from the totalset. In this dataframe, only 'R', 'directionLogRatio', and
    # 'iRatio' had NA values. Unclear if this should be kept. Should discuss with Carol
    a=sapply(tst1.totalset, function(x) sum(is.na(x)));a[a>0]
    rp = rep(0,length(a[a>0]))
    names(rp)<-names(a[a>0])
    as.list(rp)
    tst1.totalset =replace_na(tst1.totalset, as.list(rp))
    
    # Checking for nas in the test set
    try({
      tst1.totalset = na.omit(tst1.totalset)
    })
  }
  
  # PARTITIONS DATA FOR LEARNING SET
  {
    ##########################################
    # Creates learning set for random forest training
    ntrgtClass =  nrow(tst1.totalset[trainingTarget== TRUE,])
    
    # If a preset train and test is provided, use that
    # Otherwise, learning sets are randomly partitioned from the entire tst1.tst set
    if(!(is.na(url.preset.train) | is.na(url.preset.test))){
      print('Using preset test and preset train')
      tst <- preset.train %>% inner_join(preset.test, by=c('src', 'trgt'))
      tst1.tst =  preset.test %>% dplyr::select(src,trgt) %>%
        inner_join(tst1.totalset,by=c("src"="src","trgt"="trgt"))
      tst1.train = tst1.totalset %>% anti_join(tst1.tst,by=c("src"="src","trgt"="trgt"))
      
    } else{
      while(TRUE){
        set.seed(as.integer(runif(1,1,10000)))
        spltIdx =as.vector( caret::createDataPartition(1:nrow(tst1.totalset),p=(1-tstPercent),list=FALSE,times=1))
        tst1.train = tst1.totalset[spltIdx,]
        tst1.tst = tst1.totalset[-spltIdx,]
        
        print(paste0("nrow(tst1.tst[trainingTarget]", nrow(tst1.tst[trainingTarget== TRUE,])))
        if(nrow(tst1.tst[trainingTarget== TRUE,]) >= (tstPercent-0.01)* ntrgtClass) break
      }
    }
    
    tst1.totalset.tfs <- tst1.totalset[tst1.totalset$class2 == T,]
    tst1.totalset.tfs <- as.character(unique(tst1.totalset.tfs$src))
    
    # Record train and test sets for learning
    tmp = table(tst1.train$class2); names(tmp)<-paste0('train.',names(tmp))
    in_data_obj=c(in_data_obj, tmp)
    tmp = table(tst1.tst$class2); names(tmp)<-paste0('test.',names(tmp))
    in_data_obj=c(in_data_obj, tmp)
    
    # Save train and test sets
    if(exportTrainAndTest){
      try({
        if(!dir.exists(url.outputFolder)) dir.create(url.outputFolder)
        tst1.train %>% dplyr::select(src , trgt, class2, class3) %>% 
          rename(Gene1 = src, Gene2 = trgt, Type=class2,Association=class3) %>% 
          fwrite(file=paste0(url.outputFolder,'train.csv'),row.names = F, sep='\t')
        tst1.tst %>% dplyr::select(src , trgt, class2, class3) %>% 
          rename(Gene1 = src, Gene2 = trgt, Type=class2,Association=class3) %>% 
          fwrite(file=paste0(url.outputFolder,'test.csv'),row.names = F, sep='\t')
      })
    }
  }
  
  # TRAINS RANDOM FOREST MODEL
  {
    # TRAINING CARET  For Prediction and feature selection ----
    # MAX_MEM_SIZE not defined
    # H2OCnn = h2o.init(nthreads = parallel::detectCores(), enable_assertions = TRUE,
    #                   max_mem_size = '100G',strict_version_check=FALSE, port = port)
    
    if (split.val.tfs) {
      src.groups <- groupKFold(group=tst1.totalset$src, k=10)
      fit_control <- trainControl(
        method = "cv",
        index = src.groups)
    } else {
      fit_control <- trainControl(
        method = "cv",
        number = 10)
    }
    
    # Set up df for caret
    caret.y = as.factor(as.data.frame(tst1.totalset)[,trainingTarget])
    caret.x = as.data.frame(tst1.totalset)[,mdlColNames]
    rownames(caret.x) <- tst1.totalset$shared_name
    
    # Run model on caret
    tst1.rf <- caret::train(caret.x, caret.y, method=method, trControl = fit_control)
    out_data_obj$rf_outputs <- tst1.rf
    
    rownames(tst1.tst) <- tst1.tst$shared_name
    caret.test <- predict(tst1.rf, tst1.tst, type="prob")

    out_data_obj$predicted_edges <- caret.test
    
    # TODO
    out_data_obj$model_assessment <- NULL
    print('Data produced successfuly ==================================')
    return(out_data_obj)
    # 25000 random sample & bestCutoff calculations -----------------
    { 
      if(F) {
        #Predicts in smaller chunks also adds random prediction for comparison
        msg = c('================================================',
                "Reporting results on unseen sample")
        cat(paste0(msg,collapse='\n') )
        write( msg,  file = url.logfile, append = TRUE, sep = '\n')
        
        # d.new is a sample of t2.complement, which is all of the unlabeled data. Most will be labeled
        # as IR with very few CAUSAL relationships
        d.new = t2.complement %>% sample_n(size = min(50000,nrow(t2.complement) * maxunseenTest.ratio))
        print('d.new$predicate')
        print(table(d.new$predicate))
        
        print('Use tst1.tst for assessment instead')
        d.new = tst1.tst
        
        msg = sprintf("Learning set= %s | training= %s | validation= %s | unseen sample = %s | t2 + comp = %s | all Edges= %s",  
                      nrow(tst1.totalset),nrow(tst1.train),nrow(tst1.tst),nrow(d.new), nrow(t2)+nrow(t2.complement),nrow(edge_features))
        
        cat(msg)
        if(url.logfile!="noLog") write(msg, file = url.logfile, append = TRUE)
        in_data_obj$unseensmpl_stats = msg
        
        
        #d.new = edge_features
        splitcount = 2
        splts = split(1:nrow(d.new),          
                      cut(seq_along(1:nrow(d.new)),
                          splitcount,labels = FALSE))
        
        #Creates and saves random predictions for all edges
        randomPredictions = edge_features %>% dplyr::select(src,trgt) %>% mutate(rndPred = ifelse(runif(nrow(edge_features)) >=.5,1,0))
        
        stop()
        d.new.tmp = lapply(splts,
                           function(thesplit){
                             print(last(thesplit))
                             d.new.slice = d.new[thesplit,]
                             predTest.d.h2o = as.h2o(d.new.slice) # 
                             h2o.pred = as.data.frame(h2o.predict(tst1.mdl,predTest.d.h2o,keep_cross_validation_predictions=F))
                             h20.prediction=as.numeric(as.character(h2o.pred[,3]))
                             predictions =h20.prediction #pred #ens.predictions #
                             outcomes =unlist(setDT(d.new.slice)[,trainingTarget]) #outcome # ens.outcome# 
                             set.seed(runif(1,1,1000))
                             # = ifelse(runif(nrow(d.new.slice)) >=.5,1,0) #Assign a random classifier results
                             prd_outcomes = d.new.slice %>% dplyr::select(src,trgt,any_of('Weight')) %>%
                               cbind(predictions,outcomes) %>% as.data.frame()
                             prd_outcomes
                           })
        
        d.new1= rbindlist(d.new.tmp) %>% left_join(randomPredictions,by=c('src'='src','trgt'='trgt'))
        
        
        d.new1.rv =d.new1 %>% dplyr::rename(revWeight = Weight,rvpred=predictions,src1=src,trgt1=trgt) %>% 
          dplyr::select(-outcomes,-rndPred)
        pred_outcome = d.new1 %>% left_join(d.new1.rv, by=c("src"="trgt1", "trgt"="src1") )
        pred_outcome.back = pred_outcome
        
        
        
        
        pred_outcome.top2c = pred_outcome[, head(.SD, 2), by = "trgt"]
        pred_outcome.top2c=pred_outcome.top2c[,`:=`(predictions=NULL,outcomes=NULL,rvpred=NULL , is.causal1 =1)]
        d.new2 = merge(d.new1,pred_outcome.top2c,all.x=TRUE, by=c('src','trgt'))
        table(d.new2$outcomes,d.new2$is.causal1)
        
        
        pred_outcome.top1c = pred_outcome[, head(.SD, 1), by = "trgt"]
        pred_outcome.top1c=pred_outcome.top1c[,`:=`(predictions=NULL,outcomes=NULL,rvpred=NULL , is.causal1 =1)]
        
        d.new2 = merge(d.new1,pred_outcome.top1c,all.x=TRUE, by=c('src','trgt'))
        table(d.new2$outcomes,d.new2$is.causal1)
      }
      
      
      
      # Testing performance -----
      {
        prd.varimp=h2o.varimp(tst1.mdl)
        
        in_data_obj$varimp = prd.varimp[1:20,] %>% as.data.frame()
        
        h2o.performance(tst1.mdl,valid=T)
        
        msg = c('================================================',
                "Reporting model performance on validation set",
                capture.output(h2o.performance(tst1.mdl,valid=T)))
        cat(paste0(msg,collapse='\n') )
        
        setDF(tst1.totalset)
        
        require(verification)
        require(pROC)
        require(ROCR)
        require( OptimalCutpoints)
        require(precrec )
        
        reportAUC<-function(x)
        {
          a= attr(x,'aucs');
          b=attr(x,'paucs') %>% rename(aucs=paucs,standardized = spaucs) %>% 
            mutate(curvetypes = paste0('p',curvetypes))
          
          rbindlist(list(a,b),fill=T,use.names = T) %>% dplyr::select(-modnames, -dsids) %>%
            mutate(aucs = round(aucs,3), standardized = round(standardized,3) ) 
        }
        
        #TODO remove the reverse edges before assessment, just keep the one with higher prediction score
        assespreds = pred_outcome #%>% dplyr::filter(predictions>=rvpred)
        print(head(assespreds))
        print(table(assespreds$outcomes))
        theROC <- roc(assespreds$outcomes, assespreds$predictions, percent = TRUE);theROC
        
        sscurves <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes)
        pr.crv=sscurves$prcs[1][[1]];  pr.auc =  attr(pr.crv,'auc')
        
        theauc = precrec::auc(sscurves); 
        msg = paste0(theauc[[3]],"=", round(theauc[[4]],3))
        
        in_data_obj$unseensmpl_roc_pr = msg
        
        #partial precision-recall
        sscurves.part <- part(sscurves, xlim = c(0, 0.2))
        in_data_obj$unseensmpl_part = as.data.frame(reportAUC(sscurves.part) )
        
        pr.prtcrv=sscurves.part$prcs[1][[1]];  pr.prtauc =  attr(pr.prtcrv,'pauc')
        
        #random classifier
        randomClassifierCurves <- evalmod(scores = assespreds$rndPred, labels = assespreds$outcomes)
        rnd.crv=randomClassifierCurves$prcs[1][[1]]; rnd.auc =  attr(rnd.crv,'auc')
        
        rndmClscurves.part <- part(randomClassifierCurves, xlim = c(0, 0.2))
        reportAUC(rndmClscurves.part) 
        rnd.prtcrv=rndmClscurves.part$prcs[1][[1]]; rnd.prtauc =  attr(rnd.prtcrv,'pauc')
        in_data_obj$unseensmpl_rndm = as.data.frame(reportAUC(rndmClscurves.part) )
                
        mmpoins <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes, mode = "basic")
        
        # the relative cost of of a false negative classification (as compared with a false positive classification)
        # the prevalence, or the proportion of cases in the population (n.cases/(n.controls+n.cases)).
        relativeCostfn_fp = 1/2
        in_data_obj$relativeCostfn_fp=relativeCostfn_fp
        
        prv = table(setDT(tst1.totalset)[,get(trainingTarget)])
        best.weights=c(relativeCostfn_fp, prv[2]/prv[1])  
        bestcutoff =as.double(coords(theROC, "best", best.method = "closest.topleft", best.weights=best.weights,
                                     ret="threshold", transpose = FALSE));bestcutoff
        
        assespreds =assespreds %>% mutate(thresholdpreds= ifelse(assespreds$predictions>bestcutoff,assespreds$predictions,0)) 
        theROC <- roc(assespreds$outcomes, assespreds$thresholdpreds, percent = TRUE);theROC
        
        
        
        sscurves <- evalmod(scores = assespreds$thresholdpreds, labels = assespreds$outcomes)
        pr.crv=sscurves$prcs[1][[1]];  pr.auc =  attr(pr.crv,'auc')
        
        msg=""
        if(url.logfile!="noLog") {
          print(sprintf("Best threshold: %s for fn/fp cost ratio: %s " ,bestcutoff,relativeCostfn_fp ))
          write(sprintf("Best threshold: %s for fn/fp cost ratio: %s " ,bestcutoff,relativeCostfn_fp ),  file = url.logfile, append = TRUE)
          
          msg = sprintf("With FP-FN ratio 0.5 => AUCPR Ratio CICT to Random= %s,  top 20 percent AUCPR Ratio CICT to Random= %s ", 
                        round(pr.auc/rnd.auc,2), round(pr.prtauc/rnd.prtauc,2))
          write(msg,  file = url.logfile, append = TRUE)
          
        }
        
        print(msg)
        in_data_obj$costFnFp_cutoff_roc = msg
      }
    }
  }
  
  # PREDICTS ON ALL UNSEEN EDGES
  if(runOnAllEdges) {
    write('Edge ranking',  file = url.logfile, append = TRUE)
    d.new = edge_features
    
    if (remove.tfs == TRUE) {
      d.new <- d.new[!d.new$src %in% tst1.totalset.tfs,]
    }
    
    splitcount = max(floor(nrow(edge_features) / 30000),2)
    splts = split(1:nrow(d.new),          
                  cut(seq_along(1:nrow(d.new)),
                      splitcount,labels = FALSE))
    
    d.new.tmp = lapply(
      splts,
      function(thesplit){
        print(last(thesplit))
        d.new.slice = setDT(d.new)[thesplit,]
        caret.pred = as.data.frame(predict(tst1.mdl,predTest.d.h2o,keep_cross_validation_predictions=F))
        h2o.pred = as.data.frame(h2o.predict(tst1.mdl,predTest.d.h2o,keep_cross_validation_predictions=F))
        h20.prediction=as.numeric(as.character(h2o.pred[,3]))
        predictions =h20.prediction
        set.seed(runif(1,1,1000))
        rndPred = ifelse(runif(nrow(d.new.slice)) >=.5,1,0) #Assign a random classifier results
        prd_outcomes = d.new.slice %>% dplyr::select(src,trgt,any_of('Weight')) %>%
          cbind(predictions,rndPred) %>% as.data.frame()
        prd_outcomes
      })
    d.new1= rbindlist(d.new.tmp)
    
    write('Edges ranking finished',  file = url.logfile, append = TRUE)
    
    if (F) {
      d.new1.rv =d.new1 %>% dplyr::rename(revWeight = Weight,rvpred=predictions,src1=src,trgt1=trgt) %>% 
        dplyr::select(-rndPred)
      pred_outcome = d.new1 %>% left_join(d.new1.rv, by=c("src"="trgt1", "trgt"="src1"))
    }
    pred_outcome = d.new1
    print('any(tst1.totalset.tfs %in% pred_outcome$src)')
    print(any(tst1.totalset.tfs %in% pred_outcome$src))
    
    print('head(pred_outcome)')
    print(head(pred_outcome))
    
    if(exportRankedEdges) {
      pred_outcome.e=pred_outcome %>% rename(Gene1=src,	Gene2=trgt,	EdgeWeight = predictions) %>%
        dplyr::select(Gene1,Gene2,EdgeWeight,everything()) %>% arrange(desc(EdgeWeight))
      
      if(!file.exists(url.rankedEdges) | forceOutput){
        file.remove(url.rankedEdges)
        try({fwrite(pred_outcome.e,url.rankedEdges,row.names = F, sep='\t')},silent=T)
        
        in_data_obj$rankededges_count = nrow(pred_outcome.e)
        in_data_obj$rankededges_gated_count = nrow(pred_outcome.e)
        
        
        truncated_pred_outcome = pred_outcome.e %>% 
          mutate(EdgeWeight=ifelse(EdgeWeight>bestcutoff,EdgeWeight,0))
        in_data_obj$rankededges_gated_count = nrow(truncated_pred_outcome)
        in_data_obj$rankededges_gated_desc = paste0('EdgeWeight > ',bestcutoff )
        in_data_objrankededges_url =url.rankedEdges 
        
        try({fwrite(truncated_pred_outcome,url.rankedEdgesGated,row.names = F, sep='\t')})
      }
      
    }
  }
  print('Data produced successfuly ==================================')
  return(out_data_obj)
}
