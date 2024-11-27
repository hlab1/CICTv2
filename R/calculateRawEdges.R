# NEW CRE


#Mutual information steady state Multiple measures parallel

  #Parallel partition for each measure, multiple measures

calculateRawEdges <- function(n.workers=5, in_data_obj=NULL, raw_edges=NULL, gene_expression_matrix=NULL, cict_raw_edge_col = 'Spearman',in_format = "separate") {
	nParallelThreads = 12
    edgeTypes <- cict_raw_edge_col

library(dplyr)

library(parallel)

library("infotheo")

library("minet")

library("dtw")

library(tidyr)

library("WGCNA")

library("foreach")

library("doParallel")



########## all functions from before start here


################################################################################@
# This is a modified edition of the code provided by Kuzmanovski et al. to
#  produce raw edge measurements:

# Kuzmanovski V, Todorovski L, Džeroski S. Extensive evaluation of the
#  generalized relevance network approach to inferring gene regulatory networks.
#  Gigascience. 2018 Nov 1;7(11):giy118. doi: 10.1093/gigascience/giy118.
#  PMID: 30239704; PMCID: PMC6420648.
################################################################################@
# Function to build a mutual information matrix (MIM) from a dataset
my.build.mim <- function (
  dataset,
  estimator = "spearman",
  disc = "none",
  nbins = sqrt(NROW(dataset))
){ # Calculate the mutual information matrix
  mim <- calculate_mim(dataset, estimator)
  mim[mim < 0] <- 0
  mim
	# Using pairwise.complete.obs to handle missing values
	mim <- cor(dataset, method = estimator, use = "pairwise.complete.obs")^2
}
# Function to discretize the dataset based on the specified method
discretize_dataset <- function(dataset, disc, nbins) {
  if (disc %in% c("equalfreq", "equalwidth", "globalequalwidth")) {
	dataset <- infotheo::discretize(dataset, disc, nbins)
  }
  dataset
}
# Function to calculate the mutual information matrix using the specified estimator
calculate_mim <- function(dataset, estimator) {
  if (estimator %in% c("pearson", "spearman", "kendall")) {
	mim <- cor(dataset, method = estimator, use = "pairwise.complete.obs")^2
	diag(mim) <- 0
	mim[mim > 0.999999] <- 0.999999
	mim <- -0.5 * log(1 - mim)
  } else {
	estimator <- match_estimator(estimator)
	mim <- infotheo::mutinformation(dataset, method = estimator)
	diag(mim) <- 0
  }
  mim
}

# Function to match the estimator to the appropriate method
match_estimator <- function(estimator) {
  switch(estimator,
		 "mi.mm" = "mm",
		 "mi.empirical" = "emp",
		 "mi.sg" = "sg",
		 "mi.shrink" = "shrink",
		 stop("unknown estimator"))
}

################ Get and setup data ###############################
# Function to read the expression matrix from a file
getExpresionMatrix <- function(filePath, nrows, ncolumns, byrow = TRUE, what = 'character') {
	print("Reading expression data")
	if (file.exists(filePath)) {
		print("File confirmed!")
	}
	E <- scan(file = filePath, what = what, sep = ",")
	EE <- matrix(E, nrows, ncolumns, byrow = byrow)
	return(EE)
}
# This function reads a directed gold standard network from a file and converts it into a matrix format.
getDirectedGoldStandard <- function(ExpresionMatrix, nrows, filePath,byrow = TRUE){
	print("Reading directed gold standard")

	LNetz <- length(scan(file = filePath, what = 'character'))/3
	Netz <- matrix(scan(file = filePath, what = 'character'),LNetz,3,byrow = byrow)
	Netznr <- array(0,c(LNetz,3))
	for (i in 1:LNetz)
	{
 		Netznr[i,1] <- which(Netz[i,1] == ExpresionMatrix[,1])-1
 		Netznr[i,3] <- which(Netz[i,3] == ExpresionMatrix[,1])-1
 		if (Netz[i,2] == "ac")
 		{
  			Netznr[i,2] <- 1
 		}
 		if (Netz[i,2] == "re")
 		{
  			Netznr[i,2] <- -1
 		}
 		if (Netz[i,2] == "du")
 		{
  			Netznr[i,2] <- 0
 		}
	}
	NetzMat <- array(0,c((nrows-1),(nrows-1)))
	for(k in 1:LNetz)
	{
 		NetzMat[Netznr[k,1],Netznr[k,3]] <- 1
	}
	grorg <- graph.adjacency(NetzMat,mode=c("DIRECTED"))

	return(NetzMat)
}

# This function reads an undirected gold standard network from a file and creates a symmetric adjacency matrix.
getUndirectedGoldStandard <- function(ExpresionMatrix, nrows,filePath,byrow = TRUE){
	print("Reading undirected gold standard")

	LNetz <- length(scan(file = filePath, what = 'character'))/3
	Netz <- matrix(scan(file = filePath, what = 'character'),LNetz,3,byrow = byrow)
	Netznr <- array(0,c(LNetz,3))
	for (i in 1:LNetz)
	{
 		Netznr[i,1] <- which(Netz[i,1] == ExpresionMatrix[,1])-1
 		Netznr[i,3] <- which(Netz[i,3] == ExpresionMatrix[,1])-1
 		if (Netz[i,2] == "ac")
 		{
  			Netznr[i,2] <- 1
 		}
 		if (Netz[i,2] == "re")
 		{
  			Netznr[i,2] <- -1
 		}
 		if (Netz[i,2] == "du")
 		{
  			Netznr[i,2] <- 0
 		}
	}
	NetzMatsym <- array(0,c((nrows-1),(nrows-1)))
	for(k in 1:LNetz)
	{
 		NetzMatsym[Netznr[k,1],Netznr[k,3]] <- 1
 		NetzMatsym[Netznr[k,3],Netznr[k,1]] <- 1
	}
	grorgundir <- graph.adjacency(NetzMatsym,mode=c("UNDIRECTED"))

	return(NetzMatsym)
}

# This function calculates the mean expression levels of technical replicates.
getCalculatedReplicates<-function(ExpressionMatrix, ngenes, nexpr, nrepls){
	print("Calculate mean of technical replicates")

	EE_m <- array(0,c((ngenes-1),(nexpr+1)))
	EE_m[,1] <- ExpressionMatrix[2:ngenes,1]
	b1 <- 2
	j <- 1
	for(k in 1:nexpr)
	{
 		for(i in 2:ngenes)
 		{
  			EE_m[i-1,j+1] <- mean(as.numeric(ExpressionMatrix[i,b1:(b1+nrepls-1)]))
 		}
 		j <- j+1
 		b1 <- b1 + nrepls
	}
	EXRATE <- array(0,c((ngenes-1),nexpr))
	EXRATE[,1:nexpr] <- as.numeric(EE_m[,2:(nexpr+1)])

	return(EXRATE)

}


################ Relevance Networks Inferring ###############################
# This function calculates the mutual information between two variables X and Y using a specified method and discretization technique.
mutualinformation <- function(X,Y,methode,discretizers="equalwidth"){
library("infotheo")

 Xd <- unlist(discretize(X,disc=discretizers))
 Yd <- unlist(discretize(Y,disc=discretizers))
 XYd <- array(0,c(length(X),2))
 XYd[,1] <- Xd
 XYd[,2] <- Yd

 I <- entropy(Xd,method=methode) + entropy(Yd,method=methode) - entropy(XYd,method=methode)
 return(I)
}

# This function normalizes an expression matrix using either max or minimax normalization.
getNormalizedMatrix <- function(ExpressionMatrix,normalization="max"){
	##.. normalization: max, minimax

  if(normalization=="max")
	{
		NormalizedMatrix <- ExpressionMatrix/max(ExpressionMatrix)
	}
	else if(normalization=="minimax")
	{
	  minMat <- min(ExpressionMatrix)
	  maxMat <- max(ExpressionMatrix)
		NormalizedMatrix <- (ExpressionMatrix-minMat)/(maxMat-minMat)
	}

	return(NormalizedMatrix)
}
# This function assigns row and column names to a matrix with a specified prefix.
getRowColNamedMatrix <- function(ExpressionMatrix, prefixs="V"){
	rownames(ExpressionMatrix) <- rownames(ExpressionMatrix, do.NULL = FALSE, prefix = prefixs)
	colnames(ExpressionMatrix) <- colnames(ExpressionMatrix, do.NULL = FALSE, prefix = prefixs)

	return(ExpressionMatrix)
}

# This function computes a similarity matrix using various estimators.
getSimilarityMatrix_MI <- function(ExpressionMatrix, nrows, estimators="pearson", subestimators="mm", discretization = FALSE, discretizator = "equalwidth", diagr=0){
  library("minet")
  require(WGCNA)
	##.. estimators[correlation]: pearson, spearman, kendall
	##.. estimators[mutual information]: mi.empirical, mi.mm, mi.shrink, mi.sg
	##.. estimators[other]: coarse.grained, granger
	##.. subestimators[coarse.grained]: ML, MM, Jeffreys, Laplace, SG, minimax, CS, NSB, shrink
	##.. discretizator: equalwidth, equalfreq
	##.. diagr = replacement value of the main diagonal (default: diagr=0)

	if(estimators == "granger")
	{
		source("gc1.R")
		gc_simp <- array(0,c((nrows-1),(nrows-1)))
		for(i in 1:(nrows-1))
		{
 			if(i < (nrows-1))
 			{
  				for(j in (i+1):(nrows-1))
  				{
					T1 <- ExpressionMatrix[i,]
     					T2 <- ExpressionMatrix[j,]

					l1 <- VARselect(T1,lag.max = 8)$selection[[1]]
     					l2 <- VARselect(T2,lag.max = 8)$selection[[1]]
     					if(is.finite(l1) == 0)  l1 <- NA
     					if(is.finite(l2) == 0)  l2 <- NA
     					LAG <- floor(mean(c(l1,l2),na.rm = TRUE))
     					if(is.na(LAG)) LAG <- 1
      				gc_simp[i,j] <- granger(cbind(T2,T1), L=LAG)
      				gc_simp[j,i] <- granger(cbind(T1,T2), L=LAG)
   				}
 			}
		}

		mim <- gc_simp
		diag(mim) <- diagr

	}
	else if(estimators == "coarse.grained")
	{
		DATA <- t(ExpressionMatrix)
		L <- length(DATA[,1])
		Ixy <- array(0,c((nrows-1),(nrows-1)))
		tau_max <- (L-1)
		for(i in 1:(nrows-1))
		{

 			for(j in 1:(nrows-1))
 			{
  				for(tau in -tau_max: tau_max)
  				{

   					if(tau < 0)
   					{
    						X <- DATA[(-tau+1):L,i]
    						Y <- DATA[1:(L+tau),j]
    						I <- mutualinformation(X,Y,subestimators,discretizers = discretizator) #build.mim(data.frame(X,Y), estimator=subestimators)
   					}
   					if(tau > 0)
   					{
    						X <- DATA[1:(L-tau),i]
    						Y <- DATA[(tau+1):L,j]
    						I <- mutualinformation(X,Y,subestimators,discretizers = discretizator) #build.mim(data.frame(X,Y), estimator=subestimators)
   					}
   					if(tau == 0)
   					{
    						I <- 0
   					}
   					Ixy[i,j] <- Ixy[i,j] + I
  				}
  				Ixy[i,j] <- Ixy[i,j]/(2*tau_max)
 			}
		}
		mim <- Ixy
		diag(mim) <- diagr

	}
  else if(estimators =='pearsonFALSE' ){
    #NEEDS debugging
    require(WGCNA)
    allowWGCNAThreads()
    system.time({ mim <- WGCNA::corFast(ExpressionMatrix)})
  }
	else
	{
	  if(discretization == FALSE)
		{
			mim <- build.mim(t(ExpressionMatrix), estimator=estimators)
		}
		else
		{
		  #mim22 <- build.mim(discretize(t(ExpressionMatrix), discretizator), estimator = estimators)
		  mim <- build.mim(t(ExpressionMatrix), disc = discretizator, estimator = estimators)
		}
		diag(mim) <- diagr
	}
	return(mim)
}

# This function calculates a similarity matrix based on specified distance norms.
getSimilarityMatrix_DISTANCES <- function(ExpressionMatrix, nrows, norms=10, diagr=0){
  ##.. norms: 1-> Manhattan, 2-> Euclidian, m->Lm-Norm (default=10)
	# DistanceMatrix <- array(0,c((nrows-1),(nrows-1)))
	#
	# for(di in 1:(nrows-1))
	# {
	#  for(dj in 1:(nrows-1))
	#  {
	# 	DistanceMatrix[di,dj] <- (sum(abs(ExpressionMatrix[di,] - ExpressionMatrix[dj,])^norms))^(1/norms)
	#  }
	# }

  distMatObject <- dist(ExpressionMatrix,method="minkowski",p=norms)

  DistanceMatrix <- as.matrix(distMatObject)
  diag(DistanceMatrix) <- diagr

	return(DistanceMatrix)
}

# This function calculates a similarity matrix using Dynamic Time Warping (DTW) for a given expression matrix.
getSimilarityMatrix_DTW <- function(ExpressionMatrix, distmethod="Euclidean", steppatern=asymmetric){
library("dtw")
	##.. distmethod: Manhattan, Euclidian (default),...
	##.. steppattern: asymmetric(default), symmetric1, symmetric2,...
  DTWMatrix <- dtwDist(ExpressionMatrix,method="DTW", keep.internals=TRUE,step.pattern=steppatern)

	return(DTWMatrix)
}

# This function calculates a similarity matrix based on symbolic representation of an expression matrix.
getSimilarityMatrix_SYMBOLIC <- function(ExpressionMatrix, nrows, npoints, simmethod="sym", npatterns=4, patterns = NULL, diagr=0, discretization = TRUE, discretizator = "equalwidth", mitype="mm", numCores=1){
source("symbolvector.R")
library("minet")
library("parallel")
library("foreach")
library("doParallel")

	##.. simmethod: sym, sym.mu, avg.sym.mi
	##.. npatterns: 1,2,3... number that maximizes no of combination (=npoints/2)
	##.. discretizator: equalwidth, equalfreq
	##.. diagr = replacement value of the main diagonal (default: diagr=0)
    SVEC <- symbolvector(ExpressionMatrix,npoints,npatterns)


	A <- SVEC$A
	A2 <- SVEC$A2
	l_pattern <- SVEC$l_pattern

	if(simmethod=="sym")
	{
        ##setup parallel backend to use numCores cores
        # cl<-makeCluster(numCores, outfile="")
        # registerDoParallel(cl)
        # print("Cluster registered...")

		P1 <- array(0,c((nrows-1),(nrows-1)))
		P2 <- array(0,c((nrows-1),(nrows-1)))
		P3 <- array(0,c((nrows-1),(nrows-1)))

		for(i in 1:(nrows-1))
		{
 			if((i+1) <= (nrows-1))
 			{

                P4i = foreach(j = (i+1):(nrows-1), .combine='c') %do%
                {

                    p1 <- 0
   					p2 <- 0
   					for(pl in 1:l_pattern)
   					{
    						p1 <- p1 + ((length(which(A[,i]==pl & A[,j]==pl)))/length(A[,i]))
    						p2 <- p2 + ((length(which(A[,i]==pl & A2[,j]==pl)))/length(A2[,i]))
   					}

                    P4tmp <- max(c(p1,p2))
                    as.numeric(P4tmp)
                }


                P4pre <- array(NA,c(1,i))
                P3[i,] <- c(P4pre,P4i)
 			}

		}

        # stopCluster(cl)
        # print("Cluster stopped!")

        for(i in 1:(nrows-1)){
            for(j in i:(nrows-1)) {
                if(i == j) {
                    P3[i,j] <- 0
                } else {
                    P3[j,i] <- P3[i,j]
                }
            }
        }

		SimMilarityMatrix <- P3
	}
	else if(simmethod=="sym.mi")
	{
		if(discretization==TRUE)
		{
			MI_A <- mutinformation(discretize(A,disc=discretizator),method=mitype)
		}
		else
		{
			MI_A <- mutinformation(A,method=mitype)
		}
		SimMilarityMatrix <- MI_A
	}
	else if(simmethod=="avg.sym.mi")
	{
        ##setup parallel backend to use numCores cores
        # cl<-makeCluster(numCores, outfile="")
        # registerDoParallel(cl)
        # print("Cluster registered...")

		### Order pattern
		P1 <- array(0,c((nrows-1),(nrows-1)))
		P2 <- array(0,c((nrows-1),(nrows-1)))
		P3 <- array(0,c((nrows-1),(nrows-1)))
		for(i in 1:(nrows-1))
		{
 			if((i+1) <= (nrows-1))
 			{
  				#for(j in (i+1):(nrows-1))
  				#{
   				#	p1 <- 0
   				#	p2 <- 0
   				#	for(pl in 1:l_pattern)
   				#	{
    			#			p1 <- p1 + ((length(which(A[,i]==pl & A[,j]==pl)))/length(A[,i]))
    			#			p2 <- p2 + ((length(which(A[,i]==pl & A2[,j]==pl)))/length(A2[,i]))
   				#	}
   				#	P1[i,j] <- p1
   				#	P1[j,i] <- p1
   				#	P2[i,j] <- p2
   				#	P2[j,i] <- p2
   				#	P3[i,j] <- max(c(P1[i,j],P2[i,j]))
   				#	P3[j,i] <- max(c(P1[j,i],P2[j,i]))
				#
  				#}

                P4i = foreach(j = (i+1):(nrows-1), .combine='c') %do%
                {

                    p1 <- 0
   					p2 <- 0
   					for(pl in 1:l_pattern)
   					{
    						p1 <- p1 + ((length(which(A[,i]==pl & A[,j]==pl)))/length(A[,i]))
    						p2 <- p2 + ((length(which(A[,i]==pl & A2[,j]==pl)))/length(A2[,i]))
   					}

                    P4tmp <- max(c(p1,p2))
                    as.numeric(P4tmp)
                }


                P4pre <- array(NA,c(1,i))
                P3[i,] <- c(P4pre,P4i)
 			}
		}

        # stopCluster(cl)
        # print("Cluster stopped!")

        for(i in 1:(nrows-1)){
            for(j in i:(nrows-1)) {
                if(i == j) {
                    P3[i,j] <- 0
                } else {
                    P3[j,i] <- P3[i,j]
                }
            }
        }



		### Order pattern + mi
		if(discretization==TRUE)
		{
			MI_A <- mutinformation(discretize(A,disc=discretizator),method=mitype)
		}
		else
		{
			MI_A <- mutinformation(A,method=mitype)
		}

		### Finall Avg. Order pattern+mi
		SimMilarityMatrix <- ((P3+MI_A)/2)
	}

	diag(SimMilarityMatrix) <- diagr

	return(SimMilarityMatrix)
}

# This function calculates a similarity matrix using a qualitative approach.
getSimilarityMatrix_QUAL <- function(ExpressionMatrix, nrows, npoints){
source("d_qual.R")
	##.. nrows (a1) is number of genes + 1
	##.. npoints is number of time points within the time series

	D_Q_MATRIX <- Dq_F_MATRIX(ExpressionMatrix,npoints,(nrows-1))
	D_Q <- Dq(D_Q_MATRIX,(nrows-1))

	return(D_Q)
}

getScorredMatrix <- function(SimilarityMatrix, scorrer="MRNET", aracne_eps=0){
library("minet")
	##.. scorrers: mrnet(default), clr, aracne, awe
	##.. aracne_eps is aracne parameter (see minet package manual)
	SimilarityMatrix <- getNormalizedMatrix(SimilarityMatrix,normalization="minimax")

	if(scorrer=="MRNET")
	{
	  ScorredMatrix <- tryCatch({
	    return(mrnet(SimilarityMatrix))
	  },error=function(err){
	    print("Error thrown in MRNET!")
	    print(err);
	    return(NULL);
	  },warning=function(warn){
	    print("Warning thrown in MRNET!")
	    print(warn);
	    return(NULL);
	  })

    #ScorredMatrix <- mrnet(SimilarityMatrix)
	}
	else if(scorrer=="CLR")
	{
	  ScorredMatrix <- tryCatch({
	    return(clr(SimilarityMatrix))
	  },error=function(err){
	    print(SimilarityMatrix)
	    print("Error thrown in CLR!")
	    print(err);
	    return(NULL);
	  },warning=function(warn){
	    print("Warning thrown in CLR!")
	    print(warn);
	    return(NULL);
	  })
	  #ScorredMatrix <- clr(SimilarityMatrix)
	}
	else if(scorrer=="ARACNE")
	{
	  ScorredMatrix <- tryCatch({
	    return(aracne(SimilarityMatrix,eps=aracne_eps))
	  },error=function(err){
	    print("Error thrown in ARACNE!")
	    print(err);
	    return(NULL);
	  },warning=function(warn){
	    print("Warning thrown in ARACNE!")
	    print(warn);
	    return(NULL);
	  })
		#ScorredMatrix <- aracne(SimilarityMatrix,eps=aracne_eps)
	}
	else if(scorrer=="AWE")
	{
		ScorredMatrix <- SimilarityMatrix
		sm_length <- ncol(SimilarityMatrix)
		for(i in 1:sm_length)
		{
 			ScorredMatrix[,i] <- SimilarityMatrix[,i]/sum(SimilarityMatrix[,i])
		}
	}
	else
	{
		ScorredMatrix <- mrnet(SimilarityMatrix)
	}

	return(ScorredMatrix)
}


########### end here
    # TODO: allow config and throw error if in_format is not valid
    if(in_format == "separate") {
      dt_edge <- raw_edges
      dt_geneexp <- gene_expression_matrix
    }
    else if (in_format == "data_obj") {
      dt_edge <- in_data_obj$raw_edges
      dt_geneexp <- in_data_obj$gene_expression_matrix
    }

    print(" - Processing in parallel - ")

    #setwd(url.CICT_algo)
    #source(paste0(url.CICT_algo, 'requirements/rnR_Framework.R'))

    #s0m3
    library(doParallel);

    #outFolder = dirname (url.data)
    actualDataset <- dt_geneexp
    genecol = stringr::str_subset(colnames(actualDataset),'X|^(G|g)ene$')
    if(length(genecol)>0) actualDataset =actualDataset %>% column_to_rownames(genecol)
    actualDataset = actualDataset %>% select_if(is.numeric) #genes in rows and cells in columns  #  stop('Set correct data source') #  all.tdt


    actualDatasetNNodes <- nrow(actualDataset) + 1;
    actualDatasetNObservations <- ncol(actualDataset);
    actualDatasetName <- "actualDataset";
    actualDatasetSymbolicPatterns=0;
    actualDatasetPatterns=0

    simsGroup=c("MI","CORR","DIST")
    availableGroups <- c("MI","CORR","DIST","SYM");
    availableSimilarities <- vector("list", length(availableGroups));
    names(availableSimilarities) <- availableGroups;

    availableSimilarities$MI <- c("ewMImm","ewMIempirical","ewMIshrink","efMIempirical", "efMIshrink") #,"efMImm");
    availableSimilarities$CORR <- c("Pearson","Kendall","Spearman") #);
    availableSimilarities$DIST <- c("Euclidean","Granger"); # ,"L10Norm""Manhattan",
    #availableSimilarities$SYM <- c("efSym","efSymMI","ewSymMI","efAvgSymMI","ewAvgSymMI","Qual");
    #sims <- unlist(lapply(simsGroup, function(group){availableSimilarities[[group]]}), recursive=TRUE, use.names=FALSE);
    sims = intersect(edgeTypes,unlist(availableSimilarities))

    #try({ processingCluster <-getMPIcluster()}) #Uses parallelly loaded by doFuture
    #if(is.null(processingCluster)) processingCluster <-parallelly::makeClusterMPI(n.workers, autoStop = TRUE)
    # processingCluster <-parallel::makeCluster(n.workers) #, autoStop = TRUE)
    # try({
    #   clusterEvalQ(processingCluster, library(infotheo));
    #   clusterEvalQ(processingCluster, library(igraph));
    #   clusterEvalQ(processingCluster, library(minet));
    #   #clusterEvalQ(processingCluster, my.build.mim(actualDataset));
    #   # clusterEvalQ(processingCluster, discretize_dataset(actualDataset, disc = "none", nbins = sqrt(NROW(dataset))));
    #   # clusterEvalQ(processingCluster, calculate_mim(actualDataset, edgeTypes));
    #   # clusterEvalQ(processingCluster, match_estimator(edgeTypes));
    #   # clusterEvalQ(processingCluster, mutualinformation);
    #   # clusterEvalQ(processingCluster, getNormalizedMatrix(actualDataset));
    #   # clusterEvalQ(processingCluster, getRowColNamedMatrix(actualDataset));
    # })

    #If has more cores available then the number of simulations, let them be used on parallelizing subprocesses
    n.workers.subprocess = min(n.workers,length(sims))

    #similarityMatrices <- sapply(
    similarityMatrices <- sapply(
      sims, simplify = FALSE, USE.NAMES = TRUE,
      FUN= function(sim, actualDataset, actualDatasetNNodes, actualDatasetNObservations,
                    actualDatasetName, actualDatasetSymbolicPatterns, patterns, numCores)
        {

          print(paste("[",Sys.time(),"]","Processing",sim,"over",actualDatasetName,"...",sep=" "));
          ##Perform similarity/distance step
          firstStepMatrix <- tryCatch({

            locMatrix <- switch(sim,
                                efMImm = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.mm", discretization = TRUE, discretizator = "equalfreq"),
                                ewMImm = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.mm", discretization = TRUE, discretizator = "equalwidth"),
                                efMIempirical = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.empirical", discretization = TRUE, discretizator = "equalfreq"),
                                ewMIempirical = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.empirical", discretization = TRUE, discretizator = "equalwidth"),
                                efMIshrink = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.shrink", discretization = TRUE, discretizator = "equalfreq"),
                                ewMIshrink = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.shrink", discretization = TRUE, discretizator = "equalwidth"),
                                Pearson = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="pearson"),
                                Spearman = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="spearman"),
                                Kendall = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="kendall"),
                                Manhattan = getSimilarityMatrix_DISTANCES(actualDataset,actualDatasetNNodes,norms=1),
                                Euclidean = getSimilarityMatrix_DISTANCES(actualDataset,actualDatasetNNodes,norms=2),
                                L10Norm = getSimilarityMatrix_DISTANCES(actualDataset,actualDatasetNNodes,norms=10),
                                DTWasym = getSimilarityMatrix_DTW(actualDataset,steppatern=asymmetric),
                                DTWsym1 = getSimilarityMatrix_DTW(actualDataset,steppatern=symmetric1),
                                DTWsym2 = getSimilarityMatrix_DTW(actualDataset,steppatern=symmetric2),
                                efSym = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations, simmethod="sym", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalfreq", patterns=patterns, numCores),
                                ewSym = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations, simmethod="sym", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalwidth", patterns=patterns, numCores),
                                efSymMI = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations,  simmethod="sym.mi", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalfreq", patterns=patterns, numCores),
                                ewSymMI = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations,  simmethod="sym.mi", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalwidth", patterns=patterns, numCores),
                                efAvgSymMI = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations,  simmethod="avg.sym.mi", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalfreq", patterns=patterns, numCores),
                                ewAvgSymMI = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations,  simmethod="avg.sym.mi", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalwidth", patterns=patterns, numCores),
                                #Qual = getSimilarityMatrix_QUAL(actualDataset,actualDatasetNNodes,actualDatasetNObservations)
                                #       efMIml_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="ML",discretizator = "equalfreq"),
                                #       ewMIml_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="ML",discretizator = "equalwidth"),
                                #       efMImm_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="mm",discretizator = "equalfreq"),
                                #       ewMImm_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="mm",discretizator = "equalwidth"),
                                #       efMIshrink_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="shrink",discretizator = "equalfreq"),
                                #       ewMIshrink_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="shrink",discretizator = "equalwidth"),
                                Granger = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="granger")
            )

            row.names(locMatrix) <- row.names(actualDataset);
            colnames(locMatrix) <- row.names(actualDataset);
            print(paste("[",Sys.time(),"]","DONE processing",sim,"over",actualDatasetName,"...",sep=" "));
            return(locMatrix);

          },error=function(err){
            print("Error thrown in thread!")
            print(err);
            return(NULL);
          },warning=function(warn){
            print("Warning thrown in thread!")
            print(warn);
            return(NULL);
          });


          firstStepMatrix
          #getNormalizedMatrix(firstStepMatrix,normalization="minimax");

        },
      actualDataset, actualDatasetNNodes, actualDatasetNObservations, actualDatasetName,
      actualDatasetSymbolicPatterns, actualDatasetPatterns, numCores = n.workers.subprocess
    ); #max(as.numeric((detectCores()-no_cores)/6)-1,1)

    #similarityMatrices=list(res)
    n.itm.e = data.frame()

    #Merge similarityMatrices into a data.frame and return it
    for(i in 1:length(similarityMatrices)){
      tmp=tmp.1=NULL

      sm = similarityMatrices[i]
      if(is.null(sm)) next
      sm.name = names(sm)


      tmp = as.data.frame(sm[[1]])

      if(!is.data.frame(tmp)) next
      if(nrow(tmp)==0) next

      ncol(tmp);nrow(tmp)
      if(all(  rownames(tmp) %>%as.numeric() %>% is.numeric())){
          colnames(tmp) = rownames(actualDataset)
          rownames(tmp) = rownames(actualDataset)
          tmp =tmp %>% mutate(src =rownames(actualDataset)  ) %>% select(src,everything())
      } else tmp = tmp %>% tibble::rownames_to_column()
      # Ensure row names are character
tmp <- tmp %>% mutate(src = rownames(tmp) %>% as.character()) %>% select(src, everything())

# Use pivot_longer to transform data from wide to long format
tmp.1 <- pivot_longer(tmp, cols = colnames(tmp)[2:ncol(tmp)], names_to = 'trgt')

# Rename columns
names(tmp.1) <- c('src', 'trgt', sm.name)

# Filter out rows with NA values in src or trgt
tmp.1 <- tmp.1 %>% dplyr::filter(!is.na(src) & !is.na(trgt))

# Merge with n.itm.e
if (nrow(n.itm.e) == 0) {
  n.itm.e <- tmp.1
} else {
  n.itm.e <- merge(tmp.1, n.itm.e, all.y = TRUE, by = c("src", "trgt"))
}

# Stop the cluster
#parallel::stopCluster(processingCluster)


    #new cols paste0(colnames(n.itm.e),collapse="','")
    #c('efMImm','ewMImm','efMIempirical','ewMIempirical','efMIshrink','ewMIshrink','Pearson','Spearman','Kendall','Manhattan','Euclidean','L10Norm')

    return (list(gene_expression_matrix = actualDataset,
                           ground_truth = NULL,
                           raw_edges = n.itm.e,
                           edge_features = NULL,
                           model = NULL,
                           model_assessment = NULL,
                           predicted_edges = NULL))
  }
  }
