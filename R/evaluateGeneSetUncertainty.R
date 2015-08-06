evaluateGeneSetUncertainty <- function(
		...,
		dat,
		geneSet,
		analysis,
		numSamplesUncertainty,
		k						= seq(0.01, 0.99, by=0.01),
		signLevel 				= 0.05,
		preprocessGeneSet		= FALSE,
		cluster 				= NULL){

	if(analysis$significance != "significance.sampling"){
		stop("'evaluateGeneSetUncertainty': Use an analysis with significance assessment 'significance.sampling'.")
	}

	##########################################
	# init seeds
	##################
	## adapted from
	## package:		parallel
	## function:	clusterSetRNGStream
	##########################################
	if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)){
		oldseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
	}else{
		oldseed <- NULL
	}
	iseed <- round(runif(n=1) * .Machine$integer.max)
	RNGkind("L'Ecuyer-CMRG")
	set.seed(iseed)

	seedList <- vector(mode = "list", length= length(k))
	seedList[[1]] <- .Random.seed
	if(length(k) > 1){
		for (i in 2:length(k)) {
			seedList[[i]] <- nextRNGStream(seedList[[i-1]])
		}
	}
	##########################################

	if(preprocessGeneSet){
		geneSet <- preprocessGs(dat, list(originalGeneSet = geneSet))
		rownames(dat) <- tolower(rownames(dat))
	}else{
		geneSet <- list(originalGeneSet = geneSet)
	}

	##########################################
	# original gene set
	##########################################
	originalGs <- geneSetAnalysis(
		...,
		dat						= dat,
		geneSets				= geneSet,
		analysis				= analysis,
		signLevel				= signLevel,
		adjustmentMethod		= "none")

	transformation <- originalGs$res.all[[1]]$geneSetValues$transformation
	gls <- originalGs$res.all[[1]]$geneSetValues$gls

	##########################################
	# parallel: yes/no
	##########################################
	if(is.null(cluster)){
		myApply <- mapply
	}else{
		clusterCall(cluster, function() library("GiANT", character.only=TRUE))
		clusterExport(cluster,
			c("dat","geneSet","analysis","args","signLevel"),
			envir = environment())
		#export functions possibly defined by the user
		clusterExport(cluster, c(analysis$gls,
			analysis$transformation,
			analysis$gss,
			analysis$globalStat,
			analysis$significance))

		myApply <- function(...){
			clusterMap(cl = cluster, ...)
		}
	}

	##########################################
	# values for #numSamplesUncertainty different
	# (partially) random genesets
	##########################################
	if(!is.null(cluster)){
		clusterExport(cluster,
			c("gls", "transformation", "originalGs"),
			envir = environment())
	}

	gs.stat <- myApply(function(x,s){
			assign(".Random.seed", s, envir = .GlobalEnv)
			
			geneSet <- unlist(geneSet)
			Ngenes <- round(x*length(geneSet))

			# get names which are not in the gene set
			nms <- rownames(dat)
			selectable <- nms[!nms%in%geneSet] 

			allgenesets <- replicate(numSamplesUncertainty, c(geneSet[sample(x=length(geneSet),size=Ngenes)], selectable[sample(x=length(selectable), size=length(geneSet)-Ngenes)]))

			colnames(allgenesets) <- paste("uncertainGeneSet",1:numSamplesUncertainty,sep ="")
			rownames(allgenesets) <- NULL

			gssValues <- apply(allgenesets, 2, doGSS,
				dat = dat,
				analysis = analysis,
				parameters = list(...),
				transformation = transformation)

			names(gssValues) <- paste("iter",1:numSamplesUncertainty,sep ="")

			confidenceInterval <- quantile(c(gssValues,originalGs$res.all[[1]]$geneSetValues$gss), probs = c(signLevel, 0.5, 1-signLevel))
			names(confidenceInterval) <- paste(c(signLevel, 0.5, 1-signLevel)*100, "%-quantile", sep ="")

			return(list(
				confidenceValues = confidenceInterval,
				gssValues = gssValues,
				uncertainGeneSets = allgenesets,
				k = x))
		}, k, seedList, SIMPLIFY=FALSE)

	names(gs.stat) <- paste(k*100, "%-originalGenes", sep ="")

	conf <- t(sapply(gs.stat, "[[", 1))
	rownames(conf) <- paste(k*100, "%-originalGenes", sep ="")

	quant <- c(signLevel, 0.5, 1-signLevel)

	nullDistr <- quantile(originalGs$res.all[[1]]$significanceValues$gssValues, probs = quant)

	ind <- which(abs(conf[,1]-nullDistr[3]) == min(abs(conf[,1]-nullDistr[3])))
	if(originalGs$analysis$testAlternative == "greater"){
		if((conf[,1]-nullDistr[3])[ind] > 0){
			uncertainty <- k[ind]
		}else{
			uncertainty <- k[ind+1]
		}
		#uncertainty <- k[which((abs(conf[,1]-nullDistr[3]) == min(abs(conf[,1]-nullDistr[3]))) & (conf[,1]-nullDistr[3]) > 0)]
	}else if(originalGs$analysis$testAlternative == "less"){
		if((conf[,1]-nullDistr[3])[ind] < 0){
			uncertainty <- 1-k[ind]
		}else{
			uncertainty <- 1-k[ind+1]
		}
		#uncertainty <- k[which((abs(conf[,3]-nullDistr[1]) == min(abs(conf[,3]-nullDistr[1]))) & (conf[,3]-nullDistr[1]) < 0)]
	}

	sstat <- list(
		uncertainty = uncertainty,
		confidenceValues = conf,
		uncertaintyEvaluations = gs.stat,
		signLevel = signLevel,
		originalGeneSetValues = originalGs)

	class(sstat) <- "uncertaintyResult"
	return(sstat)
}
