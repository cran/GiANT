##########################################################################
#gene level statistics
#
# dat:			A matrix of expression values. One gene per row.
# labs:			The labels coresponding to the samples. One sample per
#				column.
##########################################################################
gls.cor <- function(dat, labs, method = "pearson"){
	tmp <- apply(dat, 1, function(x){
			cor(x=x, y=labs, method=method)
		})
	names(tmp) <- rownames(dat)
	return(tmp)
}

gls.regression <- function(dat, labs){
	tmp <- apply(dat, 1, function(x){
			coefficients(lm(x~labs))[2]
		})
	names(tmp) <- rownames(dat)
	return(tmp)
}

gls.foldChange <- function(dat, labs, logMeasurements = TRUE){
	cl <- sort(unique(labs))
	cl_ind <- lapply(cl, function(a){which(labs == a)})

	tmp <- apply(dat, 1, function(x){
			means <- sapply(cl_ind, function(a){sum(x[a])})

			if(logMeasurements){
				foldChange <- means[1]-means[2]
			}else{
				foldChange <- means[1]/means[2]
			}
			
		})
	names(tmp) <- rownames(dat)
	return(tmp)
}

gls.tStatistic <- function(dat, labs, pValue = FALSE, alternative = "two.sided"){
	cl <- unique(labs)

	tmp <- apply(dat, 1, function(x){
			res <- t.test(x = x[labs == cl[1]], y = x[labs == cl[2]], alternative = alternative)

			ifelse(pValue, res$p.value, res$statistic)
		})
	names(tmp) <- rownames(dat)
	return(tmp)
}

gls.moderateTStatistic <- function(dat, labs){
	requireNamespace("st")
	tmp <- st::modt.stat(X=t(dat), L=labs)
	names(tmp) <- rownames(dat)
	return(tmp)
}

gls.nBinomTest <- function(dat, labs,
		returnValue = c("pval", "qval", "foldChange", "log2FoldChange"),
		dispersionMethod = "blind",
		dispersionSharingMode = "fit-only",
		dispersionFitType = "local"){

	returnValue <- match.arg(returnValue)
	cl <- sort(unique(labs))

	requireNamespace("DESeq2")

	cds <- DESeq2::DESeqDataSet(dat, labs)   
	cds <- DESeq2::estimateSizeFactors(cds)
	cds <- DESeq2::estimateDispersions(cds,
		method=dispersionMethod,
		sharingMode=dispersionSharingMode,
		fitType=dispersionFitType)

	tmpRes <- as.matrix(DESeq2::results(cds))
	
	names(tmpRes)[which(names(tmpRes) == "pvalue")] <- "pval"
	tmpRes$foldChange <- 2^tmpRes$log2FoldChange
	
	switch(returnValue, 
	       pval={
	         tmp <- tmpRes$pval
	       },
	       qval={
	         tmp <- tmpRes$padj
	       },
	       foldChange={
	         tmp <- tmpRes$foldChange
	       },
	       log2FoldChange={
	         tmp <- tmpRes$log2FoldChange
	       })
	
	
	return(tmp)
}



