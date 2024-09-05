### R code from vignette source 'giant_package_vignette.Snw'

###################################################
### code chunk number 1: giant_package_vignette.Snw:45-46 (eval = FALSE)
###################################################
## install.packages("GiANT")


###################################################
### code chunk number 2: giant_package_vignette.Snw:51-55 (eval = FALSE)
###################################################
## # CRAN packages
## install.packages(c("st", "fdrtool"))
## # Bioconductor packages
## BiocManager::install(c("GlobalAncova", "limma", "DESeq2"))


###################################################
### code chunk number 3: giant_package_vignette.Snw:60-61
###################################################
library(GiANT)


###################################################
### code chunk number 4: giant_package_vignette.Snw:79-81
###################################################
set.seed(42)
load("./vignetteData.RData")


###################################################
### code chunk number 5: giant_package_vignette.Snw:84-95 (eval = FALSE)
###################################################
## 
## resGsea <- geneSetAnalysis(
## 	labs = labels,
## 	method = "pearson",
## 	0,
## 	dat = countdata,
## 	geneSets = pathways,
## 	analysis = analysis.gsea(),
## 	adjustmentMethod = "fdr",
## 	signLevel=0.1)
## 


###################################################
### code chunk number 6: giant_package_vignette.Snw:100-101
###################################################
summary(resGsea)


###################################################
### code chunk number 7: giant_package_vignette.Snw:107-108 (eval = FALSE)
###################################################
## tab <- createSummaryTable(resGsea)


###################################################
### code chunk number 8: giant_package_vignette.Snw:110-111
###################################################
tab


###################################################
### code chunk number 9: giant_package_vignette.Snw:114-115 (eval = FALSE)
###################################################
## signtab <- createSummaryTable(resGsea, significantOnly=TRUE, orderBy="geneSetName")


###################################################
### code chunk number 10: giant_package_vignette.Snw:117-118
###################################################
signtab


###################################################
### code chunk number 11: giant_package_vignette.Snw:131-132 (eval = FALSE)
###################################################
## hist(resGsea, subset = 3, aggregate = TRUE)


###################################################
### code chunk number 12: giant_package_vignette.Snw:134-137 (eval = FALSE)
###################################################
## pdf("gsea.pdf")
## hist(resGsea, subset = 3, aggregate = TRUE)
## dev.off()


###################################################
### code chunk number 13: giant_package_vignette.Snw:144-160 (eval = FALSE)
###################################################
## library(parallel)
## mc <- 2 #number of cpus to use
## cl <- makeCluster(mc) #initialize a cluster
## 
## resGsea <- geneSetAnalysis(
## 	labs = labels,
## 	method = "pearson",
## 	numSamples = 1000,
## 	dat = vantVeer,
## 	geneSets = pathways,
## 	analysis = analysis.gsea(),
## 	adjustmentMethod = "fdr",
## 	signLevel=0.1,
## 	cluster = cl)
## 
## stopCluster(cl)


###################################################
### code chunk number 14: giant_package_vignette.Snw:167-168
###################################################
set.seed(132)


###################################################
### code chunk number 15: giant_package_vignette.Snw:171-183 (eval = FALSE)
###################################################
## stat <- abs(apply(vantVeer,1,cor,y = labels))
## coreSet <- rownames(vantVeer)[tail(order(stat), 25)]
## 
## resOverrep <- geneSetAnalysis(
## 	dat = vantVeer,
## 	geneSets = pathways[1:4],
## 	analysis = analysis.customOverrepresentation(),
## 	coreSet = coreSet,
## 	adjustmentMethod = "fdr",
## 	signLevel=0.1)
## 
## summary(resOverrep)


###################################################
### code chunk number 16: giant_package_vignette.Snw:196-199 (eval = FALSE)
###################################################
## pdf("overrepresentation.pdf")
## plotOverrepresentation(resOverrep, subset = 1:4, aggregate = TRUE)
## dev.off()


###################################################
### code chunk number 17: giant_package_vignette.Snw:202-203 (eval = FALSE)
###################################################
## plotOverrepresentation(resOverrep, aggregate = TRUE)


###################################################
### code chunk number 18: giant_package_vignette.Snw:217-231 (eval = FALSE)
###################################################
## resUncertainty <- evaluateGeneSetUncertainty(
## 	#parameters in ...
## 	labs = labels,
## 	numSamples = 1000,
## 	#parameters for evaluateGeneSetUncertainty
## 	dat = vantVeer,
## 	geneSet = pathways[[3]],
## 	analysis = analysis.averageCorrelation(),
## 	numSamplesUncertainty = 100,
## 	k = seq(0.1,0.9,by = 0.1))
## 
## plot(resUncertainty,
## 	main = names(pathways[3]),
## 	addMinimalStability = TRUE)


###################################################
### code chunk number 19: giant_package_vignette.Snw:233-243 (eval = FALSE)
###################################################
## resUncertainty <- evaluateGeneSetUncertainty(
## 	#parameters in ...
## 	labs = labels,
## 	numSamples = 1000,
## 	#parameters for evaluateGeneSetUncertainty
## 	dat = vantVeer,
## 	geneSet = pathways[[3]],
## 	analysis = analysis.averageCorrelation(),
## 	numSamplesUncertainty = 100,
## 	k = seq(0.1,0.9,by = 0.1))


###################################################
### code chunk number 20: giant_package_vignette.Snw:245-250 (eval = FALSE)
###################################################
## pdf("uncertainty.pdf")
## plot(resUncertainty,
## 	main = names(pathways[3]),
## 	addMinimalStability = TRUE)
## dev.off()


###################################################
### code chunk number 21: giant_package_vignette.Snw:269-274
###################################################
myGLS <- function(dat, labs, method = "pearson"){
	return(apply(dat, 1, function(x){
			cor(x = x, y = labs, method = method)
		}))
}


###################################################
### code chunk number 22: giant_package_vignette.Snw:281-284
###################################################
myGSS <- function(x, geneSetIndices){
    return(mean(x[geneSetIndices]))
}


###################################################
### code chunk number 23: giant_package_vignette.Snw:292-306 (eval = FALSE)
###################################################
## myAnalysis <- function(){
## 	return(gsAnalysis(name = "myAnalysis",
## 		gls = "myGLS", 
## 		glsParameterNames = c("labs", "method"),
## 		transformation = "abs", 
## 		transformationParameterNames = NULL,
## 		gss = "myGSS", 
## 		gssParameterNames = NULL,
## 		globalStat = NULL,
## 		globalStatParameterNames = NULL, 
## 		significance = "significance.sampling",
## 		significanceParameterNames = c("numSamples"),
## 		testAlternative = "greater"))
## }


###################################################
### code chunk number 24: giant_package_vignette.Snw:314-322 (eval = FALSE)
###################################################
## myResult <- geneSetAnalysis(
## 	labs = labels,
## 	method = "pearson",
## 	numSamples = 100,
## 	dat = vantVeer,
## 	geneSets = pathways,
## 	analysis = myAnalysis(),
## 	adjustmentMethod = "fdr")


###################################################
### code chunk number 25: giant_package_vignette.Snw:324-325 (eval = FALSE)
###################################################
## hist(myResult)


