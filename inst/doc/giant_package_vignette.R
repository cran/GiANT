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
### code chunk number 5: giant_package_vignette.Snw:84-100 (eval = FALSE)
###################################################
## # load data
## require(GlobalAncova)
## data(vantVeer)
## data(phenodata)
## data(pathways)
## 
## 
## resGsea <- geneSetAnalysis(
## 	labs = phenodata$metastases,
## 	method = "pearson",
## 	0,
## 	dat = vantVeer,
## 	geneSets = pathways,
## 	analysis = analysis.gsea(),
## 	adjustmentMethod = "fdr",
## 	signLevel=0.1)


###################################################
### code chunk number 6: giant_package_vignette.Snw:105-106
###################################################
summary(resGsea)


###################################################
### code chunk number 7: giant_package_vignette.Snw:112-113 (eval = FALSE)
###################################################
## tab <- createSummaryTable(resGsea)


###################################################
### code chunk number 8: giant_package_vignette.Snw:115-116
###################################################
tab


###################################################
### code chunk number 9: giant_package_vignette.Snw:119-120 (eval = FALSE)
###################################################
## signtab <- createSummaryTable(resGsea, significantOnly=TRUE, orderBy="geneSetName")


###################################################
### code chunk number 10: giant_package_vignette.Snw:122-123
###################################################
signtab


###################################################
### code chunk number 11: giant_package_vignette.Snw:136-137 (eval = FALSE)
###################################################
## hist(resGsea, subset = 3, aggregate = TRUE)


###################################################
### code chunk number 12: giant_package_vignette.Snw:139-142 (eval = FALSE)
###################################################
## pdf("gsea.pdf")
## hist(resGsea, subset = 3, aggregate = TRUE)
## dev.off()


###################################################
### code chunk number 13: giant_package_vignette.Snw:149-165 (eval = FALSE)
###################################################
## library(parallel)
## mc <- 2 #number of cpus to use
## cl <- makeCluster(mc) #initialize a cluster
## 
## resGsea <- geneSetAnalysis(
## 	labs = phenodata$metastases,
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
### code chunk number 14: giant_package_vignette.Snw:172-173
###################################################
set.seed(132)


###################################################
### code chunk number 15: giant_package_vignette.Snw:176-188 (eval = FALSE)
###################################################
## stat <- abs(apply(vantVeer,1,cor,y = phenodata$metastases))
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
### code chunk number 16: giant_package_vignette.Snw:201-204 (eval = FALSE)
###################################################
## pdf("overrepresentation.pdf")
## plotOverrepresentation(resOverrep, subset = 1:4, aggregate = TRUE)
## dev.off()


###################################################
### code chunk number 17: giant_package_vignette.Snw:207-208 (eval = FALSE)
###################################################
## plotOverrepresentation(resOverrep, aggregate = TRUE)


###################################################
### code chunk number 18: giant_package_vignette.Snw:222-236 (eval = FALSE)
###################################################
## resUncertainty <- evaluateGeneSetUncertainty(
## 	#parameters in ...
## 	labs = phenodata$metastases,
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
### code chunk number 19: giant_package_vignette.Snw:238-248 (eval = FALSE)
###################################################
## resUncertainty <- evaluateGeneSetUncertainty(
## 	#parameters in ...
## 	labs = phenodata$metastases,
## 	numSamples = 1000,
## 	#parameters for evaluateGeneSetUncertainty
## 	dat = vantVeer,
## 	geneSet = pathways[[3]],
## 	analysis = analysis.averageCorrelation(),
## 	numSamplesUncertainty = 100,
## 	k = seq(0.1,0.9,by = 0.1))


###################################################
### code chunk number 20: giant_package_vignette.Snw:250-255 (eval = FALSE)
###################################################
## pdf("uncertainty.pdf")
## plot(resUncertainty,
## 	main = names(pathways[3]),
## 	addMinimalStability = TRUE)
## dev.off()


###################################################
### code chunk number 21: giant_package_vignette.Snw:274-279
###################################################
myGLS <- function(dat, labs, method = "pearson"){
	return(apply(dat, 1, function(x){
			cor(x = x, y = labs, method = method)
		}))
}


###################################################
### code chunk number 22: giant_package_vignette.Snw:286-289
###################################################
myGSS <- function(x, geneSetIndices){
    return(mean(x[geneSetIndices]))
}


###################################################
### code chunk number 23: giant_package_vignette.Snw:297-311 (eval = FALSE)
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
### code chunk number 24: giant_package_vignette.Snw:319-327 (eval = FALSE)
###################################################
## myResult <- geneSetAnalysis(
## 	labs = phenodata$metastases,
## 	method = "pearson",
## 	numSamples = 100,
## 	dat = vantVeer,
## 	geneSets = pathways,
## 	analysis = myAnalysis(),
## 	adjustmentMethod = "fdr")


###################################################
### code chunk number 25: giant_package_vignette.Snw:329-330 (eval = FALSE)
###################################################
## hist(myResult)


