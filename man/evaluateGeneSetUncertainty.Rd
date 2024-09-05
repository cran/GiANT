\name{evaluateGeneSetUncertainty}
\alias{evaluateGeneSetUncertainty}
\title{
Quantify gene set uncertainty
}
\description{
A robustness measure that quantifies the uncertainty of a gene set  by performing a resampling experiment and can be used in the \code{robustness} parameter of \code{gsAnalysis}.
}
\usage{
evaluateGeneSetUncertainty(
	...,
	dat,
	geneSet,
	analysis,
	numSamplesUncertainty,
	blockSize = 1,
	k = seq(0.01, 0.99, by=0.01),
	signLevel = 0.05,
	preprocessGeneSet = FALSE,
	cluster = NULL)
}
\arguments{
\item{...}{
	Additional parameters for the different steps of the analysis pipeline, depending on the concrete configuration supplied in \code{analysis}.
}
\item{dat}{
	A numeric matrix of gene expression values for all analyzed genes. Here, each row corresponds to one gene, and each column corresponds to one sample. The rows must be named with the gene names used in the gene sets.
}
\item{geneSet}{
	A vector containing the names of genes in a gene set. All genes set must correspond to the row names of \code{dat}.
}
\item{analysis}{
	The parameters of the analysis that is applied to the perturbed copies of the gene set. These parameters are described by an object of class \code{gsAnalysis} as returned by the function \code{\link{gsAnalysis}} or the predefined analysis descriptors in \code{\link{predefinedAnalyses}}.
}
\item{numSamplesUncertainty}{
	The number of resampling experiments which should be applied to estimate the robustness of \code{geneSet}.
}
\item{blockSize}{
	Number of genes in one resampled block.
}
\item{k}{
	A \code{vector} of percentages of genes in the randomized gene sets that should be taken from the original gene set. The remaining genes are chosen randomly. For each value a resampling experiment is performed.
}
\item{signLevel}{
	The significance level for the significance assessment of the gene sets (defaults to \code{0.05}). 
}
\item{preprocessGeneSet}{
	Specifies whether the gene sets in \code{geneSets} should be preprocessed or not. If set to \code{TRUE}, all genes that are not part of the data set (i.e. not in \code{rownames(dat)}) are removed from the gene sets.
}
\item{cluster}{
	If the analyses should be applied in parallel for the different values of \code{k}, this parameter must hold an initialized cluster as returned by \code{\link[parallel]{makeCluster}}. If this parameter is \code{NULL}, the analyses are performed sequentially.
}
}
\details{
The uncertainty analysis repeatedly replaces parts of the original gene sets by random genes and calculating the gene set statistics for these randomized gene sets. This yields a distribution of gene set statistic values for slightly modified variants of the original gene set.
}
\value{
Returns a list (of class \code{uncertaintyResult}) with the following elements:

\item{uncertainty}{
	The calculated stability of the original gene set.
}
\item{confidenceValues}{
	A matrix of quantiles of \code{gssValues} (signLevel, 0.5, 1-signLevel). One row for each value in \code{k}.
}
\item{uncertaintyEvaluations}{
	A list with one entry per value in \code{k} containing the following elements: 
	\itemize{
	\item Quantiles of \code{gssValues}: signLevel, 0.5, 1-signLevel.
	\item gssValues: A vector of gene set statistic values, one for each randomly sampled gene set.
	\item uncertainGeneSets: A matrix containing all partially random gene sets.
	\item k: The percentage of genes in the randomized gene sets taken from the original gene set.
	}
	
}
\item{signLevel}{
	The significance level used for this analysis.
}
\item{originalGeneSetValues}{
	Result of \code{geneSetAnalysis} for the original \code{geneSet}.
}
}



\seealso{
\code{\link{geneSetAnalysis}}, \code{\link{gsAnalysis}}, \code{\link{gls}}, \code{\link{transformation}}, \code{\link{gss}}, \code{\link{plot.uncertaintyResult}}
}
\examples{
data(exampleData)

  res <- evaluateGeneSetUncertainty(
  	# parameters for evaluateGeneSetUncertainty
  	dat = countdata,
  	geneSet = pathways[[1]],
  	analysis = analysis.averageCorrelation(),
  	numSamplesUncertainty = 10,
  	k = seq(0.1,0.9, by=0.1),
  	# additional parameters for analysis.averageCorrelation
  	labs = labels,
  	numSamples = 10)
}


%\keyword{}
