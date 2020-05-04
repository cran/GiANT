\name{Gene set analysis}
\alias{geneSetAnalysis}
\title{
Main interface for enrichment analyses.
}
\description{
The main function of the package that performs a gene set analysis for a list of gene sets.
}
\usage{
geneSetAnalysis(
	...,
	dat,
	geneSets,
	analysis,
	signLevel = 0.05,
	preprocessGeneSets = FALSE,
	adjustmentMethod = p.adjust.methods,
	cluster = NULL)
}

\arguments{
\item{...}{
  Additional parameters for the different steps of the analysis pipeline, depending on the concrete configuration supplied in \code{analysis}.
}
\item{dat}{
	A numeric matrix of gene expression values for all analyzed genes. Here, each row corresponds to one gene, and each column corresponds to one sample. The rows must be named with the gene names used in the gene sets.
}
\item{geneSets}{
	A list of gene sets, where each gene set is a vector of gene names corresponding to the row names of \code{dat}.
}
\item{analysis}{
	An object of type \code{gsAnalysis} as returned by \code{\link{gsAnalysis}} or by the predefined configurations (see \code{\link{predefinedAnalyses}}).
}
\item{signLevel}{
	The significance level for the significance assessment of the gene sets (defaults to \code{0.05}). 
}
\item{preprocessGeneSets}{
	Specifies whether the gene sets in \code{geneSets} should be preprocessed or not. If set to \code{TRUE}, all genes that are not part of the data set (i.e. not in \code{rownames(dat)}) are removed from the gene sets.
}
\item{adjustmentMethod}{
	The method to use for the adjustment for multiple testing (see \code{method} parameter of \code{\link{p.adjust}} for possible values).
}
\item{cluster}{
  If the analyses should be applied in parallel for the gene sets, this parameter must hold an initialized cluster as returned by \code{\link[parallel]{makeCluster}}. If this parameter is \code{NULL}, the analyses are performed sequentially.
}
}

\details{
This is the main interface function of the package for gene set enrichment analyses. Analyses usually consist of a pipeline of steps. Often, the first step is the calculation of a summary statistic for the relation of each gene to the class labels. These values or transformations thereof are employed to calculate a gene set statistic for each of the supplied gene sets. The significance of gene set enrichments can be determined according to different methods, and the robustness of gene sets can be evaluated by slightly modifying the gene sets. To provide a flexible mechanism for the plethora of different approaches arising from the different choices, basic pipeline configurations are encapsulated in \code{gsAnalysis} objects which can be created using the \code{\link{gsAnalysis}} function. Ready-to-use configuration objects for certain well-known methods are included in the package (see \code{\link{predefinedAnalyses}}). Parameters of the chosen analysis pipeline can be set in the \code{...} parameter.
}
\value{
An object of the type \code{gsaResult} with the following elements:

\item{adjustedPValues}{A vector of p-values, one for each gene set. These values are already adjusted for multiple testing according to the \code{adjustmentMethod} parameter.
}

\item{rawPValues}{The raw unadjusted p-values, one for each gene set.
}

\item{res.all}{
	A list comprising the detailed results for each gene set. Each element of this list is another list with the following components:

	\code{pValue}: The raw (unadjusted) p-value for the gene set.

	\code{geneSetValues}: If \code{analysis} is a global analysis, this is the object returned by the method for the corresponding gene set. For an analysis pipeline, this holds the values of the gene-level statistic, the transformed values and the values of the gene set statistic (see also \code{\link{gsAnalysis}}).

	\code{significanceValues}: Gene set statistics for each randomly drawn gene set for significance assessment and a \code{list} of this gene sets. Only set for analysis of type 'geneSetAnalysis'. \code{NULL} for 'global' analysis.

	\code{geneSet}: The supplied gene set.
}

\item{signLevel}{The significance level used for this analysis.}

\item{analysis}{The performed analysis (of type \code{gsAnalysis}).}

\item{analysisType}{A character string identifying the analysis as an enrichment analysis pipeline (\code{"geneSetAnalysis"}) or as a global analysis (\code{"global"}).}

\item{adjustmentMethod}{The method used to adjust the p-values in \code{adjustedPValues}}

}
\references{
Ackermann, M., Strimmer, K. (2009) A general modular framework for gene set enrichment analysis. \emph{BMC Bioinformatics}, \bold{10}(1), 47.
}
%%\note{
%%---
%%  ~~further notes~~
%%}
%\section{Warning}{
%Methods in the package can be freely combined and most of the combinations will lead to any kind of result. The authors do not assume any responsibility whether these results are meaningful / correct or not. This package do not replace the expert knowledge about the used methods it just simplifies the application of them.
%}

\seealso{
\code{\link{gsAnalysis}}, \code{\link{gls}}, \code{\link{transformation}}, \code{\link{gss}}, \code{\link{global}}, \code{\link{significance}}, \code{\link{evaluateGeneSetUncertainty}}, \code{\link{hist.gsaResult}}, \code{\link{preprocessGs}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# load data
require(GlobalAncova)
data(vantVeer)
data(phenodata)
data(pathways)

# apply predefined analysis for gene set enrichment analysis
res <- geneSetAnalysis(	
	# parameters for geneSetAnalysis
	dat = vantVeer,
	geneSets = pathways[1],
	analysis = analysis.averageCorrelation(),
	adjustmentMethod = "fdr",
	# additional parameters for analysis.averageCorrelation
	labs = phenodata$metastases,
	method = "pearson",
	numSamples = 10)
}

%\keyword{}
