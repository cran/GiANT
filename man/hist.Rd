\name{hist}
\alias{hist.gsaResult}
\title{
Null distribution histogram and statistic of the input set for enrichment analyses.
}
\description{
Plots the distribution of gene set statistic values obtained in different resampling settings of an enrichment analysis, and draws the statistic value of the input set as a vertical line.
}
\usage{
\method{hist}{gsaResult}(x,
	signLevel = x$signLevel,
	subset = NULL,
	ask = FALSE,
	addLegend = TRUE,
	...)
}

\arguments{
\item{x}{
	A result of a call to \code{geneSetAnalysis} (see also Details). 
}

\item{signLevel}{
	The significance level that should be applied for the plots. Default is the significance level used for the analysis in \code{x}.
}
\item{subset}{
	Indices for the results that should be included in the diagram. 
}
\item{ask}{
	If set to true, the plot function will prompt for a user input for each new plot that is shown on an interactive device (see \code{\link[=par]{par("ask")}}). 
}
\item{addLegend}{
	If set to true (default), a \code{\link{legend}} is added to the plot.
}
\item{...}{
	Other parameters which can be used for histograms (see \code{\link[graphics]{hist}}).
}
}
\details{
The function plots the distribution of gene set statistic values under the null hypothesis. It requires the 
significance assessment step of the enrichment analysis configuration (parameter \code{significance} or \code{\link{gsAnalysis}}) to be a computer-intensive testing procedure that yields a distribution of gene set statistic p-values under the null hypothesis. Predefined configurations for which this plot works are \code{\link{analysis.gsea}}, \code{\link{analysis.averageCorrelation}} and \code{\link{analysis.averageTStatistic}}.

A histogram is plotted for the analysis in \code{x}. If \code{x} includes the analyses for several gene sets, one histogram is plotted for each of the gene sets.

The statistic value of the input set is depicted as a vertical line.

The most common graphical parameters can be supplied as vectors (one entry per analyzed gene set) to vary them between the different analyses. These parameters are: \code{main}, \code{xlab}, \code{ylab}.
}

\seealso{
\code{\link{geneSetAnalysis}}, \code{\link{predefinedAnalyses}}, \code{\link{gsAnalysis}}, \code{\link{evaluateGeneSetUncertainty}}, \code{\link{plot.uncertaintyResult}}
}
\examples{
# load data
require(GlobalAncova)
data(vantVeer)
data(phenodata)
data(pathways)

res <- geneSetAnalysis(
	# global parameters
	dat = vantVeer,
	geneSets = pathways[3], 
	analysis = analysis.averageCorrelation(), 
	# additional parameters for analysis.averageCorrelation
	labs = phenodata$metastases,
	p = 1,
	numSamples = 10)

# plot the histogram for the cell cycle control gene set
hist(res, main = names(pathways[3]))
}

%\keyword{}
