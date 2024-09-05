\name{plot}
\alias{plot.uncertaintyResult}
\title{
Plots the results of an uncertainty analysis.
}
\description{
For each percentage of original gene set genes, the quantiles of the distribution obtained by a resampling simulation are plotted. Significance threshold (quantile of the Null distribution) and the test statistic of the original gene set are drawn as horizontal lines.
}
\usage{
\method{plot}{uncertaintyResult}(x,
	signLevel = x$signLevel,
	addLegend = TRUE,
	addMinimalStability = FALSE,
	...)
}
\value{No return value, called for creating a plot.}
\arguments{
\item{x}{
	A result of a call to \code{evaluateGeneSetUncertainty} (see also Details). 
}

\item{signLevel}{
	Only results with significance level smaller than the given value are plotted.
}
\item{addLegend}{
	If set to true (default), a \code{\link{legend}} is added to the plot.
}
\item{addMinimalStability}{
	If set to true, a line is added to the plot giving the minimal stability.
}
\item{...}{
	Other parameters which can be used for histograms (see \code{\link[graphics]{plot}}).
}
}
\details{
The function plots the quantiles of the resampling distributions for evaluated degrees of fuzziness. It requires the 
significance assessment step of the enrichment analysis configuration (parameter \code{significance} or \code{\link{gsAnalysis}}) to be a computer-intensive testing procedure that yields a distribution of gene set statistic values under the null hypothesis. Predefined configurations for which this plot works are \code{\link{analysis.gsea}}, \code{\link{analysis.averageCorrelation}} and \code{\link{analysis.averageTStatistic}}.

Three lines, corresponding to the different qunatiles with one dot per fuzziness evaluation (\code{k}) are plotted for the analysis in \code{x}. The significance threshold is shown as a green horizontal line. The statistic value of the original input set is depicted as a red horizontal line.

If \code{addMinimalStability} is \code{TRUE}, the lower bound of the stability is ploted as a dotted line.
}

\seealso{
\code{\link{geneSetAnalysis}}, \code{\link{predefinedAnalyses}}, \code{\link{gsAnalysis}}, \code{\link{evaluateGeneSetUncertainty}}
}
\examples{
# load data
data(exampleData)

  res <- evaluateGeneSetUncertainty(
  	# parameters for evaluateGeneSetUncertainty
  	dat = countdata,
  	geneSet = pathways[[1]],
  	analysis = analysis.averageCorrelation(),
  	numSamplesUncertainty = 10,
  	N = seq(0.1,0.9, by=0.1),
  	# additional parameters for analysis.averageCorrelation
  	labs = labels,
  	numSamples = 10)
  
  # plot the results for the cell cycle control gene set
  plot(res, addMinimalStability = TRUE)
}


%\keyword{}
