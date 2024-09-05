\name{summary.gsaResult}
\alias{summary.gsaResult}
\title{
Summarize gene set analysis results
}
\description{
Prints a summary of a gene set analysis result object.
}
\usage{
\method{summary}{gsaResult}(object,
		mode = c("summary", "table"),
		orderBy = c("adjustedPValues", "rawPValues", "geneSetName"),
		significantOnly = FALSE,
		signLevel = object$signLevel,
		...)
}
\value{Returns a data.frame with a summary of the produced results such as data set name, p-values, size of the given gene set.}
\arguments{
\item{object}{
	A result object as returned by \code{\link{geneSetAnalysis}}.
}
\item{mode}{
	Specifies the type of information that is displayed: By default (\code{mode="summary"}), a brief summary of the number of significant and insignificant gene sets is printed. For \code{mode="table"}, \code{\link{createSummaryTable}} is called, and a detailed table of adjusted and unadjusted p-values and the number of genes for each gene set is printed.
}

\item{orderBy}{
	If \code{mode="table"}, this specifies which field should be used for the row ordering. By default, rows are ordered according to the adjusted p-values.
}
\item{significantOnly}{
	If \code{mode="table"}, this specifies whether all gene sets (\code{significantOnly=FALSE}) or only the statistically significant gene sets (\code{significantOnly=TRUE}) should be included in the table. 
}
\item{signLevel}{
  If \code{mode="table"} and \code{significantOnly=TRUE}, this specifies the significance level for the results that should be included in the table. By default, the original significance level of the analysis is used.
}

\item{...}{
Currently unused
}
}

\seealso{
\code{\link{geneSetAnalysis}}, \code{\link{hist.gsaResult}}, \code{\link{createSummaryTable}}
}
\examples{
data(exampleData)
  # perform gene set analyses for several pathways
  res <- geneSetAnalysis(
  	# global parameters
  	dat = countdata,
  	geneSets = pathways, 
  	analysis = analysis.averageCorrelation(), 
  	# additional parameters for analysis.averageCorrelation
  	labs = labels,
  	numSamples = 10)
  
  #summarize the analyses
  summary(res, mode = "summary")
  
  summary(res, mode = "table", orderBy = "rawPValues")
}


%\keyword{ }
