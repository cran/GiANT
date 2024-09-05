\name{createSummaryTable}
\alias{createSummaryTable}
\title{
Create an overview table for an analysis
}
\description{
Creates a data frame summarizing an analysis. This table has one row per gene set, each comprising the adjusted and unadjusted p-values and the number of genes for the set.
}
\usage{
createSummaryTable(object,
		orderBy = c("adjustedPValues", "rawPValues", "geneSetName"),
		significantOnly = FALSE,
		signLevel = object$signLevel)
}

\arguments{
\item{object}{
	A result object as returned by \code{\link{geneSetAnalysis}}.
}

\item{orderBy}{
	Specifies which field should be used for the row ordering. By default, rows are ordered according to the adjusted p-values.
}
\item{significantOnly}{
	Specifies whether all gene sets (\code{significantOnly=FALSE}) or only the statistically significant gene sets (\code{significantOnly=TRUE}) should be included in the table. 
}
\item{signLevel}{
  If \code{significantOnly=TRUE}, this specifies the significance level for the results that should be included in the table. By default, the original significance level of the analysis is used.
}
}

\value{
A data frame with one row for each included gene set and the columns "adjustedPValues", "rawPValues", "geneSetName" and "geneSetSize". For overrepresentation analyses, there is an additional column "intersectSize" specifying the size of the intersection of the core set and the corresponding gene set.
}

\seealso{
\code{\link{geneSetAnalysis}}, \code{\link{hist.gsaResult}}, \code{\link[=summary.gsaResult]{summary}}
}
\examples{
# load data
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

  tab <- createSummaryTable(res)
}


%\keyword{ }
