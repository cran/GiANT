\name{plotOverrepresentation}
\alias{plotOverrepresentation}
\title{
Plot overlap of gene sets and core set
}
\description{
Plots a Venn diagramm of the overlaps of the core set and gene sets in an overrepresentation analysis.
}
\usage{
plotOverrepresentation(
	object,
	signLevel = object$signLevel,
	subset = NULL,
	aggregate = FALSE,
	ask = FALSE,
	...)
}
\value{No return value, called for creating a plot.}
\arguments{
\item{object}{
	A result of a call to \code{geneSetAnalysis} using the predefined analysis \if{latex}{\cr}\code{\link[=predefinedAnalyses]{analysis.customOverrepresentation}} or \code{\link[=predefinedAnalyses]{analysis.overrepresentation}}.
}
\item{signLevel}{
	Only results with significance level smaller than the given value are included in the venn diagram.
}
\item{subset}{
	Indices for the results that should be included in the diagram. 
}
\item{aggregate}{
	 Specifies whether all gene sets should be plotted in a single Venn diagram (which is possible for at most four gene sets) or whether there should be one Venn diagram for each gene set.
}
\item{ask}{
If set to true, the plot function will prompt for a user input for each new plot that is shown on an interactive device (see \code{\link[=par]{par("ask")}}). If \code{aggregate = TRUE}, \code{ask} is ignored.
}
\item{...}{
	Further parameters to be passed to \code{\link[limma]{vennDiagram}}.
}

}

\seealso{
\code{\link{geneSetAnalysis}}, \code{\link{predefinedAnalyses}}, \code{\link{gsAnalysis}}
}
\examples{
	data(exampleData)
  # use the absolute correlation as a gene-level statistic
  stat <- abs(apply(countdata,1,cor,y = labels))
  # define the core set as the 25% genes with the highest correlation
  coreSet <- rownames(countdata)[tail(order(stat), 25)]
  
  # perform an overrepresentation analysis
  resOverrep <- geneSetAnalysis(
  	dat = countdata,
  	geneSets = pathways,
  	analysis = analysis.customOverrepresentation(),
  	coreSet = coreSet,
  	adjustmentMethod = "fdr")
  
  # plot a Venn diagram
  plotOverrepresentation(resOverrep, subset = 1:3, aggregate = TRUE)
}


%\keyword{}
