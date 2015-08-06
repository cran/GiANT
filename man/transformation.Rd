\name{Transformations}
\alias{transformation}
\alias{transformation.abs}
\alias{transformation.square}
\alias{transformation.localFdr}
\alias{transformation.binarize}
\alias{transformation.rank}
\alias{transformation.adjust}
\alias{transformation.adjustAndBinarize}
\title{
Transformations
}
\description{
Functions to transform the gene-level statistic values prior to the calculation of the gene set statistics, as used in the \code{transformation} parameter of \code{gsAnalysis}. Most of the functions wrap existing R functions.
}
\usage{
transformation.abs(x)

transformation.square(x)

transformation.localFdr(x,
	statistic="pvalue",
	cutoff.method="fndr",
	pct0=0.75)

transformation.binarize(x, quant)

transformation.rank(x)

transformation.adjust(x, adjMethod = "fdr")

transformation.adjustAndBinarize(x, adjMethod = "fdr", threshold = 0.05)
}
\arguments{
\item{x}{
	A numeric vector of gene-level statistic values, one per gene. These values are calculated by the previous step (see \code{\link{gls}}).
}
\item{statistic}{
	Specifies the null model for \code{transformation.localFdr} (see \code{statistic} parameter of \code{\link[fdrtool]{fdrtool}} for possible values).
}
\item{cutoff.method}{
	Type of cut-off method used in \code{transformation.localFdr} (see \code{cutoff.method} parameter of \code{\link[fdrtool]{fdrtool}} for possible values).
}
\item{pct0}{
	Fraction of data used by \code{transformation.localFdr} if \code{cutoff.method="pct0"} (see \code{\link[fdrtool]{fdrtool}} for a detailed description).
}

\item{quant}{
	For \code{transformation.binarize}, this numeric value in the interval [0,1] defines the percentage of gene-level statistic values which should be set to zero. The remaining values are set to one.
}
\item{adjMethod}{
	The method to use for the adjustment for multiple testing (see \code{method} parameter of \code{\link{p.adjust}} for possible values).
}
\item{threshold}{
	The threshold for differential expression of a gene (defaults to \code{0.05}). Values smaller than these threshold are set to 1 , others to 0.
}
}
\details{
Standard transformation functions for gene-level statistics (to be used in an analysis pipeline defined by \code{gsAnalysis}):
\itemize{

\item{\code{transformation.abs}: }
{
Calculates the absolute values of the elements in \code{x} (a wrapper for \code{\link{abs}}).
}

\item{\code{transformation.square}: }
{
Squares all elements in \code{x}.
}

\item{\code{transformation.localFdr}: }
{
Calculates the local fdr for the elements in \code{x}. This is a wrapper for \code{\link[fdrtool]{fdrtool}}.
}

\item{\code{transformation.binarize}: }
{
Binarizes the values in \code{x} by using the \code{quant} quantile as a threshold.
}

\item{\code{transformation.rank}: }
{
Ranks the values in \code{x} and returns the rank vector.
}

\item{\code{transformation.adjust}: }
{
Adjusts for multiple testing according to the adjustment method specified in \code{adjMethod}.
}

\item{\code{transformation.adjustAndBinarize}: }
{
Adjusts for multiple testing according to the adjustment method specified in \code{adjMethod} and binarizes the resulting p-values according to \code{threshold} (values smaller than the threshold become 1 others 0).}
}
}
\value{
All functions return a vector of transformed values having the same length as \code{x}.
}

\seealso{
\code{\link{geneSetAnalysis}}, \code{\link{gsAnalysis}}, \code{\link{gss}}, \code{\link{gls}}
}

%\keyword{}
