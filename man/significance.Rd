\name{SignificanceAssessment}
\alias{significance}
\alias{significance.sampling}
\alias{significance.permutation}
\alias{significance.restandardization}
\title{
Significance assessment
}
\description{
Functions to assess the significance of the gene-level statistics, as used in the \code{significance} parameter of \code{\link{gsAnalysis}}. These functions are based on applying the same analysis to randomly modified data sets or gene sets and comparing their statistic values to the original gene set statistic value.
}
\usage{
significance.sampling(
	...,
	dat,
	geneSet,
	analysis,
	glsValues,
	numSamples = 1000)

significance.permutation(
	...,
	dat,
	geneSet,
	analysis,
	glsValues,
	numSamples = 1000,
	labs)

significance.restandardization(
	...,
	dat,
	geneSet,
	analysis,
	glsValues,
	numSamples = 1000,
	labs)
}

\arguments{
\item{...}{
  Additional parameters for the different steps of the analysis pipeline, depending on the concrete configuration supplied in \code{analysis}.
}
\item{dat}{
	The original data set as a numeric matrix of gene expression values for all analyzed genes. Here, each row corresponds to one gene, and each column corresponds to one sample. The rows must be named with the gene names used in the gene sets.
}
\item{geneSet}{
  The original gene set in form of a vector of gene names corresponding to the row names of \code{dat}.
}
\item{analysis}{
  The analysis applied to the original gene set (that should also be applied to the modified gene sets). This is an object of type \code{gsAnalysis} as produced by the function \code{\link{gsAnalysis}}.
}
\item{glsValues}{
	A vector containing the (possibly transformed) gene-level statistic values for each gene in the original data set \code{dat}.
}
\item{numSamples}{
	The number of random samples that should be taked to calculate the null distribution for the significance assessment. Default is \code{1000} for each test.
}
\item{labs}{
	A vector of class labels for the samples in \code{dat} for \code{significance.permutation} and \code{significance.restandardization}.
}
}

\details{
Standard methods for the significance assessment of a gene set statistic (to be used in an analysis pipeline defined by \code{gsAnalysis}):

\itemize{
\item{\code{significance.sampling}: }
{
This function repeatedly draws random gene sets. Their gene set statistic values form the null distribution.
}

\item{\code{significance.permutation}: }
{
This function repeatedly permutes the labels of the data set. The gene set statistic values for the original gene set on the permuted data set form the null distribution.
}

\item{\code{significance.restandardization}: }
{
This function applies both a gene set sampling and a label permutation. The permutation statistic values are standardized by their mean and standard deviation and then restandardized based on the gene set sampling statistic values. These restandardized values form the null distribution (Efron and Tibshirani).
}

}
}
\value{
\code{significance.sampling} returns a list with the following elements:
\item{gssValues}{A vector of gene set statistic values, one entry per sample.}
\item{randomGeneSets}{A matrix containing the gene sets which were sampled randomly from the set of all genes.}

\code{significance.permutation} returns a list with the following elements:
\item{gssValues}{A vector of gene set statistics, one entry per sample.}
\item{permutations}{A matrix, where each column contains the indices of one permutation.}

\code{significance.restandardization} returns a list with the following elements:
\item{gssValues}{A vector of gene set statistics, one entry per sample.}
\item{samplingValues}{A list of sub-lists, each containing one sampling result as defined above.}
\item{permutationValues}{A list of sub-lists, each containing one permutation result as defined above.}
}
\references{
Efron, B., Tibshirani, R. (2007) On testing the significance of sets of genes. \emph{Annals of Applied Statistics}, \bold{1}, 107-129.
}
\seealso{
\code{\link{geneSetAnalysis}}, \code{\link{gsAnalysis}}, \code{\link{hist.gsaResult}}
}

%\keyword{}
