\name{GeneLevelStatistics}
\alias{gls}
\alias{gls.cor}
\alias{gls.regression}
\alias{gls.foldChange}
\alias{gls.tStatistic}
\alias{gls.moderateTStatistic}
\alias{gls.nBinomTest}
\title{
Gene-level statistics
}
\description{
Functions to calculate the gene-level statistic, as used in the \code{gls} parameter of \code{gsAnalysis}.
A gene-level statistic calculates a measure of correlation between the expression of a gene and the class labels.
}
\usage{
gls.cor(dat, labs, method = "pearson")

gls.regression(dat, labs)

gls.foldChange(dat, labs, logMeasurements = TRUE)

gls.tStatistic(dat, labs, pValue = FALSE, alternative = "two.sided")

gls.moderateTStatistic(dat,labs)

gls.nBinomTest(dat, labs,
	returnValue = c("pval", "qval", "foldChange", "log2FoldChange"),
	dispersionMethod = "blind",
	dispersionSharingMode = "fit-only",
	dispersionFitType = "local")
}

\arguments{
\item{dat}{
	A numeric matrix of gene expression values for all analyzed genes. Here, each row corresponds to one gene, and each column corresponds to one sample. The rows must be named with the gene names used in the gene sets.
}
\item{labs}{
	A vector of class labels for the samples in \code{dat}. 
}
\item{logMeasurements}{
	For \code{gls.foldChange}, whether the values in \code{dat} are logarithmized (\code{logMeasurements=TRUE}) or not (\code{logMeasurements=FALSE}).
}
\item{method}{
	For \code{gls.cor}, the correlation method to be used (see \code{\link{cor}}).
}
\item{pValue}{
	For \code{gls.tStatistic}, this specifies whether the p-value (\code{pValue=TRUE}) or the test statistic (\code{pValue=FALSE}) of the t test should be returned.
}
\item{alternative}{
	For \code{gls.tStatistic}, this specifies the alternative of the t-test. See also \code{\link{t.test}}.
}
\item{returnValue}{
	For \code{gls.nBinomTest}, this determines the type of values values that should be returned. \code{"pval"} returns the raw p-values, \code{"qval"} returns the p-values adjusted by the FDR, \code{"foldChange"} returns the fold changes, and \code{"log2FoldChange"} returns the log2 fold changes. For more details, see \code{\link[DESeq2]{results}}.
}
\item{dispersionMethod}{
	For \code{gls.nBinomTest}, this specifies how the empirical dispersion is computed (see \code{\link[DESeq2]{estimateDispersions}}).
}
\item{dispersionSharingMode}{
	For \code{gls.nBinomTest}, this specifies which values should be used by \code{\link[DESeq2]{results}} (fitted values or empirical values, see \code{\link[DESeq2]{estimateDispersions}} for more details).
}
\item{dispersionFitType}{
	For \code{gls.nBinomTest}, this determines the method for fitting the dispersion-mean relation (see \code{\link[DESeq2]{estimateDispersions}}).
}
}
\details{
Standard functions for the calculation of gene-level statistics (to be used in an analysis pipeline defined by \code{gsAnalysis}): 
\itemize{

\item{\code{gls.cor}: }
{
Calculates the correlation of the gene expression values to the class labels.
}

\item{\code{gls.regression}: }
{
Calculates the slope of a linear regression of the gene expression values and the class labels.
}

\item{\code{gls.foldChange}: }
{
Calculates the (standard or log2) fold change between the measurements for the two classes.
}

\item{\code{gls.tStatistic}: }
{
Calculates the p-value or the statistic of a two-sample t test for the measurements of the two classes
}

\item{\code{gls.moderateTStatistic}: }
{
Calculates the moderate t statistic for the measurements of the two classes
}

\item{\code{gls.nBinomTest}: }
{
Applies the negative binomial test for sequencing data based on the \pkg{DESeq2} package to test for differences between two classes (see \code{\link[DESeq2]{results}}).
}

}
}
\value{
Each of these function returns a numeric vector of gene-level statistics (one entry per gene).
}

\seealso{
\code{\link{geneSetAnalysis}}, \code{\link{gsAnalysis}}, \code{\link{gss}}, \code{\link{transformation}}
}

%\keyword{}
