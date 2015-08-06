\name{gsAnalysis}
\alias{gsAnalysis}
\title{
Gene set analysis.
}
\description{
Defines the configuration of an analysis that can be performed using \code{\link{geneSetAnalysis}}, and returns it as a wrapper object.
}
\usage{
gsAnalysis(name,
	gls = NULL,
	glsParameterNames = NULL,
	transformation = NULL,
	transformationParameterNames = NULL,
	gss = NULL,
	gssParameterNames = NULL,
	globalStat = NULL,
	globalStatParameterNames = NULL,
	significance = NULL,
	significanceParameterNames = NULL,
	testAlternative = c("greater", "less"))

}

\arguments{
\item{name}{
	A character string describing the analysis.
}
\item{gls}{
	The name of the function that calculates the gene-level statistic for a given dataset. If set to \code{NULL}, it is assumed that the input data already comprises gene-level statistic values, and the input is directly supplied to \code{transformation}. The first (fixed) parameter of a gls function is the dataset.
}
\item{glsParameterNames}{
	A character vector of names of the parameters used by the gene-level statistic defined in \code{gls}.
}
\item{transformation}{
	The name of the function that transforms the gene-level statistics values. If set to \code{NULL}, the values supplied by \code{gls} are directly handed over to \code{gss}. The only fixed parameter for transformation is the gene-level statistic (supplied as first parameter).
}
\item{transformationParameterNames}{
  	A character vector of names of the parameters used by the transformation defined in \code{transformation}.
}
\item{gss}{
	The name of the function that calculates the gene set statistics from untransformed or transformed gene-level statistic values. If set to \code{NULL}, the values supplied by \code{transformation} are directly handed over to \code{significance}. Fixed parameters are the transformed values (first parameter) and \code{geneSetIndices} containing the (row-) indices of the current gene set genes in the dataset.
}
\item{gssParameterNames}{
  	A character vector vector of names of the parameters used by the gene set statistic defined in \code{gss}.
}
\item{globalStat}{
	If the gene set analysis consists of a global test and cannot be divided into the calculation of a gene-level statistic, a transformation and the calculation of the gene set statistics, this parameter defines the name of the function that performs this global analysis. In this case, the parameters \code{gls}, \code{transformation} and \code{gss} are ignored. Fixed parameters are \code{dat} containing the whole dataset and \code{geneSet} containing the current gene set. 
}
\item{globalStatParameterNames}{
  	A character vector names of the parameters used by the global analysis defined in \code{globalStat}.
}
\item{significance}{
	The name of a method that performs a significance assessment for the gene set statistic values. If set to \code{NULL}, \code{geneSetAnalysis} does not return p-values, but returns the statistics supplied by \code{gss} or \code{globalStat}. Fixed parameters are \code{dat} containing the whole dataset, \code{geneSet} containing the current gene set, \code{analysis} with the supplied \code{gsAnalysis} and \code{glsValues} with (depending on whether a transformation is supplied or not) transformed gene-level statistics for each gene in the dataset.
}
\item{significanceParameterNames}{
  	A character vector of names of the parameters used by the significance assessment method defined in \code{significance}.
}
\item{testAlternative}{
	Specifies the alternative hypothesis of the significance test for the gene set enrichment, which may be dependent on the chosen gene set statistic. Must be one of \code{"greater"} or \code{"less"}.
}
}

\details{
The function provides a way of flexibly defining the steps of the gene set analysis pipeline. This pipeline consists of a subset of the following steps, each of which may have specific parameters:
\itemize{

\item{Gene-level analysis: }{
 A gene-level statistic scores the relationship between the measurements for a specific gene and the class labels. Typical measures include correlation coefficients, the t statistic or the fold change between the groups (see \code{\link{gls}} for gene-level statistics included in the package). 
  }

\item{Transformation of gene-level statistics: }{
  Optionally, the gene-level statistic values can be postprocessed, e.g. by taking the absolute value or the square for correlation values or by binarizing or ranking values. See \code{\link{transformation}} for transformations included in the package.
}

\item{Gene set analysis: }
{
Based on the (possibly transformed) gene-level statistics, the gene set(s) of interest is/are scored. Examples are the median, the mean or the enrichment score of the gene-level statistic values in the gene set(s). See \code{\link{gss}} for gene set statistics included in the package.
}

\item{Significance assessment: }
{
To assess the significance of the gene set statistic value(s) with respect to a null distribution, computer-intensive tests are performed. These tests repeatedly sample random label vectors or gene sets and calculate their gene set statistic values. These values can then be compared to the true gene set statistics. See \code{\link{significance}} for significance assessment methods included in the package.
}

\item{Global analysis: }
{
As an alternative to the above pipeline steps, it is possible to define a single, global method that directly calculates an enrichment p-value for a supplied data set and gene set. See \code{\link{global}} for the global analysis tests included in the package.
}
}

Several state-of-the-art analyses have predefined configuration objects in which the above steps are defined accordingly (see \code{\link{predefinedAnalyses}}).
}
\value{
An object of class \code{gsAnalysis} with components corresponding to the above parameters.
}

\seealso{
\code{\link{predefinedAnalyses}}, \code{\link{geneSetAnalysis}}, \code{\link{evaluateGeneSetUncertainty}}, \code{\link{gls}}, \code{\link{transformation}}, \code{\link{gss}}, \code{\link{global}}, \code{\link{significance}}
}
\examples{
# defines an analysis that corresponds to gsAna1()
gsa <- gsAnalysis(
	name = "averageCorrelation",
	gls = "gls.cor",
	glsParameterNames = c("labs","method"),
	transformation = "transformation.abs",
	transformationParameterNames = NULL,
	gss = "gss.mean",
	gssParameterNames = NULL,
	globalStat = NULL,
	globalStatParameterNames = NULL,
	significance = "significance.sampling",
	significanceParameterNames = c("numSamples"),
	testAlternative = "greater")
print(gsa)

# load data
require(GlobalAncova)
data(vantVeer)
data(phenodata)
data(pathways)

# apply the previously defined analysis
res <- geneSetAnalysis(
	# global parameters
	dat = vantVeer,
	geneSets = pathways[1], 
	analysis = gsa, 
	# parameters for the specific analysis gsAna1
	labs = phenodata$metastases,
	numSamples = 100)

}

%\keyword{}
