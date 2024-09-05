\name{GiANT-package}
\alias{GiANT}
\docType{package}
\title{
Enrichment analysis
}
\description{
	Toolbox for gene set analysis of uncertain gene sets.
}
\details{
\tabular{ll}{
Package: \tab GiANT\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2020-04-29\cr
License: \tab Artistic-2.0 \cr
LazyLoad: \tab yes\cr
}

This package provides an approach for evaluating the fuzziness of a gene set. This is done by repeatedly performing gene set analyses on slightly modified versions of the gene set and comparing their enrichment scores. A utility for such uncertainty tests is provided in the \code{\link{evaluateGeneSetUncertainty}} function.

The package also comprises a generic framework for different types of enrichment analyses (Ackermann and Strimmer). It establishes a customizeable pipeline that typically consists of the following steps:
\itemize{
  \item Calculation of gene-level statistics: A gene-level statistic scores the relationship between the measurements for a specific gene and the class labels. Typical measures include correlation coefficients, the t statistic or the fold change between the groups (see \code{\link{gls}} for gene-level statistics included in the package). 
  
  \item Transformation of gene-level statistic values: Optionally, the gene-level statistic values can be postprocessed, e.g. by taking the absolute value or the square for correlation values or by binarizing or ranking values. \if{latex}{\cr}See \code{\link{transformation}} for transformations included in the package.
  
  \item Calculation of gene set statistics: Based on the (possibly transformed) gene-level statistics, the gene set(s) of interest is/are scored. Examples are the median, the mean or the enrichment score of the gene-level statistic values in the gene set(s). See \code{\link{gss}} for gene set statistics included in the package.

  \item Significance assessment: To assess the significance of the gene set statistic value(s) with respect to a null distribution, computer-intensive tests are performed. These tests repeatedly sample random label vectors or gene sets and calculate their gene set statistic values. These values can then be compared to the true gene set statistics. See \code{\link{significance}} for significance assessment methods included in the package.
}
The package represents such analysis pipelines as configuration objects that can be created using the function \code{\link{gsAnalysis}}. For predefined state-of-the-art methods, such as Gene Set Enrichment Analysis (Subramanian et al), Overrepresentation Analysis or Global Ancova (Hummel et al), it provides predefined configurations (see \code{\link{predefinedAnalyses}}).

The main function for standard gene set analyses, \code{\link{geneSetAnalysis}}, performs enrichment analyses based on pipeline configuration objects.
}
\author{
Florian Schmid, Christoph Müssel, Johann M. Kraus, Hans A. Kestler

Maintainer: Hans A. Kestler <hans.kestler@uni-ulm.de>
}
\references{
Ackermann, M., Strimmer, K. (2009) A general modular framework for gene set enrichment analysis. \emph{BMC Bioinformatics}, \bold{10}(1), 47.

Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S., Mesirov, J. P. (2005) Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. \emph{Proceedings of the National Academy of Science of the United States of America}, \bold{102}, 15545-15550.

Hummel, M., Meister, R., Mansmann, U. (2008) Globalancova: exploration and assessment of gene group effects. \emph{Bioinformatics}, \bold{24}(1), 78--85.

}

%\keyword{}

\examples{
\donttest{
  data(exampleData)
  ##################################
  # Example 1: gene set analysis   #
  ##################################
  res <- geneSetAnalysis(  
    # parameters for geneSetAnalysis
    dat = countdata,
    geneSets = pathways[1],
    analysis = analysis.averageCorrelation(),
    adjustmentMethod = "fdr",
    # additional parameters for analysis.averageCorrelation
    labs = labels,
    method = "pearson",
    numSamples = 50)
  
  summary(res, mode="table")
  
  ####################################
  # Example 2: uncertainty analysis  #
  ####################################
  resUncertainty <- evaluateGeneSetUncertainty(  
    # parameters for evaluateGeneSetUncertainty
    dat = countdata,
    geneSet = pathways[[3]],
    analysis = analysis.averageCorrelation(),
    numSamplesUncertainty = 5,
    blockSize = 1,
    k = seq(0.1,0.9,by=0.1),
    # additional parameters for analysis.averageCorrelation
    labs = labels,
    numSamples = 5)
  
  plot(resUncertainty, main = names(pathways[3]))
}
}
