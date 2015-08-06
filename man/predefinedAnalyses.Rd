\name{predefinedAnalyses}
\alias{predefinedAnalyses}
\alias{analysis.gsea}
\alias{analysis.overrepresentation}
\alias{analysis.customOverrepresentation}
\alias{analysis.averageCorrelation}
\alias{analysis.averageTStatistic}
\alias{analysis.globalTest}
\alias{analysis.globalAncova}
\title{
Predefined enrichment analyses
}
\description{
Predefined analysis configurations that can be used in \code{\link{geneSetAnalysis}}
}
\usage{
analysis.gsea()
analysis.overrepresentation()
analysis.customOverrepresentation()
analysis.averageCorrelation()
analysis.averageTStatistic()
analysis.globalTest()
analysis.globalAncova()
}

\details{

The above functions return configurations for state-of-the-art analysis pipelines that can be used in \code{\link{geneSetAnalysis}}. All configurations are preconfigured collections of standard methods for the different pipeline steps. The following lists the methods chosen for the different steps and their parameters. For more detailed descriptions of these methods, please refer to the linked manual pages.

\itemize{
\item{\code{analysis.gsea} }{defines the Gene Set Enrichment Analysis (GSEA) method by Subramanain et al. 
Here, the gene-level statistic the absolute correlation calculated by \code{\link{gls.cor}} with the associated parameters \code{labs, method} and a preprocessing by \code{\link{transformation.abs}}. As a gene set statistic, the enrichment score (function \code{\link{gss.enrichmentScore}} with parameter \code{p}) is calculated. The significance is assessed in a permutation test using \code{\link{significance.permutation}} with \code{testAlternative = "greater"} and free parameter \code{numSamples, labs}.}

\item{\code{analysis.overrepresentation} }{calculates an overrepresentation analysis using the gene-level statistic \code{\link{gls.tStatistic}} with parameters \code{pValue}(should be \code{TRUE}), \code{alternative} and \code{labs}. The resulting values are then transformed via \code{transformation.adjustAndBinarize} (parameters are the \code{adjMethod} and \code{threshold}). Finally \code{gss.fisherExactTest} is used as gene set statistic.
}

\item{\code{analysis.customOverrepresentation} }{calculates an overrepresentation analysis using a user-defined core set \code{coreSet}. That is, instead of calculating this core set internally based on differential expression as the standard overrepresentation analysis, this function allows for defining custom core sets. It internally uses the global analysis \code{\link{global.overrepresentation}}.
}

\item{\code{analysis.averageCorrelation} }{calculates the gene-level statistic as the absolute correlation using \code{\link{gls.cor}} (with parameters \code{labs, method}) and \code{\link{transformation.abs}}. The gene set statistic is the mean correlation calculated by \code{\link{gss.mean}}. The significance is assessed by comparing the gene set statistic to randomly sampled gene sets using \code{\link{significance.sampling}} (with the parameter \code{numSamples} and the preset parameter \code{testAlternative = "greater"}).
}

\item{\code{analysis.averageTStatistic} }{uses the absolute t statistic as the gene-level statistic by applying \code{\link{gls.tStatistic}} (with parameters \code{labs, pValue, alternative}) and \if{latex}{\cr}\code{\link{transformation.abs}}. The gene set statistic is the mean t statistic in the gene set as returned by \code{\link{gss.mean}}. The significance is assessed by comparing the gene set statistic to randomly sampled gene sets using \code{\link{significance.sampling}} (with the parameter \code{numSamples} and the preset parameter \code{testAlternative = "greater"}).
}

\item{\code{analysis.globalTest} }{performs a global gene set enrichment analysis by Goeman et al. by applying the \code{\link{global.test}} function which in turn wraps the \code{\link[globaltest]{gt}} function in the \pkg{globaltest} package.}

\item{\code{analysis.globalAncova} }{applies the global ANCOVA method by Hummel et al. using the global method \code{\link{global.ancova}} which wraps the \code{\link[GlobalAncova]{GlobalAncova}} function in the \pkg{GlobalAncova} package.}

}
}
\value{
All functions return an object of class \code{gsAnalysis} that specifies the corresponding analysis parameters for \code{\link{geneSetAnalysis}}.
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S., Mesirov, J. P. (2005) Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. \emph{Proceedings of the National Academy of Science of the United States of America}, \bold{102}, 15545--15550.

Hummel, M., Meister, R., Mansmann, U. (2008) GlobalANCOVA: exploration and assessment of gene group effects. \emph{Bioinformatics}, \bold{24}(1), 78--85.

Goeman, J. J., van de Geer, S. A., de Kort, F., van Houwelingen, H. C. (2004) A global test for groups of genes: testing association with a clinical outcome. \emph{Bioinformatics}, \bold{20}(1), 93--99.
}

\seealso{
\code{\link{geneSetAnalysis}}, \code{\link{gsAnalysis}}
}
\examples{
# load data
require(GlobalAncova)
data(vantVeer)
data(phenodata)
data(pathways)

# apply a gene set analysis based on the average absolute correlation
resAvCor <- geneSetAnalysis(
	# parameters for geneSetAnalysis
	dat = vantVeer,
	geneSets = pathways[1:2],
	analysis = analysis.averageCorrelation(),
	adjustmentMethod = "fdr",
	# additional parameters for analysis.averageCorrelation
  labs = phenodata$metastases,
	method = "pearson",
	numSamples = 100)

# apply an overrepresentation analysis
resOverrep <- geneSetAnalysis(
	# parameters for geneSetAnalysis
	dat = vantVeer,
	geneSets = pathways,
	analysis = analysis.overrepresentation(),
	adjustmentMethod = "fdr",
	# additional parameters for analysis.overrepresentation
	pValue = TRUE,
	threshold = 0.1,
	labs = phenodata$metastases
)

# apply a global analysis using GlobalAncova
resGA <- geneSetAnalysis(
	# parameters for geneSetAnalysis
	dat = vantVeer,
	geneSets = pathways[1:2],
	analysis = analysis.globalAncova(),
	adjustmentMethod = "fdr",
	# additional parameters for analysis.globalAncova
	labs = phenodata$metastases,
	method = "approx")
}

%\keyword{}
