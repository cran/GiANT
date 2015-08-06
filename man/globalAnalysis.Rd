\name{GlobalAnalysis}
\alias{global}
\alias{global.overrepresentation}
\alias{global.ancova}
\alias{global.test}
\title{
Global analyses
}

\description{
Functions to perform global gene set analyses, as used in the \code{globalStat} parameter of \code{\link{gsAnalysis}}.
}
\usage{
global.overrepresentation(dat,
	geneSet,
	coreSet)

global.ancova(dat,
	geneSet,
	labs,
	...)

global.test(dat,
	geneSet,
	labs,
	...)
}

\arguments{
\item{dat}{
	A numeric matrix of gene expression values for all analyzed genes. Here, each row corresponds to one gene, and each column corresponds to one sample. The rows must be named with the gene names used in the gene sets.
}
\item{geneSet}{
  A gene set in form of a vector of gene names corresponding to the row names of \code{dat}.
}
\item{coreSet}{
  A gene set of interest resulting from an analysis of \code{dat} that should be compared to \code{geneSet} in the overrepresentation analysis. This is also a vector of gene names corresponding to the row names of \code{dat}.
}
\item{labs}{ 
  A vector of class labels for the samples in \code{dat}. 
}
\item{...}{ 
  Further parameters for  \code{\link[GlobalAncova]{GlobalAncova}} and \code{\link[globaltest]{gt}}, as defined in the corresponding manual pages. The parameters \code{xx}, \code{test.genes} and \code{group} are set automatically by \code{global.ancova}, and the parameters \code{alternative}, \code{subsets} and \code{response} are set automatically by \code{global.test}.
}
}

\details{
Wrapper functions for global gene set analyses.

\itemize{
\item{\code{global.overrepresentation}: }
{
	This function performs an overrepresentation analysis by rating the overlap of \code{geneSet} and \code{coreSet} with respect to the set of all genes using Fisher's exact test.
}
\item{\code{global.ancova}: }
{
	This function performs a global gene set enrichment analysis using the global ANCOVA method by Hummel et al. It wraps the \code{\link[GlobalAncova]{GlobalAncova}} function in the \pkg{GlobalAncova} package.
}
\item{\code{global.test}: }
{
	This function performs a global gene set enrichment analysis using a global test by Goeman et al. It wraps the \code{\link[globaltest]{gt}} function in the \pkg{globaltest} package.
}

}
}

\references{
Hummel, M., Meister, R., Mansmann, U. (2008) GlobalANCOVA: exploration and assessment of gene group effects. \emph{Bioinformatics}, \bold{24}(1), 78--85.

Goeman, J. J., van de Geer, S. A., de Kort, F., van Houwelingen, H. C. (2004) A global test for groups of genes: testing association with a clinical outcome. \emph{Bioinformatics}, \bold{20}(1), 93--99.
}

\value{
A list containing the following items:
\item{pValue}{
	The p-value for the significance of \code{geneSet}.
}
\item{intersectGeneSetCoreSet}{
	This element is only returned in case of an overrepresentation analysis and consists of a vector of genes included in both sets (\code{geneSet} and \code{coreSet}).
}
\item{res.all}{
	The full result object returned by \code{\link[stats]{fisher.test}}, \code{\link[GlobalAncova]{GlobalAncova}} or \code{\link[globaltest]{gt}} respectively.
}
}

\seealso{
\code{\link{geneSetAnalysis}}, \code{\link{gsAnalysis}}
}
%\keyword{}
