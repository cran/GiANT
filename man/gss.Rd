\name{GeneSetStatistics}
\alias{gss}
\alias{gss.mean}
\alias{gss.sum}
\alias{gss.wilcoxonRankTest}
\alias{gss.maxmean}
\alias{gss.median}
\alias{gss.enrichmentScore}
\alias{gss.fisherExactTest}
\alias{gss.gsz}
\title{
Gene set statistics
}
\description{
Functions to calculate a gene set statistic, as used in the \code{gss} parameter of \code{gsAnalysis}. A gene set statistic summarizes a single gene set.
}
\usage{
gss.mean(x, geneSetIndices)

gss.sum(x, geneSetIndices)

gss.wilcoxonRankTest(x, geneSetIndices)

gss.maxmean(x, geneSetIndices)

gss.median(x, geneSetIndices)

gss.enrichmentScore(x, geneSetIndices, p = 1)

gss.fisherExactTest(x, geneSetIndices)

gss.gsz(x, geneSetIndices, w1 = 0.2, w2 = 0.5, preVar = 0, varConstant = 10)
}

\arguments{
\item{x}{
	A vector comprising one numeric value for each gene in the data set. This vector is usually obtained from the previous step, the gene-level statistic (see \code{\link{gls}}) or the transformed gene-level statistic (see \code{\link{transformation}}).
}
\item{geneSetIndices}{
	A vector containing the indices of the genes in the gene set with respect to the full gene set (i.e., the indices of the rows containing the measurements for these genes in \code{dat}).
}
\item{p}{
  Factor for \code{gss.enrichmentScore} that specifies the way hits are weighted. For \code{p = 0}, the enrichment score is a Kolmogorov-Smirnov statistic. For \code{p = 1} (the default), hits are weighted by their correlation.
}

\item{w1}{
  Weight for the median of the variance estimates for a gene set of size \code{varConstant}. Should be between 0 and 1. Default is \code{w1 = 0.2}.
}
\item{w2}{
  Weight for the median of the variance estimates for a gene set across the whole gene list. Should be between 0 and 1. Default is \code{w2 = 0.5}.
}
\item{preVar}{
  Parameter for incorporating the uncertainty of the observations. This is omitted by default (\code{preVar = 0}).
}
\item{varConstant}{
  Reference gene set size used for variance estimates. Default is \code{varConstant = 10}.
}
}
\details{
Standard functions for the calculation of gene set statistics (to be used in an analysis pipeline defined by \code{gsAnalysis}):
\itemize{
  \item{\code{gss.mean}: }
  {
  Calculates the mean of the (transformed) gene-level statistic values for the genes in the set.
  }

  \item{\code{gss.sum}: }
  {
  Calculates the sum of the (transformed) gene-level statistic values for the genes in the set.
  }
  
  \item{\code{gss.wilcoxonRankTest}: }
  {
  Calculates a Wilcoxon test comparing the (transformed) gene-level statistic values for the genes in the set versus those of the genes not in the set.
  }
  
  \item{\code{gss.maxmean}: }
  {
  Calculates the maximum of the means of positive and negative statistic values, weighted by the overall proportion of positive/negative values (e.g. for correlation scores where the sign denotes the direction). Described in Efron and Tibshirani.
  }
  
  \item{\code{gss.median}: }
  {
  Calculates the median of the (transformed) gene-level statistic values for the genes in the set.
  }
  
  \item{\code{gss.enrichmentScore}: }
  {
  Calculates the enrichment score of the (transformed) gene-level statistic values for the genes in the set, as described in Subramanian et al.
  }

  \item{\code{gss.fisherExactTest}: }
  {
  Performs Fisher's exact test to check gene sets for overrepresentation in the differential genes. This should be used used in combination with the transformation \if{latex}{\cr}\code{\link[=transformation]{transformation.adjustAndBinarize}}.
  }

  \item{\code{gss.gsz}: }
  {
  Calculates the Gene Set Z-score of the (transformed) gene-level statistic values for the genes in the set. Described in Toronen et al.
  }

}
}
\value{
Each method returns a single numeric value, the gene set statistic for the supplied gene set.
}
\references{
Efron, B., Tibshirani, R. (2007) On testing the significance of sets of genes. \emph{Annals of Applied Statistics}, \bold{1}, 107-129.

Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S., Mesirov, J. P. (2005) Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. \emph{Proceedings of the National Academy of Science of the United States of America}, \bold{102}, 15545--15550.

Toronen, P., Ojala, P. J., Marttinen, P., Holm L. (2009) Robust extraction of functional signals from gene set analysis using a generalized threshold free scoring function. \emph{BMC Bioinformatics}, \bold{10}(1), 307.
}

\seealso{
\code{\link{geneSetAnalysis}}, \code{\link{gsAnalysis}}, \code{\link{gls}}, \code{\link{transformation}}
}

%\keyword{}
