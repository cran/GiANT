\name{mergeProbesForGenes}
\alias{mergeProbesForGenes}
\title{
Merge multiple probes for one gene
}
\description{
Merges all probes belonging to the same gene by identifying duplicate row names in a data matrix.
}
\usage{
mergeProbesForGenes(dat,
	method = c("mean", "max", "min", "median"))
}
\arguments{
\item{dat}{
	A numeric matrix of gene expression values for all analyzed genes. Here, each row corresponds to one probe, and each column corresponds to one sample. The rows must be named with the gene names and may contain duplicates if multiple probes correspond to the same gene.
}
\item{method}{  
	The method which should be used to merge probes entries in \code{dat}. Depending on the chosen method, the merged value for a gene and a specific sample is defined as the mean value, the maximum value, the minimum value or the median of all probes of this sample belonging to the gene.
}
}
\value{
A matrix of the same structure as \code{dat}, but possibly with fewer rows if probes were merged.
}

\seealso{
\code{\link{geneSetAnalysis}}
}
\examples{
dat <- matrix(1:6, nrow=3, ncol=2)
rownames(dat) <- c("g1", "g2", "g1")

newDat <- mergeProbesForGenes(dat, method = "mean")
}

%\keyword{}
