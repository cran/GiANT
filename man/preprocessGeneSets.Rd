\name{preprocessGeneSets}
\alias{preprocessGeneSets}
\alias{preprocessGs}
\title{
Eliminate unknown genes from gene sets
}
\description{
This function removes all genes that are not part of the experiment (not in \code{rownames(dat)}) from the specified gene sets which. All names are set to lower case.
}
\usage{
preprocessGs(
	dat,
	geneSets)
}

\arguments{
\item{dat}{
	A numeric matrix of gene expression values for all analyzed genes. Here, each row corresponds to one gene, and each column corresponds to one sample. The rows must be named with the gene names used in the gene sets. Here, only the row names (i.e. the gene names) are used by the function.
}
\item{geneSets}{
	A list of gene sets to be processed, where each gene set is a vector of gene names corresponding to the row names of \code{dat}.
}
}

\value{
A list of preprocessed gene sets, where each gene set only contains those genes that are also present in \code{dat}
}

\seealso{
\code{\link{geneSetAnalysis}}
}
\examples{
#values are not important, only the row names are used
dat <- matrix(0,100,10)
rownames(dat) <- paste("gene",1:100,sep="")

geneSets <- list(
		gs1 = paste("GENE",1:20,sep=""),# all genes in the analyzed data
		gs2 = paste("Gene",101:110,sep=""),#no gene in the analyzed data
		gs3 = paste("gene",90:110,sep="")#some genes in the analyzed data
	)

newGeneSets <- preprocessGs(dat = dat, geneSets = geneSets)
}

%\keyword{}
