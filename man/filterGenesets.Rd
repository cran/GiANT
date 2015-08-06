\name{Filter gene sets}
\alias{filterGeneSets}
\title{
Filtering of gene sets
}
\description{
Filters gene sets according to different criteria.
}
\usage{
filterGeneSets(
	geneSets,
	includedGenes = NULL,
	minIntersectSize = length(includedGenes),
	adjMatrix = NULL,
	steps = NULL)
}
\arguments{
\item{geneSets}{
	A list of gene sets, where each gene set is a vector of gene names corresponding to the row names of \code{dat}.
}
\item{includedGenes}{
  A vector of gene names whose occurence in each of the gene sets is checked. The further parameters how these genes are used to filter gene sets
}
\item{minIntersectSize}{
  If this parameter is not \code{NULL}, only gene sets with an intersection of at least \code{minIntersectSize} genes with respect to \code{includedGenes} (or \code{includedGenes} expanded by its interactions if \code{adjMatrix} and \code{steps} are supplied) are included in the result set. By default, this is the size of \code{includedGenes}, \if{latex}{\cr}requiring \code{includedGenes} to be a subset of each gene set.
}
\item{adjMatrix}{
  An optional adjacency matrix in which an entry is \code{1} if there is a direct interaction between the corresponding genes and \code{0} otherwise. If this is non-null, the set of genes in \code{includedGenes} is expanded by adding all genes whose distance in the adjacency graph is at most \code{steps}.
}
\item{steps}{
  The maximum distance of interacting genes to the genes in \code{includedGenes} according to \code{adjMatrix} to be added to the expanded gene list. E.g., \code{steps = 1} means that all genes which are direct interaction partners of the initial genes in \code{includedGenes} are added to \code{includedGenes}.
}
}
\value{
Returns a filtered list of gene sets with the same structure as \code{geneSets}.
}

\seealso{
\code{\link{geneSetAnalysis}}, \code{\link{preprocessGs}}
}
\examples{
geneSets <- list(
		gs1 = paste("gene",1:20,sep=""),
		gs2 = paste("gene",50:60,sep=""),
		gs3 = paste("gene",90:92,sep=""),
		gs4 = paste("gene",55:65,sep="")
	)

newGeneSets1 <- filterGeneSets(
	geneSets = geneSets,
	includedGenes = c("gene55","gene60"))

newGeneSets2 <- filterGeneSets(
	geneSets = geneSets,
	includedGenes = c("gene1","gene55","gene20","gene100"),
	minIntersectSize = 2)

examplePathway <- c("gene1","gene2","gene3","gene4")
pathwayAdjMatrix <- matrix(0,100,100)
rownames(pathwayAdjMatrix) <- paste("gene",1:100,sep="")
colnames(pathwayAdjMatrix) <- paste("gene",1:100,sep="")

# gene1 interacts with gene2 and gene3
# -> step 1 if gene1 is starting point
pathwayAdjMatrix[1,2:3] <- 1
pathwayAdjMatrix[2:3,1] <- 1

# gene3 interacts with gene4
# -> step 2 if gene1 is starting point
pathwayAdjMatrix[3,4] <- 1
pathwayAdjMatrix[4,3] <- 1

newGeneSets3 <- filterGeneSets(
	geneSets = geneSets,
	includedGenes = "gene1",
	adjMatrix = pathwayAdjMatrix,
	steps = 2)
}

%\keyword{}
