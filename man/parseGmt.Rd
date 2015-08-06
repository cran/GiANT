\name{parse GMT files}
\alias{parseGmt}
\title{
Venn Euler Diagramm
}
\description{
Parses a GMT file as downloadable from MSigDB (presented in Subramanian et al.) and returns a list of gene sets.
}
\usage{
	parseGmt(file)
}

\arguments{
\item{file}{
	A file name.
}
}
\details{
Parses a GMT file and returns a list of gene sets. Each list element named according to the included gene set. The gene set files can be downloaded from http://www.broadinstitute.org/gsea/msigdb.
}
\value{
	A named \code{list} of gene sets.
}

\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S., Mesirov, J. P. (2005) Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. \emph{Proceedings of the National Academy of Science of the United States of America}, \bold{102}, 15545--15550.
}

\seealso{
\code{\link{geneSetAnalysis}}, \code{\link{predefinedAnalyses}}, \code{\link{gsAnalysis}}
}
\examples{
\donttest{
  # download the file from MSigDB first and place it in the workspace directory
	geneSets <- parseGmt("msigdb.v4.0.entrez.gmt")
}
}

%\keyword{}
