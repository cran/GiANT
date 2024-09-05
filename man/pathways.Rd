\name{pathways}
\alias{pathways}
\docType{data}
\title{
Example data for pathways to to gene set enrichtment analyses
}
\description{
A list comprising exemplary pathways. Each pathway is represented by a vectro of random genes from countdata which are associated to this pathway.
}
\usage{data(exampleData)}
\details{
The data consists of a variable \code{pathways} of class \code{list} with 9 exemplary pathways. Each element of the list corresponds to one pathway - represented by a vector of gene names which are involved in this pathway. 
}

\examples{
data(exampleData)

print(pathways)

}