\name{nMDS}
\Rdversion{1.1}
\alias{nMDS}
\title{ non-Metric Multi-Dimensional Scaling}
\description{ Given a matrix, and a distance measure, an embedding of the rows into desired Euclidean space is performed using non-Metric Multi-Dimensional Scaling.
}
\usage{
nMDS(data, embed.dim = 2, n.iters = 300, metric = "pearson",random.seed=NULL)
}
\arguments{
  \item{data}{matrix whose rows shall be embedded.}
  \item{embed.dim}{Dimensionality of Euclidean space into which embedding shall be performed.}
  \item{n.iters}{Number of iterations of the nMDS scheme}
  \item{metric}{The distance metric used to compare rows. Currently only "pearson" and "euclidean" ae supported.
}
  \item{random.seed}{A random seed used by nMDS. Use of this option allows reproducability of nMDS results}
}
\details{non-Metric Multi-Dimensional Scaling is performed using the scheme proposed by Taguchi and Oono.

If an element is missing (NA) in a particular row, all distance comparisons to that row shall ignore that particular element.
}
\value{An object of class "nMDS" containing:
\item{x}{matrix with the same number of rows and row names as \code{data} and having \code{embed.data} columns}}
\references{
Relational patterns of gene expression via non-metric multidimensional scaling analysis: Y.-h. Taguchi  and Y. Oono, Bioinformatics, 2005 21(6):730-740.
}
\author{Satwik Rajaram and Yoshi Oono}

\seealso{
\code{\link[pkg=stats]{prcomp}}
}
\examples{
#Two dimensional embedding
mtcars.nMDS<-nMDS(as.matrix(mtcars),embed.dim=2,metric="euclidean")
plot(mtcars.nMDS$x,type='n')
text(mtcars.nMDS$x,labels=rownames(mtcars.nMDS$x))
}
\keyword{multivariate}
\keyword{nonparametric}
