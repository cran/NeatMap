\name{draw.dendrogram}
\Rdversion{1.1}
\alias{draw.dendrogram}
\title{Draws a dendrogram in 2d}
\description{
Draw a dendrogram in 2d given clustering results. Leaf order may be specified, allowing comparison to clustering result.}
\usage{
draw.dendrogram(cluster, leaf.order=NULL, scale = 10, dendro.dir = "left", 
order.dir = "up", origin = as.vector(c(0.5, 0)), heights = NULL)
}
\arguments{
  \item{cluster}{Hierarchical clustering result of type \code{hclust} to be plotted.}
  \item{leaf.order}{A vector containing the order of leaf tips, if NULL, the order specified in \code{cluster} is used}
  \item{scale}{Maximum dendrogram height. Width is equal to number of leaves}
  \item{dendro.dir}{Direction in which the leaves point. One of "up", "down", "left" or "right"}
  \item{order.dir}{Direction in which leaves are ordered. Should be perpendicular to \code{dendro.dir}. One of"up", "down", "left" or "right"} 
  \item{origin}{Position of first leaf}
  \item{heights}{A vector of heights, which may be used to over-ride the height information included in \code{cluster} (which are used by default if this is NULL)}
}
\details{Similar to the \command{plot} command for cluster results of type \command{hclust}. Primarily intended for internal use, but for ease of formatting, may be of use in producing standard dendrogram plots instead of the standard \command{hclust} commands. Also the leaf order, and branch heights used here may be different from that specified in \code{cluster}.
The leaf tips are placed at the integer valued positions (with respect to the origin) specified in the leaf order. This may be used to place labels, colored points etc at the leaf tips.
}
\value{a \code{ggplot} layer}
\author{Satwik Rajaram and Yoshi Oono}
\note{
The leaf heights are scaled so that the maximum height corresponds to \code{scale}.
}

\seealso{
\code{\link[pkg:stats]{hclust}}.
}
\examples{
#compare complete and average linkage
mtcars.cluster<-hclust(dist(mtcars),method="complete")
mtcars.cluster.avg<-hclust(dist(mtcars),method="average")
ggplot.default()+draw.dendrogram(mtcars.cluster,leaf.order=mtcars.cluster.avg$order)
}
\keyword{aplot}
