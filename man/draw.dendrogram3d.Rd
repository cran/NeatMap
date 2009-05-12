\name{draw.dendrogram3d}
\Rdversion{1.1}
\alias{draw.dendrogram3d}
\title{Draws a dendrogram in 3d}
\description{
Draw a dendrogram in 3d given clustering results and leaf positions. Labels at leaf tips may be specified}
\usage{
draw.dendrogram3d(cluster, positions, direction = c(0, 0, -1), scale = NULL,
heights = NULL, labels = NULL, label.colors = NULL, label.size = 3)
}
\arguments{
  \item{cluster}{Hierarchical clustering result of type \code{hclust} to be plotted.}
  \item{positions}{Positions (in 3d) of leaf tips. If 2d positions are given, these are placed in the z=0 plane.}
  \item{direction}{A vector specifying the direction in which the leaves point}
  \item{scale}{Maximum hieght of dendrogram}
  \item{heights}{A vector of heights, which may be used to over-ride the height information included in \code{cluster} (which are used by default if this is NULL)}
  \item{labels}{Vector of text labels to be placed at leaf tips.}
  \item{label.colors}{Vector of colors used used for labels.}
  \item{label.size}{Text size for labels.}
}
\details{
May be used to validate 2D embedding results with those of cluster analysis. When cluster results and the corresponding 2D embedding are specified, this function places the positions in a rotatable environment with the cluster analysis result superposed on it. This allows the user to understand the relationship between the clustering and embedded results. Labels (or desired colors) may be placed at the leaf tips.
}
\value{This function is called for the side effect it produces. It returns the id number of the last object drawn.
}
\author{Satwik Rajaram and Yoshi Oono}
\seealso{
\code{\link{draw.dendrogram}}
}
\examples{
mtcars.nMDS<-nMDS(mtcars,metric="euclidean")
mtcars.cluster<-hclust(dist(mtcars),method="complete")
draw.dendrogram3d(mtcars.cluster,mtcars.nMDS$x,labels=rownames(mtcars),label.size=0.5)
}
\keyword{hplot}
