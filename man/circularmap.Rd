\name{circularmap}
\Rdversion{1.1}
\alias{circularmap}
\alias{make.circularmap}
\title{Heatmap arranged as an annular region}
\description{These functions display a matrix as an annular heatmap reflecting the periodicity of the row ordering}
\usage{
circularmap(pos, profiles, column.order=NULL, cluster.result = NULL, 
cluster.heights = NULL, Rin = 10, Rout = 30, thickness = 3, label.names = NULL,
Rlabel = 32, label.size = 1.5, normalize.profiles = T)

make.circularmap(profiles, method = "nMDS", column.method="none",
cluster.method = "average.linkage", metric = "pearson", column.metric="pearson",
Rin = 10, Rout = 30, thickness = 3, label.names = NULL, Rlabel = 32,
label.size = 1.5, normalize.profiles = T, row.random.seed=NULL,column.random.seed=NULL)

}
\arguments{
  \item{pos}{The positions for the rows as produced by some dimensional reduction technique. Can either be list of angles, or of 2d positions.
}
  \item{profiles}{A matrix containing the values to be displayed as a heatmap.}
  \item{column.order}{A vector containing the order in which columns should be shown}
  \item{cluster.result}{The hierarchical clustering result of type \code{hclust} for the rows, superposed on the heatmap for validation. If NULL no clustering result is shown.}
  \item{cluster.heights}{An optional vector of heights to over-ride the heights specified in \code{cluster.result}.}
  \item{Rin}{Inner radius of annulus.}
  \item{Rout}{Outer radius of annulus.}
  \item{thickness}{Thickness of the individual row expression level bands.}
  \item{label.names}{list of row labels.}
  \item{Rlabel}{Radius at which row labels are placed.}
  \item{label.size}{Font size for labels.}
  \item{normalize.profiles}{logical: if true the profiles are normalized (mean zero, unit variance) before display.}
  \item{method}{The dimensional reduction method used by \code{make.circularmap} to produce angular positions for the rows. One of "nMDS" or "PCA".}
  \item{column.method}{The dimensional reduction method used by \code{make.circularmap} to order columns. One of "none","nMDS","PCA", "average.linkage" or "complete.linkage".}
  \item{cluster.method}{clustering method used for superposed row cluster. Either "complete.linkage" or "average.linkage".}
  \item{metric}{the distance function used for row embedding. Can be either "pearson" or "euclidean".}
  \item{column.metric}{the distance function used for column embedding. Can be either "pearson" or "euclidean".}
  \item{row.random.seed}{Random seed to be used if nMDS is used to construct row embedding}
  \item{column.random.seed}{Random seed to be used in nMDS is used to generate column ordering}
}
\details{ These are function used to construct heatmap like plots arranged in an annular ring. The assumption is that the 2d embedding result for the rows of the given matrix is circular in shape (the user should confirm this is indeed the case). The angular positions of each point (i.e., row) in this embedding is the angular position of its profile in the circular map. \command{circularmap} is the base function that takes the expression matrix and the result of dimensional reduction (expressed as a 2d positions or a list of angles) and produces the annular heatmap. \command{make.circularmap} is a convenience wrapper function that takes the given matrix, performs dimensional reduction using either "nMDS" or "PCA", produces cluster analysis on the row, and passes these results to \command{circularmap} for plotting.

Apart from the circular (vs linear) ordering, this function differs from \command{heatmap1} in that the profiles are no longer equally spaced, and may reflect the non-uniformity of the angular positions. On the other hand, \command{heatmap1} treats rows and columns on an equal footing while these functions focus primarily on the relations between rows. 

The cluster analysis makes use of the same distance measure (specified by \code{metric}) as the dimensional reduction method.
}
\value{A \pkg{ggplot2} plot of class \code{ggplot}.}
\author{Satwik Rajaram and Yoshi Oono}
\seealso{\code{\link{image}},\code{\link{heatmap}},\code{\link{heatmap1}}.}
\examples{
make.circularmap(as.matrix(mtcars),metric="euclidean",cluster.method="complete.linkage",
normalize.profiles=FALSE)

#is equivalent to
mtcars.nMDS<-nMDS(as.matrix(mtcars),metric="euclidean")
mtcars.cluster<-hclust(dist(mtcars),method="complete")
circularmap(mtcars.nMDS$x,as.matrix(mtcars),normalize.profiles=FALSE,
cluster.result=mtcars.cluster)

#To change coloring scheme etc
make.circularmap(as.matrix(mtcars),metric="euclidean",cluster.method="complete.linkage",
normalize.profiles=FALSE)+scale_colour_gradient2(low="yellow",high="blue",
mid="black",midpoint=200);
}
\keyword{hplot}
