\name{profileplot3d}
\Rdversion{1.1}
\alias{profileplot3d}
\alias{make.profileplot3d}
\title{Make a 3D rotatable plot showing data profiles}
\description{Make a 3d rotatable plot depicting the intensity levels of a matrix, while showing the relations between rows in two dimensions and that of the columns in the third. Cluster analysis results for rows and columns may be superposed.}
\usage{
profileplot3d(pos, profiles, normalize.rows = T, column.order = NULL, 
row.cluster = NULL, column.cluster = NULL, labels = NULL, col = NULL,
color_scaling_function=NULL, point.size = 3, label.colors = NULL, 
label.size = 0.5)

make.profileplot3d(profiles, row.method = "nMDS", normalize.rows=T, 
column.method = "average.linkage", row.metric = "pearson",
column.metric = "pearson", row.cluster.method = "average", 
column.cluster.method = "average", point.size = 3, col=NULL, 
color_scaling_function=NULL, labels = NULL, label.colors = NULL, 
label.size = 0.5, row.random.seed=NULL,column.random.seed=NULL)
}
\arguments{
  \item{pos}{matrix: the 2d positions for the rows in profiles as produced by any dimensional reduction scheme.}
  \item{profiles}{matrix: containing the data to be plotted.}
  \item{normalize.rows}{logical: If TRUE, then the rows shall be normalized before plotting.}
  \item{column.order}{The ordering of the columns, as would be the case in a typical heatmap, produced using some dimensional reduction scheme. If it is NULL, then the ordering in profiles is used.}
  \item{row.cluster}{hierarchical clustering result (of type \code{hclust}), of the rows, for superimposing the clustering result on the 3d profile plot. If it is NULL, no cluster result will be plotted.} 
  \item{column.cluster}{similar to \code{row.cluster} except for the clustering of columns.}
  \item{labels}{labels for the rows. If it is set to NULL, no labels will be plotted.}
  \item{col}{A list of colors such as that generated by \code{rainbow} used in depicting low to high intensities as in a heat-map.}
  \item{color_scaling_function}{A function mapping [0:1] onto [0:1] used for scaling the color levels. If Null, linear scaling is performed}
  \item{point.size}{The size of intensity points.}
  \item{label.colors}{A list of colors used for the row labels}
  \item{label.size}{Initial size of row labels. The sizes can be changed dynamically after plotting using the 3rd mouse button}
  \item{row.method}{dimensional reduction method for embedding rows, currently only "PCA" and "nMDS" are supported.}
  \item{column.method}{dimensional reduction method for ordering columns. Can be any one of "nMDS", "PCA", "average.linkage" and "complete.linkage".}
  \item{row.metric}{the distance function used for row embedding. Can be either "pearson" or "euclidean".}
  \item{column.metric}{like \code{row.metric} except for columns.}
  \item{row.cluster.method}{clustering method used for superposed row cluster. Either "complete.linkage" or "average.linkage".}
  \item{column.cluster.method}{like \code{row.metric} except for columns.}
  \item{row.random.seed}{Random seed to be used if nMDS is used to construct row structure}
  \item{column.random.seed}{Random seed to be used in nMDS is used to generate column ordering}
}
\details{These functions display data in a 3d rotatable format. The xy positions are the result of a 2D embedding of the rows, the profiles for which are shown along the z-axis. \command{profileplot3d} is the primary function to do this, accepting the results of any dimensional reduction scheme, while \command{make.profileplot3d} is a convenience function performing both the dimensional reduction (using nMDS or PCA) and then calling \command{profileplot3d}.

The ordering of columns is similar to that for \command{heatmap1}. Thus if PCA or nMDS are used it is assumed that the embedding of columns is annular, and the order used is that of angular positions. The user should confirm that this is indeed the case.
}
\value{These functions are called for the side-effects they produce}
\author{Satwik Rajaram and Yoshi Oono}
\seealso{
\code{\link{heatmap1}}.
}
\examples{
make.profileplot3d(mtcars,row.method="PCA",column.method="average.linkage")

#is equivalent to
mtcars.PCA<-prcomp(mtcars)
mtcars.col.cluster<-hclust(dist(t(mtcars)),method="average")
mtcars.row.cluster<-hclust(as.dist(1-cor(t(mtcars))),method="average")
profileplot3d(mtcars.PCA$x,mtcars,column.order=mtcars.col.cluster$order,
row.cluster=mtcars.row.cluster,column.cluster=mtcars.col.cluster)

#use of alternate colors and color scaling
make.profileplot3d(mtcars,row.method="PCA",column.method="average.linkage",
col=c("yellow","black","blue"),
color_scaling_function=function(x){0.5+tanh(10*(x-0.5))/2})


}
\keyword{hplot}
