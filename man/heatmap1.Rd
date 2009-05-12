\name{heatmap1}
\Rdversion{1.1}
\alias{heatmap1}
\alias{make.heatmap1}
\title{Make a non-clustered heatmap}
\description{Makes a heatmap without need for cluster analysis}
\usage{
heatmap1(profiles, row.order = NULL, column.order = NULL, row.cluster = NULL,
column.cluster = NULL, column.labels = NULL, row.labels = NULL, 
column.label.size = 3, row.label.size = 3,row.normalize=F)

make.heatmap1(profiles, row.method = "nMDS", column.method = "none",
row.metric = "pearson", column.metric = "pearson", row.cluster.method = "average",
column.cluster.method = "average", column.labels = NULL, row.labels = NULL,
row.label.size = 3, column.label.size = 3,row.normalize=F)
}
\arguments{
  \item{profiles}{matrix: containing the data to be plotted.}
  \item{row.order}{vector containing order of rows such as produced by \code{order}}
  \item{column.order}{vector containing order of columns such as produced by \code{order}}
  \item{row.cluster}{hierarchical clustering result for rows of type \code{hclust}}
  \item{column.cluster}{hierarchical clustering result for columns of type \code{hclust}}
  \item{column.labels}{vector of labels for columns}
  \item{row.labels}{vector of labels for rows}
  \item{column.label.size}{size for column label text}
  \item{row.label.size}{size for row label text}
  \item{row.normalize}{logical: If true the rows are normalized to zero mean and unit variance}
  \item{row.method}{dimension reduction method used by \code{make.heatmap1} to order rows. One of "none","nMDS","PCA","complete.linkage" and "average.linkage".}
  \item{column.method}{dimension reduction method used by \code{make.heatmap1} to order columns. One of "none","nMDS","PCA","complete.linkage" and "average.linkage".}
  \item{row.metric}{Distance metric used by \code{row.method} either "pearson" or "euclidean"}
  \item{column.metric}{Distance metric used by \code{column.method} either "pearson" or "euclidean"}
  \item{row.cluster.method}{Clustering algorithm used for clustering rows. Either "average.linkage" or "complete.linkage". If NULL, no row cluster results are shown}
  \item{column.cluster.method}{Clustering algorithm used for clustering columns. Either "average.linkage" or "complete.linkage". If NULL, no column cluster results are shown}
}
\details{The traditional heatmap uses clustering to order rows and columns. These functions allow us to use alternate schemes for this ordering. They use the same format as the traditional heatmap, and are therefore similar to \command{image} and \command{heatmap}. \command{heatmap1} assumes the user has already ordered the rows and columns according to the scheme of their choice. \command{make.heatmap1} is a convenience wrapper which performs ordering using "nMDS", "PCA", or hierarchical clustering. This ordering is the passed on to \code{heatmap1} for plotting.
If "nMDS" or "PCA" are chosen as the ordering method, then it is assumed that their 2D embedding is annular in shape. This is often observed if \acronym{PCA} or nMDS (with euclidean distance) is applied to normalized data or Pearson correlation is used with nMDS. Angles measured at the centre of mass are then used for ordering. The user is therefore advised to confirm that such an annular structure is indeed present. Note that the two opposite ends of the ordering are typically separated by 360 degrees and are therefore very similar. To avoid artifacts produced by this, the user could consider using \command{circularmap} instead.
The metric specified by \code{row.metric} and \code{column.metric} are also used by the clustering algorithms.}
\value{
A \pkg{ggplot2} plot of class \code{ggplot}.
}
\author{Satwik Rajaram and Yoshi Oono}
\seealso{
\code{\link[graphics]{image}},\code{\link[stats]{heatmap}},\code{\link{circularmap}}.
}
\examples{
make.heatmap1(mtcars,row.method="PCA",column.method="average.linkage")

#is equivalent to
mtcars.PCA<-prcomp(mtcars)
mtcars.PCA.order<-order(apply(mtcars.PCA$x[,1:2],1,function(x){atan2(x[1],x[2])}))
mtcars.column.cluster<-hclust(as.dist(1-cor(mtcars)),method="average")
mtcars.row.cluster<-hclust(as.dist(1-cor(t(mtcars))),method="average")
heatmap1(mtcars,row.order=mtcars.PCA.order,column.order=mtcars.column.cluster$order,
row.cluster=mtcars.row.cluster,column.cluster=mtcars.column.cluster)

#Changing Color Scheme
make.heatmap1(mtcars,row.method="PCA",column.method="average.linkage")+
scale_fill_gradient2(low="yellow",high="blue",mid="black",midpoint=200)

#Adding labels (the scale function ensures that labels are not clipped)
make.heatmap1(mtcars,row.method="PCA",column.method="average.linkage",
row.labels=rownames(mtcars),column.labels=colnames(mtcars))+
scale_x_continuous(lim=c(-1,15))
}
\keyword{hplot}
