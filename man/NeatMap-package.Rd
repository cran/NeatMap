\name{NeatMap-package}
\Rdversion{1.1}
\alias{NeatMap-package}
\alias{NeatMap}
\docType{package}
\title{Non-clustered Heatmaps}
\description{
The \pkg{NeatMap} package is a set of functions to create heatmap like plots in two and three dimensions, without the need for cluster analysis. Like the heatmap, the plots created by \pkg{NeatMap} display both a dimensionally reduced representation of the data as well as the data itself. They are intended to be used in conjunction with dimensional reduction techniques such as PCA (as opposed to cluster analysis for the standard clustered heatmap).
}
\details{
\tabular{ll}{
Package: \tab NeatMap\cr
Type: \tab Package\cr
Version: \tab 0.3.1\cr
Date: \tab 2009-09-18\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}
The traditional clustered heatmap makes use of cluster analysis to  re-order rows and columns such that similar elements are placed together. However, cluster analysis is a poor choice for ordering method since it does not provide a unique ordering. The cluster analysis results are meant to be read along the tree, not in terms of the leaf order. The leaf order may be changed, while preserving the tree structure, by swinging the leaves at the bifurcations.  In fact there are methods that exploit this degree of freedom to improve the appearance of the heatmap. However, these are not standardized and could potentially place unrelated elements together. Apart from this, for the results of cluster analysis to be reliable the separation of groups need to be very pronounced. When this is not the case other dimensional reduction methods are likely to provide a better representation of the data. The \pkg{NeatMap} package is intended to be used in conjunction with such methods to display both the dimensional reduction result as well as the data underlying it. It includes plots in 2 and 3 dimensions. The two dimensional plots are built using the \pkg{ggplot2} package while the three dimensional plots use \pkg{rgl}.

The most basic functions are \command{heatmap1} and its convenience wrapper function \command{make.heatmap1}. This is virtually identical to the traditional heatmap, except the ordering of rows and columns do not use cluster analysis. The user may either supply an ordering of the rows using the method of their choice, or PCA/nMDS may be used. In the latter case, if normalized data or distance measure are used it is common to get an annular embedding in two dimensions (e.g. first two PCA components). The angular positions in this embedding are then used to order the rows and columns in the heatmap. However, the two opposite ends in this ordering will be separated by 360 degrees and are therefore very symmetric. To avoid artifacts produced by this, one may use the \command{make.circularmap} function which twists the heatmap into an annular format to reflect this periodicity (see examples below).

The \command{lineplot} takes the 2 dimensional embedding result, places it in a grid. Then for each grid cell, the profiles of all the points in that cell are displayed together as line plots. Lineplots are easier to comprehend than intensity patterns, so this format provides a good representation of the data. However, rows and columns are not treated on an equal footing, and comparison of genes is more difficult than the heatmap.

A circle like embedding was required to produce the ordering used in \command{make.heatmap1} and \command{make.circularmap}. In the more general case, where the embedding is not circular one may use the 3d plot \command{profileplot3d} (or its convenience wrapper function \command{make.profileplot3d}. The 2d embedding of rows is placed in the xy plane. For each point (i.e. row) in the xy plane, its profile, heatmap style, is shown parallel to the z axis.  It is possible to rotate and zoom this plot to focus on interesting parts. Stereo versions of these plots may be created using \command{stereo.profileplot3d} (or its convenience wrapper function \command{stereo.profileplot3d}). The stereo plots too are dynamically rotatable, and give the impression of observing truly three dimensional structure. Stereo plots could be useful in showing a 3D structure in publications (since rotations will not be possible).

There are also the dendrogram creating functions \command{draw.dendrogram} and its 3d version \command{draw.dendrogram3d} which are used internally by the functions described above. However, they may be called directly to compare the cluster analysis result to that of the dimensional reduction method.

For convenience purposes an implementation of non-Metric Multidimensional Scaling is also provided through the \command{nMDS} function.
}
\author{
Satwik Rajaram and Yoshi Oono
Maintainer: Satwik Rajaram <srajaram@illinois.edu>
}
\keyword{ package }
\seealso{
\code{\link{heatmap1}}, \code{\link{circularmap}}, \code{\link{lineplot}}, \code{\link{profileplot3d}},\code{\link{stereo.profileplot3d}}, \code{\link{draw.dendrogram}}, \code{\link{draw.dendrogram3d}}.
}
\examples{

#heatmap1 using pca (the scale_x_continuous ensures that labels can be seen)
make.heatmap1(mtcars,row.method="PCA",column.method="average.linkage",
row.labels=rownames(mtcars),column.labels=colnames(mtcars))+
scale_x_continuous(lim=c(-1,15))
#circularmap using nMDS and pearson correlation
make.circularmap(as.matrix(mtcars),metric="euclidean",cluster.method="complete.linkage",
normalize.profiles=FALSE,label.names=rownames(mtcars),label.size=3)

#lineplot using pca
mtcars.PCA<-prcomp(mtcars);
lineplot(mtcars.PCA$x,mtcars);

#profileplot3d  and stereo.profileplot3d using PCA
make.profileplot3d(mtcars,row.method="PCA",column.method="average.linkage")
make.stereo.profileplot3d(mtcars,row.method="PCA",column.method="average.linkage")

#dendrogram3d to compare nMDS result to clustering
mtcars.nMDS<-nMDS(mtcars,metric="euclidean")
mtcars.cluster<-hclust(dist(as.matrix(mtcars)),method="complete")
draw.dendrogram3d(mtcars.cluster,mtcars.nMDS$x,labels=rownames(mtcars),
label.size=0.5)
}
