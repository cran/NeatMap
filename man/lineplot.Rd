\name{lineplot}
\Rdversion{1.1}
\alias{lineplot}
\title{Line plots of profiles in a grid of embedded results}
\description{A 2d embedding of rows of given matrix is gridded and line plots of the profiles of the points within each cell are displayed}
\usage{
lineplot(pos, profiles, n.div.x = 10, n.div.y = 10, normalize = F, 
ylim=NULL, clipped=F)
}
\arguments{
  \item{pos}{2d positions of rows}
  \item{profiles}{Matrix of data to be plotted}
  \item{n.div.x}{Number of grid boxes in x direction}
  \item{n.div.y}{Number of grid boxes in y direction}
  \item{normalize}{logical: if true, the profiles are normalized before plotting}
  \item{ylim}{a vector of length 2 containing the \code{profile} values that correspond to the top and bottom of the grid boxes, if NULL the largest and smallest values in \code{profile} are used}
  \item{clipped}{logical: If TRUE, then values exceeding \code{ylim} are clipped}
  
}
\details{
        \code{pos} is assumed to be the embedding/2 dimensional representation of the rows of \code{profiles}. The \code{pos} result is then placed in a uniform grid with the number of divisions in the x and y directions specified by \code{n.div.x} and \code{n.div.y} respectively. Then in each grid cell the profiles for all the points in it are displayed together. Missing data is not plotted.
The grid extends 5 percent of the difference between the max and min point beyond these points. In each cell 90 percent of its width is used. By default, the profiles are scaled so that the maximum and minimum values (in \code{profiles}) would appear at the top and bottom of a cell. Different limits can be chosen using the \code{ylim} option. If \code{clipped} is true  values going beyond these limits will be clipped, in order to prevent the overlap of profiles in different cells.
If \code{normalize} is true, the profiles shall be normalized to have zero mean and unit variance. 
}
\value{
Returns a ggplot2 plot of class \code{ggplot}.
}
\author{Satwik Rajaram and Yoshi Oono}
\note{
\code{ylim} is applied to the profiles that will be plotted. So if normalization is turned on, the limits apply to these normalized values which are not accessible to the user. Therefore if this functionality is desired it may be better to normalize the data before invoking \command{lineplot}.
}
\examples{
#PCA and line plot. Notice how the profiles are dominated by two high value columns
mtcars.PCA<-prcomp(mtcars);
lineplot(mtcars.PCA$x,mtcars);

#Use ylim and clipping to allow us to focus on the columns with lower values
lineplot(mtcars.PCA$x,mtcars,ylim=c(0,10),clipped=TRUE)
}
\keyword{hplot}
