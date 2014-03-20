\name{facetPlot}
\alias{facetPlot}
\title{
Plot all chromosomes horizontally next to each other
}
\description{
Plots all chromosomes in horizontal alignment next to each other, 
with optionally a moving average smoother applied to the data 
}
\usage{ 
facetPlot(datalist,samples=c(1,2),slid=FALSE,combine=FALSE,k=1,file="default")
}
\arguments{
\item{datalist}{
The CAFE datalist to be analyzed, i.e. the output of \code{\link{ProcessCels}}.
}
\item{samples}{
A vector or sample numbers to be plotted
}
\item{slid}{
If \code{TRUE}, use moving average smoother
}
\item{combine}{
If \code{TRUE}, will plot the unaltered raw data in the background
}
\item{k}{
The sliding window size. Must be a positive integer, smaller than the 
length of Affy IDs on the chromosome
}
\item{file}{
Specify a file name to store output png file
}
}
\value{
Plot to file system.
Return a ggplot2 graph
}
\references{
H. Wickham. ggplot2: elegant graphics for data analysis. 
Springer New York, 2009.
}
\author{
Sander Bollen
}
\note{
Makes heavy use of the ggplot2 package
}


\seealso{
\code{\link{slidPlot}}
\code{\link{rawPlot}}
\code{\link{discontPlot}}
}
\examples{
data("CAFE_data")
facetPlot(CAFE_data,samples=9)
}
\keyword{dplot}
\keyword{hplot}

