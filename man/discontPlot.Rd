\name{discontPlot}
\alias{discontPlot}

\title{

Plot with discontinuous smoother
}
\description{

Plots chromosome plots with a discontinuous smoother
}
\usage{
discontPlot(datalist,samples=c(1,2),chromNum=1,gamma=300,idiogram=FALSE,
file="default")
}

\arguments{
\item{datalist}{
The CAFE datalist to be analyzed, i.e. the output of
\code{\link{ProcessCels}}.
}
\item{samples}{
A vector or sample numbers to be plotted
}
\item{chromNum}{
the chromosome to be plotted
}
\item{gamma}{
The \code{gamma} level can be roughly compared to the sliding window size
in a normal continuous smoother. 
The \code{gamma} level determines how strict the algorithm functions; a
higher level will correspond to fewer jumps. This can not be higher than
the total number of probesets on the to-be-analyzed chromosome. Must be a
positive integer. 
}
\item{idiogram}{
if \code{TRUE}, will overlay a chromosome idiogram over the chromosome plot
}
\item{file}{
Specify a file name to store output png file
}
}

\value{
Plot to file system;
Returns a ggplot2 graph if \code{chromNum!="ALL"}. When \code{chromNum=="ALL"},
returns a list of ggplot2 graphs.

}
\references{
Friedrich, F., Kempe, a, Liebscher, V., & Winkler, G. (2008). 
Complexity Penalized M-Estimation. Journal of Computational and Graphical
Statistics, 17(1), 201-224. doi:10.1198/106186008X285591

}
\author{

Sander Bollen
}




\seealso{
\code{\link{rawPlot}}
\code{\link{slidPlot}}
\code{\link{facetPlot}}

}
\examples{
data("CAFE_data")
discontPlot(CAFE_data,samples=9,chromNum=17,gamma=300)
}

\keyword{dplot}
\keyword{hplot}
