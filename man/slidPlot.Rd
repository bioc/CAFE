\name{slidPlot}
\alias{slidPlot}

\title{
Plot with sliding average smoother

}
\description{ 

Plots chromosome plots with a moving average smoother
}
\usage{
slidPlot(datalist,samples=c(1,2),chromNum=1,combine=FALSE,k=1,idiogram=FALSE,file="default")
}

\arguments{
\item{datalist}{
The CAFE datalist to be analyzed, i.e. the output of \code{\link{ProcessCels}}.
}
\item{samples}{
A vector of sample numbers to be plotted
}
\item{chromNum}{
The chromosome to be analyzed
}
\item{combine}{
If \code{TRUE}, will plot the unaltered raw data in the background
}
\item{k}{
The sliding window size. Must be a positive integer, smaller than the total
number of probesets on the chromosome
}
\item{idiogram}{
If \code{TRUE}, will plot a chromosome idiogram over the plot
}
\item{file}{
Specify a file name to store output png fileS
}
}

\value{

Plot to file system;
Returns a ggplot2 graph if \code{chromNum!="ALL"}. When \code{chromNum=="ALL"},
returns a list of ggplot2 graphs.
}
\references{
H. Wickham. ggplot2: elegant graphics for data analysis. Springer New York, 2009.

}
\author{
Sander Bollen

}
\note{
Makes heavy use of the ggplot2 package.

}



\seealso{
\code{\link{rawPlot}}
\code{\link{facetPlot}}
\code{\link{discontPlot}}

}
\examples{
data("CAFE_data")
slidPlot(CAFE_data,samples=9,chromNum=17,k=50,combine=TRUE)
}

\keyword{dplot}
\keyword{hplot}
