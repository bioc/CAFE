\name{rawPlot}
\alias{rawPlot}

\title{
Plot without any smoother

}
\description{

Makes chromosome plot using raw data values
}
\usage{
rawPlot(datalist,samples=c(1,2),chromNum=1,idiogram=FALSE,file="default")
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
The chromosome to be analyzed
}
\item{idiogram}{
If \code{TRUE}, will plot a chromosome idiogram over the plot
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

\author{

Sander Bollen
}



\seealso{

\code{\link{slidPlot}}
\code{\link{facetPlot}}
\code{\link{discontPlot}}
}
\examples{
data("CAFE_data")
rawPlot(CAFE_data,samples=8,chromNum=17)
}
\keyword{dplot}
\keyword{hplot}
