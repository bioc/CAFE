\name{guiSubset}
\alias{guiSubset}

\title{

Subset data with a GUI
}
\description{

Provides graphical user interface for subsetting input datasets
}
\usage{
guiSubset(datalist,alternative)
}

\arguments{
  \item{datalist}{

the dataset to be subsetted
}
\item{alternative}{
"greater" or "less"
}
}

\value{

Subset of input to variable \code{guiSelectedSet} in working directory
}

\author{

Sander Bollen
}




\seealso{

\code{\link{cliSubset}}
}
\examples{
\dontrun{
data("CAFE_data")
guiSubset(CAFE_data,alternative="greater")
}
}
\keyword{manip}
