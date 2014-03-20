\name{cliSubset}
\alias{cliSubset}

\title{
Subset data with a CLI

}
\description{

Provides command line interface for subsetting input datasets
}
\usage{
cliSubset(datalist,alternative)
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

subset of input
}
\author{

Sander Bollen
}



\seealso{

\code{\link{guiSubset}}
}
\examples{
\dontrun{
datalist <- data("CAFE_data")
sub <- cliSubset(datalist,alternative="greater")
}
}
\keyword{manip}
