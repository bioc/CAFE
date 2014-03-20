\name{chromosomeStats}
\alias{chromosomeStats}

\title{
Find aberrations with whole-chromosome resolution

}
\description{
Calculate significant chromosomes with various statistical tests

}
\usage{
chromosomeStats(datalist, chromNum=1, samples=NULL, select="cli", test="fisher",
bonferroni = TRUE, enrichment = "greater")
}

\arguments{
\item{datalist}{
The CAFE datalist to be analyzed, i.e. the output of 
\code{\link{ProcessCels}}.
}
\item{chromNum}{
The chromosome to be calculated. This can be \code{"ALL"} to calculate all
chromosomes.
}
\item{samples}{
A vector containing sample numbers to be analyzed
}
\item{select}{
Signifies which type of sample selection prompt will be shown, if
\code{samples=NULL}. Currently supported are \code{"cli"} for a command line
interface and \code{"gui"} for a tcl/tk-based graphical user interface.
}
\item{test}{
Signifies which statistical test to be used in the final calculation. 
Must be either \code{"fisher"} for an exact fisher test or \code{"chisqr"}
for a chi square test. 
}
\item{bonferroni}{
If \code{bonferroni=TRUE}, will correct the p-values of the enrichment test
with a bonferroni method.
}
\item{enrichment}{
Test for over or underexpression. Can be set to \code{"greater"} or
\code{"less"}.
}
}

\value{
A named vector containing p-values.
}

\author{
Sander Bollen
}
\note{
Technically speaking, the Fisher's exact test is better than the 
chi-square test; the Fisher's exact test gives an exact p-value, whereas the
chi-square test only gives an approximation. However, the Fisher's exact test
can get slow for large sample sizes, and the chi-square test becomes better 
with increasing sample size but does not slow down as much.
}

\seealso{
\code{\link{bandStats}}
\code{\link{armStats}}
}
\examples{
data("CAFE_data")
sam <- c(9,11)
chromosomeStats(CAFE_data,chromNum=17,samples=sam,test="fisher")
}
\keyword{multivariate}
