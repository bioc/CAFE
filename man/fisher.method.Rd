\name{fisher.method}
\alias{fisher.method}

\title{
Combines pvalues by using Fisher's method
}
\description{

Combines pvalues by using Fisher's method
}
\usage{
fisher.method(pvals)
}

\arguments{
\item{pvals}{
Vector of p values
}
}
\value{
Combined p value
}
\author{
Sander Bollen
}

\examples{
pvals <- runif(20) #generate 20 pvals
fisher.method(pvals)
}

\keyword{htest}
