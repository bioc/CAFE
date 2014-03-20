\name{discontSmooth}
\alias{discontSmooth}

\title{
A discontinuous smoother

}
\description{

Calculates discontinuous smoother
}
\usage{
discontSmooth(y,gamma)
}

\arguments{
\item{y}{
input vector
}
\item{gamma}{
The \code{gamma} level can be roughly compared to the sliding window
size in a normal continuous smoother. 
The \code{gamma} level determines how strict the algorithm functions;
a higher level will correspond to fewer jumps. This cannot be larger
than \code{length(y)}. Must be a positive integer.
}
}
\details{
Uses the potts filter algorithm described by Friedrich et al. 

}
\value{
Vector with same length as input \code{y}

}
\references{
Friedrich, F., Kempe, a, Liebscher, V., & Winkler, G. (2008). 
Complexity Penalized M-Estimation. Journal of Computational and Graphical
Statistics, 17(1), 201-224. doi:10.1198/106186008X285591

}
\author{

Sander Bollen
}
\examples{
#generate piecewise vector with gaussian noise
y <- 1:450
y[1:150] <- 2
y[151:300] <- 3
y[301:450] <- 1
y <- y + rnorm(450)

#calculate smoother
y_smooth <- discontSmooth(y,20)
}

\keyword{smooth}
