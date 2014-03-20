\name{slidSmooth}
\alias{slidSmooth}

\title{
A moving average smoother

}
\description{

Calculates moving average smoother
}
\usage{
slidSmooth(x,k)
}

\arguments{
\item{x}{
input vector
}
\item{k}{
The moving average window size. Must be an integer value greater than 0,
and no larger than \code{length(y)}.
}
}
\value{
Vector with same length as input \code{y}

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
y_smooth <- slidSmooth(y,20)
}

\keyword{smooth}
