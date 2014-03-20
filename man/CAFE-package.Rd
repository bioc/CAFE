\name{CAFE-package}
\alias{CAFE-package}
\alias{CAFE}
\docType{package}
\title{
Chromosomal Aberrations Finder in Expression data
}
\description{
CAFE attempts to find chromosomal aberrations in microarray expression (mRNA)
data. It contains several plotting functions to aid in visualizing these
aberrations. It generally recapitulates the workflow described by
Mayshar et al (see references), and implements several algorithms 
described by Friedrich et al (see references).
}
\details{
\tabular{ll}{
Package: \tab CAFE\cr
Type: \tab Package\cr
Version: \tab 0.6.9.5 \cr
Date: \tab 2013-05-16\cr
License: \tab GPLv3\cr
}
}
\author{
Sander Bollen
}
\references{
Friedrich, F., Kempe, a, Liebscher, V., & Winkler, G. (2008). Complexity
Penalized M-Estimation. Journal of Computational and Graphical Statistics, 
17(1), 201-224. doi:10.1198/106186008X285591

Mayshar, Y., Ben-David, U., Lavon, N., Biancotti, J.-C., Yakir, B.,
Clark, A. T., Plath, K., et al. (2010). Identification and classification of
chromosomal aberrations in human induced pluripotent stem cells. 
Cell stem cell, 7(4), 521-31. doi:10.1016/j.stem.2010.07.017
}

\keyword{ package }

\examples{
\dontrun{
setwd("/some/path/to/cel/files")
data <- ProcessCels() 
# process cel files
samples <- c(1,2) 
# select samples 1 and 2 to compare against the rest
chromosomeStats(data,chromNum="ALL",samples=samples) 
# check for chromosomal gains
chromosomeStats(data,chromNum="ALL",samples=samples,alternative="less") 
# check for chromosomal losses
bandStats(data,chromNum=1,samples=samples) 
# check for band gains in chr1
bandStats(data,chromNum=1,samples=samples,alternative="less") 
# check for band losses in chr1
rawPlot(data,chromNum=1,samples=samples,idiogram=TRUE) 
# plot raw data with an ideogram
slidPlot(data,chromNum=1,samples=samples,idiogram=TRUE,combine=TRUE,k=100) 
# moving average plot with ideogram
discontPlot(data,chromNum=1,samples=samples,idiogram=TRUE) 
# discontinuous plot with ideogram

}
}
