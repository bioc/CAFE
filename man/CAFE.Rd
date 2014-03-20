\name{CAFE_data}
\alias{CAFE_data}
\docType{data}
\title{
CAFE data set
}
\description{
Contains the dataset of GSE6561 and GSE10809 processed by 
\code{\link{ProcessCels}} 
}
\usage{data("CAFE_data")}
\format{
A list containing two lists
\describe{
\item{\code{whole}}{A list containing a dataframe for each sample}
\item{\code{over}}{A list containing a dataframe for each sample,
but with only those probes that are deemed overexpressed}
}
The dataframes inside the lists contain the following columns:
 
\describe{
\item{\code{ID}}{Affymetrix probe IDs}
\item{\code{Sym}}{Gene symbols}
\item{\code{Value}}{Log2 transformed expression values}
\item{\code{LogRel}}{Log2 transformed relative expression values
(to the median)}
\item{\code{Loc}}{Chromosomal locations}
\item{\code{Chr}}{Chromosome identifiers}
}
}

\source{
GSE6561: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6561 

GSE10809: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10809
}

\examples{
data("CAFE_data")
}
\keyword{datasets}
