\name{ProcessCels}
\alias{ProcessCels}

\title{

Processing CEL files
}
\description{

Normalizes and computes relative expressions for all CEL files in work directory
}
\usage{
ProcessCels(threshold.over=1.5,threshold.under=(2/3),remove_method=1,
local_file=NULL)
}
\arguments{
\item{threshold.over}{
Determines the threshold, as a multiple of median value, where probes are
considered overexpressed. Default is 1.5
}
\item{threshold.under}{
Determines the threshold, as a fraction of median value, where probes are
considered underexpressed. Default is 2/3
}
\item{remove_method}{
Determines which method is used to remove multiple probesets that are
annotated to map to the same gene. The default option, \code{1}, 
will keep 1 probeset with the following priority:
1): \code{nnn_at};
2): \code{nnn_a_at};
3): \code{nnn_s_at};
4): \code{nnn_x_at};
5): lowest \code{nnn} if multiple probes still exist

If \code{remove_method=2}, probesets will \emph{only} be removed if
several probesets of the same gene map to the exact same location. 
In the case that many probesets map to the same location, one probeset 
will be retained according to the priority of option 1 above.

If \code{remove_method=0}, no multiple probesets will be removed 
}
\item{local_file}{
Use a local - previously downloaded - UCSC file
(e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/
affyU133Plus2.txt.gz) instead of directly retrieving the file instead.
}
}

\details{

this function uses the RMA algorithm to normalize *.CEL files in work directory. 
It then computes relative expressions for every probe on every sample.
Locations for probesets are downloaded from UCSC, as the
standard BioConductor annotations do not map probeset location
(they only map the location to the corresponding gene).
Multiple probesets belonging to the same gene are removed as described above.
The function then determines which probes are overexpressed and
underexpressed relative to the median probeset values across all samples.
Finally,the relative expressions are log2-transformed.
}
\value{

list
\item{$whole}{named list, where each element is a data.frame corresponding to
a *.CEL file - containing columns:
1): \code{"ID"} (Affy ID number); 
2): \code{"Sym"} (gene Symbol);
3): \code{"Value"} (Expression values); 
4): \code{"LogRel"} (Relative expressions);
5): \code{"Loc"} (Chromosomal locations);
6): \code{"Chr"} (Chromosome number);
7): \code{"Band"} (Cytoband);
8): \code{"Arm"} (Chromosomal arm)}
\item{$over}{same as $whole, but contains only those probes
which are deemed overexpressed}
\item{$under}{
same as $whole, but contains only those probes which are deemd underexpressed
}
}

\author{

Sander Bollen
}



\examples{
\dontrun{
data <- ProcessCels()
}
}
\keyword{datagen}