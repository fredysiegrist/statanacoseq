% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CodonProbabilities.R
\name{CodonProbabilities}
\alias{CodonProbabilities}
\title{Codon Probabilities}
\usage{
CodonProbabilities(entries = mylist(whatout = 1))
}
\arguments{
\item{entries}{function for importing Entries() by default}
}
\value{
Named (codons) numerical vector with relative synonymous codon usage for the 64 codons
}
\description{
\code{CodonProbabilities} , for Each Codon, Computes the Probability that it Occurs at Least Once in a Gene
}
\details{
Should compute the same CP as in Darwin
}
\section{Original code in Darwin}{

\subsection{Compute CAI, the Codon Adaptation Index (Sharp and Li 1987)}{\preformatted{
CodonProbabilities := proc()
 res := CreateArray(1..64);
 for e in Entries() do
   occurs := CreateArray(1..64);
   dna := SearchTag('DNA', e);
   for c to length(dna) by 3 do
     cint := CodonToCInt(dna[c..c+2]);
     occurs[cint] := 1;
   od;
   res := res + occurs;
 od;
 res / DB[TotEntries]
end: } }
}
\examples{
CodonProbabilities('ATGTGGTACTCCGACTACGGAGGATAA')
CodonProbabilities(c2s(mylist(whatout=1)[[1]]))

}
\author{
Roth, A.; Friberg, M.; Siegrist, F. and Cannarozzi, G. M. \email{gina@cannarozzi.com}
}
\seealso{
\code{\link{seqinr}} \code{\link{statanacoseq}} \code{\link{readstats}}
}
\keyword{CodonBias}

