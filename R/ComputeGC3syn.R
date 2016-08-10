#' @title GC content 3rd position of synonymous codons
#'
#' @description
#' \code{ComputeGC3syn} Computes the G+C content 3rd position of synonymous codons
#'
#' @details
#' Should compute the same CP as in Darwin. By definition, GC3s values are the proportion of GC
#' nucleotides at the variable third coding position of synonymous codons.
#' This can be used to evaluate the degree of base composition bias.
#'
#' @param tD Nucleotide character string
#' @return o/n
#'
#' @author Roth, A.; Friberg, M.; Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \code{\link{readstats}}
#' @keywords CodonBias
#' @examples
#' ComputeGC3syn('ATGTGGTACTCCGACTACGGAGGATAA')
# list arguments have to be properly implemented
# ComputeGC3syn(c2s(mylist(whatout=1)[[1]]))
#'
#' @import seqinr
#'
#' @export
#' @section Original code in Darwin:
#' \subsection{G+C content 3rd position of synonymous codons.  AR (April 2007)}{\preformatted{
#' ComputeGC3syn:= proc(td:string)
#'  if member(td[-3..-1], AToCodon('$')) then d:=td[1..-4] else d:=td fi; # remove stop codon
#'  o := CreateArray(1..4);
#'  n:=0;
#'  for i to length(d) by 3 do
#'    c:=d[i..i+2];
#'    if length(IntToCInt(CodonToInt(c)))>1 then
#'      n:=n+1;
#'      oi:=BToInt(c[3]);
#'      o[oi] := o[oi]+1
#'    fi;
#'  od;
#'  o:=o/n;
#'  return(o[2]+o[3]);
#' end: } }
ComputeGC3syn <- function(tD) {
  require(seqinr)
  if(!(checkCDS(tD))) stop("non valid CDS)", call.=FALSE)
  # remove stop codon
  if (c2s(s2c(tD)[(nchar(tD)-2):nchar(tD)]) %in% substr(names(aa_ac), 5, 7)[c(59, 60, 64, 65)]) {
    d <- substring(tD, 1, nchar(tD)-2)
  }
  else {
    d <- tD
  }
  bases <- c("A", "C", "G", "T")
  o <- rep(0, times=4)
  names(o) <- bases
  n <-0
  for (i in seq(from=1, to=nchar(d), by=3)) {
    c <- substring(d, i, i+2)
    if (c %in% substr(names(aa_ac), 5, 7)[c(-59, -60, -64, -65)]) {
      n <- n+1
      oi <- which(names(o) %in% substring(c, 3))
      o[oi] <- o[oi]+1
    }
  }
  o <- o/n
  return(o[2]+o[3])
}
