#' @title GC content 3rd position of synonymous codons
#'
#' @description
#' \code{ComputeGC3syn} Computes the G+C content 3rd position of synonymous codons
#'
#' @details
#' Should compute the same CP as in Darwin
#'
#' @param td Dictionary of String and Entry # nope probably either string or Entry type
#' @param pos List of 3 positive integers for the nucleotide positions
#' @return o/n
#'
#' @author Roth, A.; Friberg, M.; Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \code{\link{readstats}}
#' @keywords CodonBias
#' @examples
#' ComputeGC3syn('ATGTGGTACTCCGACTACGGAGGATAA')
#' ComputeGC3syn(c2s(mylist(whatout=1)[[1]]))
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
  # remove stop codon
  if (member(tD[-3:-1], AToCodon('$'))) { # introduce here the cds check and check for the stop codon missing warning
    d <- td[1:-4]
  }
  else {
    d <- tD
  }
  o <- 1:4
  n <-0
  for (i in seq(from=1, to=length(d), by=3)) {
    c <- d[i:i+2]
    if (length(IntToCInt(CodonToInt(c)))>1) {
      n <- n+1
      oi <- BToInt(c[3])
      o[oi] <- o[oi]+1
    }
  }
  o <- o/n
  return(o[2]+o[3])
}
