#' @title Nucleotide Content
#'
#' @description
#' \code{NucleotideContent} Computes for a gene or a database entry the content of nucleotides
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
#' NucleotideContent('ATGTGGTACTCCGACTACGGAGGATAA')
#' NucleotideContent(c2s(mylist(whatout=1)[[1]]))
#'
#' @import seqinr
#'
#' @export
#' @section Original code in Darwin:
#' \subsection{Compute Nucleotide content.  AR (2006)}{\preformatted{
#' NucleotideContent := proc( ; tD:{string, Entry}, pos=[1,2,3]:list(posint)) -> list;
#'  o := CreateArray(1..4);
#'  n:=0;
#'  if not assigned(tD) then
#'    for z to DB[TotEntries] do
#'      d:=SearchTag(DNA, Entry(z));
#'      for i1 to length(d)-max(pos) by 3 do for i2 in pos do
#'        i:=i1+i2;
#'        n:=n+1;
#'        o[BToInt(d[i])] := o[BToInt(d[i])]+1
#'      od od;
#'    od;
#'  else
#'    if type(tD, Entry) then d:=SearchTag('DNA', tD)
#'      else d:=tD fi;
#'    for i1 to length(d)-max(pos) by 3 do for i2 in pos do
#'      i:=i1+i2;
#'      n:=n+1;
#'      o[BToInt(d[i])] := o[BToInt(d[i])]+1
#'    od od;
#'  fi;
#'  return(o/n);
#' end: } }
NucleotideContent <- function( tD=as.list("", Entry()), pos=list(c(1,2,3))) {
  o <- 1:4
  n <- 0
  if (missing(tD)) {
  #  for z to DB[TotEntries] do
  #  d:=SearchTag(DNA, Entry(z));
  #  for i1 to length(d)-max(pos) by 3 do for i2 in pos do
  #   i:=i1+i2;
  #   n:=n+1;
  #   o[BToInt(d[i])] := o[BToInt(d[i])]+1
  #   od od;
  #   od;
  }
#  else if type(tD, Entry) then d:=SearchTag('DNA', tD)
  else d <- tD
  for (i1 in seq(from=1, to=length(d)-max(pos), by=3)) {
    for (i2 in pos) {
      i <- i1+i2
      n <- n+1
      o[BToInt(d[i])] <- o[BToInt(d[i])]+1
    }
  }
  return(o/n)
}
