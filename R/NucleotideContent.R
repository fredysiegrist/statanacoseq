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
#' @return Vector of four fractions (for every base) of nucleotide content in all bases
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
NucleotideContent <- function( tD=as.list("", Entries()), pos=c(1,2,3)) {
  bases <- c("A", "C", "G", "T")
  o <- rep(0, times=4)
  names(o) <- bases
  n <- 0
  if (missing(tD)) {
    for (z in 1:length(DB)) {
      d <- toupper(c2s(entries[[z]]))
      for (i1 in seq(from=1, to=abs(nchar(d)-max(pos)), by=3)) {
        for (i2 in pos) {
          i <- i1+i2
          n <- n+1
          o[s2c(d)[i]] <- o[s2c(d)[i]]+1
        }
      }
    }
  }
# if Entry type will be defined as in Darwin use this convertion
#  else if type(tD, Entry) then d:=SearchTag('DNA', tD)
  else d <- tD
  for (i1 in seq(from=1, to=abs(nchar(d)-max(pos)), by=3)) {
    for (i2 in pos) {
      i <- i1+i2
      n <- n+1
      o[s2c(d)[i]] <- o[s2c(d)[i]]+1
    }
  }
  return(o/n)
}
