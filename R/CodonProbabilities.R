#' @title Codon Probabilities
#'
#' @description
#' \code{CodonProbabilities} , for Each Codon, Computes the Probability that it Occurs at Least Once in a Gene
#'
#' @details
#' Should compute the same CP as in Darwin
#'
#' @param entries function for importing Entries() by default
#' @return Named (codons) numerical vector with relative synonymous codon usage for the 64 codons
#'
#' @author Roth, A.; Friberg, M.; Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \code{\link{readstats}}
#' @keywords CodonBias
#' @examples
#' RSCU('ATGTGGTACTCCGACTACGGAGGATAA')
#' RSCU(c2s(mylist(whatout=1)[[1]]))
#'
#' @import seqinr
#'
#' @export
#' @section Original code in Darwin:
#' \subsection{Compute CAI, the Codon Adaptation Index (Sharp and Li 1987)}{\preformatted{
#' CodonProbabilities := proc()
#'  res := CreateArray(1..64);
#'  for e in Entries() do
#'    occurs := CreateArray(1..64);
#'    dna := SearchTag('DNA', e);
#'    for c to length(dna) by 3 do
#'      cint := CodonToCInt(dna[c..c+2]);
#'      occurs[cint] := 1;
#'    od;
#'    res := res + occurs;
#'  od;
#'  res / DB[TotEntries]
#' end: } }
CodonProbabilities <- function(entries=Entries()) {
  res <- 1:64
  for (e in entries) {
    occurs <- 1:64
    dna <- SearchTag('DNA', e)
    for (c in seq(from=1, to=nchar(dna), by=3)) {
      cint <- CodonToCInt(dna[c..c+2])
      occurs[cint] <- 1
    }
    res <- res + occurs
  }
  return(res / length(entries))
}