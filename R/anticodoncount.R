#' @title Calculate Global Anticodon Usage
#'
#' @description
#' \code{anticodoncount} sums up the anticodon usage in a mRNA and a protein fasta file
#'
#'
#' @details
#' This function should count the codon frequency in a genome based on a set of verified and matched transcript and protein fasta files.
#' This is a roxygen2 created helping file.
#'
#' @param entry Index number in SeqFastadna list myaa
#' @param accnt Vector of 65 anticodons used
#'
#' @return Codon Frequency Table
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}}
#' @keywords CodonBias
#' @examples
#' anticodoncount()
#'
#' @import seqinr
#'
#' @export
#
anticodoncount <- function(entry=3) {
  require(seqinr)
  if (!exists("myseq")) myseq <- mylist(whatout=1)
  if (!exists("myaa")) myaa <- mylist(whatout=2)
  codons <- NULL
  current_seq  <- myseq[[entry]]
  current_aa <- myaa[[attr(myseq[[entry]], "name")]]
  for (AA in 1:length(myaa[[entry]])) {
    codons <- c(codons, sprintf("%s_%s%s%s", aaa(toupper(current_aa[AA])), toupper(comp(current_seq[(AA-1)*3+3])), toupper(comp(current_seq[(AA-1)*3+2])), toupper(comp(current_seq[(AA-1)*3+1])) ))
  }

  return(table(codons))
}

