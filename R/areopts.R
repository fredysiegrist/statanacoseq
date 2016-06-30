#' @title Optimality Decision for Codons of Custom Genome
#'
#' @description
#' \code{areopts} Calculates an Optimality Score for Codon in a Protein-Transcript Pair based on tRNA-Gene Abundancy in a Custom Genome (eg \emph{Eragostis tef})
#'
#' @details
#' This function should count the codon frequency in a genome based on a set of verified and matched transcript and protein fasta files.
#' This is a roxygen2 created helping file.
#'
#' @param entry Index number in SeqFastadna list myaa
#' @param opt If boolean and not fraction of tRNA genes to optimal isocode is used
#' @param areopt Vector of anticodon optimality  0>=areopt<=1
#'
#'
#' @return Boolean or numerical vector with Optimality decision or score
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}}
#' @keywords CodonBias
#' @examples
#' areopts(entry=1, opt=TRUE)
#'
#' @import seqinr
#'
#' @export
areopts <- function(entry=3, opt=TRUE) {
  require(seqinr)
  if (!exists("myseq")) myseq <- mylist(whatout=1)
  if (!exists("myaa")) myaa <- mylist(whatout=2)
  codons <- NULL
  current_seq  <- myseq[[entry]]
  current_aa <- myaa[[attr(myseq[[entry]], "name")]]
  for (AA in 1:length(myaa[[entry]])) {
	codons <- c(codons, sprintf("%s_%s%s%s", aaa(toupper(current_aa[AA])), toupper(comp(current_seq[(AA-1)*3+3])), toupper(comp(current_seq[(AA-1)*3+2])), toupper(comp(current_seq[(AA-1)*3+1])) ))
  }
  areopt <- NULL
  for (anticodon in codons) {
    if (opt) {
      areopt <- c(areopt, Tef$isopt[paste(Tef$tRNAs, Tef$codons, sep="_")==anticodon])
    }
    else {
      areopt <- c(areopt, Tef$optfac[paste(Tef$tRNAs, Tef$codons, sep="_")==anticodon])
    }
  }
  return(areopt)
}
