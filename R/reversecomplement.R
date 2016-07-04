#' @title Reverse Complement Anticodon
#'
#' @description
#' \code{reversecomplement} calculates anticodon from reverse complement of its codon sequence
#'
#' @details
#' Reverse complement codon to isotype / anticodon sequence
#'
#' @param codon 3 character string for codon
#'
#' @return upper case 3 character string for anticodon
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \url{http://trna.ucsc.edu/tRNAscan-SE/}
#' @keywords CodonBias
#' @examples
#' reversecomplement("CAT")
#'
#' @import seqinr
#'
#' @export
reversecomplement <- function(codon) {
  toupper(c2s(rev(comp(s2c(codon)))))
}
