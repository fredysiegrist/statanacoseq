#' @title Checks Integrity of CDS Sequence
#'
#' @description
#' \code{checkCDS} Checks that the character number in the cds is in codons (modulo 3 = zero)
#'
#' @details
#' Based on early codon usage mesures proposed by Ikemura 1981. Original implementation from Alexander Roth 2005-2007,
#' first of several indices for codon usage.
#'
#' @param cds Sequence to be checked     # may be implemented as list List of Coding Sequences as read.fasta input
#' @param ignore If TRUE returns TRUE even in case of biological violations (as for stop codons in sequence).
#' @return TRUE if all checks are passed with or without warnings, FALSE if an violation is found
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \code{\link{ComputeFop}}
#' @keywords CodonBias
#' @examples
#' checkCDS(c2s('ATGTGGTACTCCGACTACGGAGGATAA'))
#'
#' @import seqinr
#'
#' @export

checkCDS <- function(cds, ignore=FALSE) {
  cdslength <- nchar(cds)
  # madatory check if cds is dividable in codons
  if(!(cdslength %% 3 == 0)) {warning("CDS is not splittable in codons (blocks of 3 characters)"); return(FALSE)}
  # optional check if cds start with the start codon
  if(!(toupper(substring(cds, 1, 3))=="ATG")) {warning("CDS does not start with ATG)")}
  # optional check if cds ends on a stop codon
  if(!(toupper(substring(cds, cdslength-2, cdslength)) %in% c('TAA', 'TAG', 'TGA'))) {warning("CDS does not end with stop codon)")}
  # mandatory check if cds contains stop codons
  if("*" %in% getTrans(s2c(cds))[1:((cdslength/3)-1)]) {warning("CDS is not a Open Reading Frame"); return(ignore)}

  # if list of DNA and AA is given it may check for translation errors (mismatches)
  # that the number of the codons represent length of protein sequence
  # that the stop codon is correctly coded and internal stop codons are at the end of the sequence and translated to an amino acid

  return(TRUE)
}
