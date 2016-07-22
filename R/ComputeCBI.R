#' @title Codon Bias Index
#'
#' @description
#' \code{ComputeCBI} Computes the Codon Bias Index for a Gene after Bennetzen and Hall 1982
#'
#' @details
#' Should be implemented because missing in Darwin
#'
#' @param cds Coding Sequence
#' @return VOID
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}}
#' @keywords CodonBias
#' @examples
#' ComputeCBI('ATGTGGTACTCCGACTACGGAGGATAA')
#' ComputeCBI(c2s(mylist(whatout=1)[[1]]))
#'
#' @import seqinr
#'
#' @export
#' @section Original code in Darwin:
#' \subsection{Codon Bias Index (CBI) (Bennetzen and Hall 1982).}{\preformatted{
#' ComputeCBI :=  proc(d:string)
#' # to be implemented
#' end: } }
  ComputeCBI <- function(cds) {
  # to be implemented
}
