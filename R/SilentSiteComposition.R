#' @title Base Composition at Silent Sites
#'
#' @description
#' \code{SilentSiteComposition} Computes the Base Composition at Silent Sites
#'
#' @details
#' Should be implemented because missing in Darwin
#'
#' @param cds Coding Sequence for Silent Site Composition calculation
#' @return VOID
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}}
#' @keywords CodonBias
#' @examples
#' SilentSiteComposition('ATGTGGTACTCCGACTACGGAGGATAA')
#' SilentSiteComposition(c2s(mylist(whatout=1)[[1]]))
#'
#' @import seqinr
#'
#' @export
#' @section Original code in Darwin:
#' \subsection{Base composition at silent sites.}{\preformatted{
#' SilentSiteComposition := proc(d:string)
#' # to be implemented
#' end: } }
SilentSiteComposition <- function(cds) {
  # to be implemented
}
