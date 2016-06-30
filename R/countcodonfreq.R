#' @title Count Codon Frequency
#'
#' @description
#' \code{countcodonfreq} Counts Codon Frequency for a Custom Genome (eg \emph{Eragostis tef})
#'
#'
#' @details
#' This function should count the codon frequency in a genome based on a set of verified and matched transcript and protein fasta files.
#' This is a roxygen2 created helping file.
#'
#' @param no Number of file in working directory to be processed
#' @param wdir Adress pointing to directory with transcript and protein fasta files
#'
#' @return Codon Frequency Table
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}}
#' @keywords CodonBias
#' @examples
#' countcodonfreq()
#'
#' @import seqinr
#'
#' @export
countcodonfreq <- function(no=1, wdir="/windows/R/") {
    require(seqinr)
    list(myseq, myaa) <- mylist(wdir)
    print(length(myseq[[no]])/3/length(myaa[[attr(myseq[[no]], "name")]]))
}
