#' @title Frequency of Optimal Codons
#'
#' @description
#' \code{ComputeFop} calculates ratio of optimal codons used to the sum of synonymous codons
#'
#' @details
#' Based on early codon usage mesures proposed by Ikemura 1981. Original implementation from Alexander Roth 2005-2007,
#' first of several indices for codon usage.
#'
#' @param cds Coding sequence in reading frame
#' @param db CDS list to generate codon usage count table if NULL
#' @param ref 'tRNAgene' takes tRNA gene copy number of isotypes out of genomic data given in codonusagec and 'mostcommon' takes
#' optimal codons as those that are the most common in a database or calculates it from a set of reference genes given as db
#' @param codonusagec Codon usage count table
#' @return Numerical value for ratio of optimal codons (total number of optimal codons / total number of codons included in the analysis)
#'
#' @author Roth, A.; Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \url{http://gtrnadb.ucsc.edu/GtRNAdb2/} \code{\link{readstats}}
#' @keywords CodonBias
#' @examples
#' ComputeFop('ATGTGGTACTCCGACTACGGAGGATAA', ref='tRNAgene', codonusagec=Tef)
#' ComputeFop(c2s(mylist(whatout=1)[[1]]), ref='mostcommon', codonusagec=CodonUsage())
#'
#' @import seqinr
#'
#' @export

ComputeFop <- function(cds, db, ref=c('tRNAgene','mostcommon'), codonusagec=NULL) {
  if(!(checkCDS(cds))) {stop("non valid CDS)", call.=FALSE)}
  else {
   if (ref=='mostcommon' && is.null(codonusagec)) {codonusagec <- CodonUsage(db) }
   try(if((is.null(codonusagec)) && ref=='tRNAgene') {stop("Please indicate tRNAgene counts as codonusagec parameter", call.=FALSE)})
   avoid <- c('TAA', 'TAG', 'TGA', 'ATG', 'TGG') # stop codon $, and amino acids with only one isoacceptor M and W are ommited by definition
   xop<-0; xnon<-0
   for (i in seq(1, nchar(cds), 3)) {
     codon <- substr(cds, i, i+2)
     if (!(codon %in% avoid)) {
       if (codonusagec[codonusagec[,2]==reversecomplement(codon), 4] == 1)  {# define CodonUsage that it behaves like that
         xop <- xop+1
     }
     else {
       xnon <- xnon+1
     }
   }

  }
  return(xop/(xop+xnon))
  }
}
