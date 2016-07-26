#' @title Table for Codon Usage and Optimality of Codons
#'
#' @description
#' \code{CodonUsage} calculates ratio of optimal codons used to the sum of synonymous codons
#'
#' @details
#' Creates a count table to decide on optimality of codons based on the usage of the imput DNA sequences.
#'
#' @param cdslist List of coding sequences in reading frame (eg. of full genome), handles also DNA string sequences.
#' @return data frame with codon usage count Table with statistics for optimality based on occurence
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{statanacoseq}} \code{\link{ComputeFop}}
#' @keywords CodonBias
#' @examples
#' CodonUsage('ATGTGGTACTCCGACTACGGAGGATAA')
#' CodonUsage(mylist(whatout=1))
#'
#' @export
CodonUsage <- function(cdslist=mylist()) {
  if (!(is.list(cdslist))) {
    if(checkCDS(cdslist)) {
      cdslist <- s2c(cdslist)
    }
  }
  ucoc <- uco(unlist(cdslist), frame = 0, index = "eff", as.data.frame = FALSE, NA.rscu = 0)
  reps <- aa_ac
  tRNAs <-  substr(names(reps), 1, 3)
  codons <- sapply(1:65, function(x) reversecomplement(unlist(sapply(names(aa_ac), function(x) substring(x, nchar(x)-2, nchar(x))))[x]) )
  counts <- c(ucoc[tolower(codons[1:64])],"???"=0)
  df <- data.frame(tRNAs=tRNAs, codons=codons, counts)
  optfactor <- function(x){optfactor <- as.numeric(x[3])/max(df[tRNAs==x[1],3])}
  rownames(df) <- 1:65
  optfact <- apply(df, 1, optfactor)
  isopt <- (optfact == 1)
  fraction <- function(x){fraction <- as.numeric(x[3])/sum(df[tRNAs==x[1],3])}
  frac <- apply(df, 1, fraction)
  return(data.frame(df, optfact, frac, isopt))
}
