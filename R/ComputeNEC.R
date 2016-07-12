#' @title Effective number of codons
#'
#' @description
#' \code{ComputeNEC} is a method for measuring synonymous codon usage bias in a gene by estimation of the “effective number of codons” Nc
#'
#' @details
#' Based on early codon usage mesures proposed by Wright 1990 and Fuglsang 2004. Original implementation from Alexander Roth 2005-2007,
#' second of several indices for codon usage. Nc is the total number of different codons used in a sequence. The values for Nc range from 20 to 61.
#' Where in the first case only one codon is used per amino acid, and in the latter all posible ones are used.
#' Highly expressed genes should use fewer codons due to selection away from equal use. The here computed overall number of effective codons for a gene
#' is a sum of average homozygosities for different redundency classes.
#'
#' @param cds Coding sequence in reading frame
#' @return Numerical value for effective number of codons Nc (by the original definition)
#'
#' @author Roth, A.; Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \url{http://gtrnadb.ucsc.edu/GtRNAdb2/} \code{\link{readstats}}
#' @keywords CodonBias
#' @examples
#' ComputeNEC('ATGTGGTACTCCGACTACGGAGGATAA')
#' ComputeNEC(c2s(mylist(whatout=1)[[1]]))
#'
#' @import seqinr
#'
#' @export
ComputeNEC <- function(cds) {
  if(!(checkCDS(cds))) {stop("non valid CDS)", call.=FALSE)}
  else {
    cds <- toupper(cds)
    cod <- rep(0, times=64)
    names(cod) <- sapply(as.character(Tef$codons), reversecomplement)
    cod <- cod[c(-59, -60, -64)]
    aa <- rep(0, times=20)
    names(aa) <- levels(Tef[,1])[c(-16,-18)]
    avoid <- c('TAA', 'TAG', 'TGA') # stop codons $ are ommited by definition
    count <- 0
    for (i in seq(1, nchar(cds), 3)) {
      codon <- substr(cds, i, i+2)
      if (!(codon %in% avoid)) {
          cod[codon] <- cod[codon]+1
          aa[aaa(translate(s2c(codon)))] <- aa[aaa(translate(s2c(codon)))]+1
          count <- count+1
      }
    }
    Nc <- 0
    for (i in names(aa)) {
      Acods <- sapply(as.character(Tef[,2]), reversecomplement)[Tef[,1]==i]
      k <- length(Acods)
      if (k<2) {
        Nc <- Nc+1
      }
      else {
        n <- sum(cod[Acods])
        S <- sum(sapply(cod[Acods], function(x) (x/n)^2))
        F <- (n*S-1) / (n-1)
        Nc <- Nc+(1/F)
      }
    }
    return(Nc)
  }
}


#ComputeNEC2 <- function(cds) {
#  if(!(checkCDS(cds))) {stop("non valid CDS)", call.=FALSE)}
#  else {
#    cod <- count(cd, word = 3, by=3)
#    aa <- table(translate(s2c(cds)))
#  }
#}
