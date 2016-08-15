#' @title GC content 1st and 2nd position of synonymous codons
#'
#' @description
#' \code{ComputeGC3syn} Computes the G+C content at the 1st and 2nd position of synonymous codons
#'
#' @details
#' Computes the GC content of the 1st and 2nd position of synonymous codons to estimate the background
#' GC content when compared with the GC content at the 3rd position. By definition, GC12s values are
#' the proportion of GC nucleotides at the fixed first and less variable second coding position of
#' synonymous codons. This can be used to evaluate the degree of base composition bias.
#'
#' @param tD Nucleotide character string
#' @return vector of (o/n)first (o/n)second
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \code{\link{ComputeGC3syn}}
#' @keywords CodonBias
#' @examples
#' ComputeGC12syn('ATGTGGTACTCCGACTACGGAGGATAA')
#' plot(mean(ComputeGC12syn(toupper(c2s(mylist(whatout=1)[[1]])))), ComputeGC3syn(toupper(c2s(mylist(whatout=1)[[1]]))))
#'
# list arguments have to be properly implemented
# ComputeGC12syn(c2s(mylist(whatout=1)[[1]]))
#'
#' @import seqinr
#'
#' @export
ComputeGC12syn <- function(tD) {
  require(seqinr)
  if(!(checkCDS(tD))) stop("non valid CDS)", call.=FALSE)
  # remove stop codon
  if (c2s(s2c(tD)[(nchar(tD)-2):nchar(tD)]) %in% substr(names(aa_ac), 5, 7)[c(59, 60, 64, 65)]) {
    d <- substring(tD, 1, nchar(tD)-2)
  }
  else {
    d <- tD
  }
  bases <- c("A", "C", "G", "T")
  o1 <- rep(0, times=4)
  names(o1) <- bases
  o2 <- o1
  n <- 0
  for (i in seq(from=1, to=nchar(d), by=3)) {
    c <- substring(d, i, i+2)
    if (c %in% substr(names(aa_ac), 5, 7)[c(-59, -60, -64, -65)]) {
      n <- n+1
      oi1 <- which(names(o1) %in% substring(c, 1, 1))
      oi2 <- which(names(o2) %in% substring(c, 2, 2))
      o1[oi1] <- o1[oi1]+1
      o2[oi2] <- o2[oi2]+1
    }
  }
  o1 <- o1/n
  o2 <- o2/n
  return(c(o1[2]+o1[3], o2[2]+o2[3]))
}
