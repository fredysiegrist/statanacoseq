#' @title Search and Return a Fasta File Name from Partial Name
#'
#' @description
#' \code{searchfafile} searches fasta file name on GtRNAdb2 homepage
#'
#'
#' @details
#' Import helper function to manage automatical downloading of fasta files.
#'
#' @param fam Taxonomy kingdom defining the folder on gtrnadb homepage
#' @param spec Species identifier as pointer to .fa filename
#'
#' @return string with .fa file name
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \url{http://gtrnadb.ucsc.edu/GtRNAdb2/}
#' @keywords CodonBias
#' @examples
#' for (i in virus) try(download.file(paste("http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/viruses/", i, "/", searchfafile(fam='viruses', spec=i), sep=""), paste("viruses/", i, "-tRNAs.fa", sep="") ))
#'
#' @import seqinr
#'
#' @export
searchfafile <- function(fam='eukaryota', spec='Hsapi19') {
  specfolder <- gsub("#", "%23", spec)
  # idxfile <- paste("http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/", fam, "/", specfolder, "/", "index.html", sep="")
  statfile <- paste("http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/", fam, "/", specfolder, "/", gsub(c("#"), "%23", veab[which(GtRNAdb2species %in% spec)]), sep="")
  con <- file(statfile, "r", blocking = FALSE) #create file connection
  d=scan(con, what="a", nlines=200, sep="<", quiet=TRUE)
  sline <- d
  pattern <-  '[-A-Za-z0-9%*_#+\\.]*[.]fa'
  m <- regexpr(pattern, sline)
  faname <- gsub("#", "%23", regmatches(sline, m))
  closeAllConnections()
  return(faname)
}
