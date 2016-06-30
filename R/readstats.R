#' @title Read tRNA Gene Count Statistics
#'
#' @description
#' \code{readstats} read-in anticodon count from tRNAscan-2.0 or similar summary files from GtRNAdb2 homepage
#'
#'
#' @details
#' Puts gene count into a table according to anticodon sequence on tRNA.
#' This functions works but on the summary files there are some not
#' nicely behaving statistics where the new filtering methods have not been applied as for the fasta files.
#'
#' @param fileno Index number of html-file in the working directory
#' @param chunksize The maximum number of lines of data to be read
#'
#' @return tRNA gene count Table
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \url{http://gtrnadb.ucsc.edu/GtRNAdb2/}
#' @keywords CodonBias
#' @examples
#' setwd("E:/R/viruses")
#' filecount <- length(dir(pattern=".stats.html$"))
#' vir_optcodon <- as.list(rep(data.frame(x=0), times=filecount))
#' for (n in 1:filecount) vir_optcodon[[n]] <- readstats(n)
#' names(vir_optcodon) <- sapply(dir(pattern="stats.html$"), function(x) strsplit(x, "-")[[1]][1])
#'
#' @import seqinr
#'
#' @export

readstats <- function(fileno=156, chunksize=1) {
  try(if(!(fileno %in% 1:length(veab))) {stop("Too few statistic files or non-numerical file number", call.=FALSE)}
      else {
  kingdoms <- rep(c("viruses","eukaryota","archaea","bacteria"), c(1,155,184,4032))
  kingdom <- kingdoms[fileno]
  species <- sapply(strsplit(veab, "-", fixed = TRUE), function(x) x[1])[fileno]
  statfile <- paste("http://gtrnadb.ucsc.edu/GtRNAdb2/genomes", kingdom, species, veab[fileno], sep="/")
  con <- file(statfile, "r", blocking = FALSE) #create file connection
  anticodon_count <- FALSE; pipe <- FALSE
  reps <- as.numeric(); tRNAs <- as.character(); codons <- as.character(); csums <- as.numeric()
  for(i in seq(1,3500,chunksize)){
    d <- scan(con,what="a",nlines=chunksize, quote="", quiet=TRUE)
    if ('</PRE>' %in% d) {
      anticodon_count <- FALSE
      closeAllConnections()
      df <- data.frame(tRNAs, codons, reps, csums)
      optfactor <- function(x){optfactor <- as.numeric(x[3])/max(df[tRNAs==x[1],3])}
      optfact <- apply(df, 1, optfactor)
      isopt <- (optfact == 1)
      fraction <- function(x){fraction <- as.numeric(x[3])/as.numeric(x[4])}
      frac <- apply(df, 1, fraction)
      return(data.frame(tRNAs, codons, reps, optfact, frac, isopt))
    }
    if ( length(d)>0) {
      if (anticodon_count==TRUE) {
        e <- gsub(":", "", d)
        e <- e[nchar(e)>0]
        acvariants <- grep("[[:upper:]]{3}", e)
        for (t in acvariants) {
          if(!is.na(e[t+1]) && length(grep("[[:upper:]]{3}", e[t+1]))==0 ) {
            reps <- c(reps, as.numeric(e[t+1]))
          }
          else {
            reps <-  c(reps, 0)
          }
          tRNAs <- c(tRNAs, e[1])
          codons <- c(codons, e[t])
          csums <- c(csums, e[2])
        }
      }
      if (d[1]=="|") pipe <- TRUE
      if (pipe && d[1]=="Isotype") anticodon_count <- TRUE
    }
  }
} )
}
