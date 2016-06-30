#' @title Read tRNA Intronic Gene Count Statistics
#'
#' @description
#' \code{readintronic} read-in anticodon count from tRNAscan-2.0 or similar summary files from GtRNAdb2 homepage
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
#' @return tRNA gene with introns count Table
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \url{http://gtrnadb.ucsc.edu/GtRNAdb2/}
#' @keywords CodonBias
#' @examples
#' filecount <- length(dir(pattern=".stats.html$"))
#' vir_optcodon <- as.list(rep(data.frame(x=0), times=filecount))
#' for (n in 1:filecount) vir_optcodon[[n]] <- readstats(n)
#' names(vir_optcodon) <- sapply(dir(pattern="stats.html$"), function(x) strsplit(x, "-")[[1]][1])
#'
#' @import seqinr
#'
#' @export

readintronic <- function(fileno=156, chunksize=1) {
  try(if(!(fileno %in% 1:length(veab))) {stop("Too few statistic files or non-numerical file number", call.=FALSE)}
      else {
        kingdoms <- rep(c("viruses","eukaryota","archaea","bacteria"), c(1,155,184,4032))
        kingdom <- kingdoms[fileno]
        species <- sapply(strsplit(veab, "-", fixed = TRUE), function(x) x[1])[fileno]
        statfile <- paste("http://gtrnadb.ucsc.edu/GtRNAdb2/genomes", kingdom, species, veab[fileno], sep="/")

  con <- file(statfile, "r", blocking = FALSE) #create file connection
  #d=scan(con,what="a",nlines=1,sep="|") #remove the header line

  for(i in seq(1,3500,chunksize)){

    d=scan(con,what="a",nlines=chunksize,sep="|", quote="",quiet=TRUE)
    #    d = t(matrix(d,nrow=11))
    #    d = data.frame(d)
    #Do stuff with d....

    if (length(d)>0 && d[1]==""){
      print(substr(d[2:(length(d)-1)], 2, 15) )
      reps <- as.numeric(gsub("[^0-9]","", d[2:(length(d)-1)]))
      tRNAs <- substr(d[2:(length(d)-1)], 2, 4)
      codons <- substr(d[2:(length(d)-1)], 6, 8)
      closeAllConnections()
      df <- data.frame(tRNAs, codons, reps)
      optfactor <- function(x=c("Val","UUU","1")){optfactor <- as.numeric(x[3])/max(df[tRNAs==x[1],3])}
      optfact <- apply(df, 1, optfactor)
      isopt <- (optfact == 1)
      return(data.frame(tRNAs, codons, reps, optfact, isopt))
    }
  }
} )}
