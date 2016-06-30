#' @title Processed tRNA Containing Fasta File Headers to Count Anticodons
#'
#' @description
#' \code{readfasta} processes fasta sequence headers from tRNA gene sequences of complete genomes and compares anticodon usage in the gene set.
#'
#'
#' @details
#' Import helper function to manage automatical downloading of fasta files.
#'
#' @param fam Taxonomy kingdom defining the folder on gtrnadb homepage
#' @param spec Species identifier as pointer to .fa filename
#' @param verbose Wheter to show the aminoacid, anticodon and tRNA Model Score from fasta entry header
#' @return tRNA gene count Table with statistics
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \url{http://gtrnadb.ucsc.edu/GtRNAdb2/} \code{\link{readstats}}
#' @keywords CodonBias
#' @examples
#' readfasta(131)
#'
#' @import seqinr
#'
#' @export
readfasta <- function(fileno=1, chunksize=1, verbose=TRUE) {
  # decide on directory on server
  try(if(!(fileno %in% 1:length(GtRNAdb2species))) {stop("Too few statistic files or non-numerical file number", call.=FALSE)}
      else {
        kingdoms <- rep(c("viruses","eukaryota","archaea","bacteria"), c(1,155,184,4032))
        kingdom <- kingdoms[fileno]
        species <- gsub("#", "%23", sapply(strsplit(veab, "-stats", fixed = TRUE), function(x) x[1])[fileno])
        fafile <- paste("http://gtrnadb.ucsc.edu/GtRNAdb2/genomes", kingdom, species, searchfafile(fam=kingdom, spec=gsub("#", "%23", GtRNAdb2species[fileno])), sep="/")

  # load file from GtRNAdb2 server
  con <- file(fafile, "r", blocking = FALSE) #create file connection
  reps <- aa_ac
  tRNAs <-  substr(names(reps), 1, 3)
  codons <- unlist(sapply(names(reps), function(x) substring(x, nchar(x)-2, nchar(x))))
  for(i in seq(1,10000,chunksize)){
    d=scan(con,what="a",nlines=chunksize, quote="",quiet=TRUE)
    if (length(d)>0 && substr(d[1], 1, 1)==">"){
      rgno <- grep("[(][?A-Z]{3}[)]", d)
      if (!verbose) print(paste(d[rgno-1], substr(d[rgno], 2, 4), d[rgno+4]) )
      if (paste(substr(d[rgno-1], 1, 3), substr(d[rgno], 2, 4), sep='_') %in% names(reps)) {
        reps[paste(substr(d[rgno-1], 1, 3), substr(d[rgno], 2, 4), sep='_')] <- reps[paste(substr(d[rgno-1], 1, 3), substr(d[rgno], 2, 4), sep='_')]+1
      }
      else { # count Unk ??? one up also for Und NNN etc
        reps[65] <- reps[65]+1
      }
    }
  }
  closeAllConnections()
  df <- data.frame(tRNAs=substr(names(reps), 1, 3), codons=substr(names(reps), nchar(names(reps))-2, nchar(names(reps))), reps)
  optfactor <- function(x){optfactor <- as.numeric(x[3])/max(df[tRNAs==x[1],3])}
  optfact <- apply(df, 1, optfactor)
  isopt <- (optfact == 1)
  fraction <- function(x){fraction <- as.numeric(x[3])/sum(df[tRNAs==x[1],3])}
  frac <- apply(df, 1, fraction)
  return(data.frame(df, optfact, frac, isopt))
} )

}
