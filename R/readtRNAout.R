#' @title Read tRNA Gene Count Statistics
#'
#' @description
#' \code{readtRNAout} read-in anticodon count from tRNAscan-2.0 output file from online web tool
#'
#'
#' @details
#' Puts gene count into a table according to anticodon sequence on tRNA.
#' This functions works on a tRNA-scan output reanalyzed by tRNAscan-SE 2.0 online tool
#' to behave similarly to the data available on GtRNAdb2.
#'
#' @param file seq.out file from tRNAscan-SE 2.0
#'
#' @return tRNA gene count Table
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \url{http://trna.ucsc.edu/tRNAscan-SE/}
#' @keywords CodonBias
#' @examples
#' readtRNAout()
#'
#' @import seqinr
#'
#' @export

readtRNAout <- function(file=system.file("extdata",	"seq30469.out",	package	=	"statanacoseq")) {
  try(if(!file.exists(file)) {stop("No such file or directory.", call.=FALSE)}
      else {
        tRNAscan.out <- read.table(file=file, skip=3)
        tRNAscan.out[,13] <- paste(substr(tRNAscan.out[,5], 1, 3), gsub("N", "?",tRNAscan.out[,6]), sep="_")
        reps <- aa_ac
        tRNAs <-  substr(names(reps), 1, 3)
        codons <- unlist(sapply(names(reps), function(x) substring(x, nchar(x)-2, nchar(x))))
        for(i in names(aa_ac) ) {
          reps[which(names(reps) %in% i)] <- sum(tRNAscan.out[,13]==i)
        }
        df <- data.frame(tRNAs=substr(names(reps), 1, 3), codons=substr(names(reps), nchar(names(reps))-2, nchar(names(reps))), reps)
        optfactor <- function(x){optfactor <- as.numeric(x[3])/max(df[tRNAs==x[1],3])}
        optfact <- apply(df, 1, optfactor)
        isopt <- (optfact == 1)
        fraction <- function(x){fraction <- as.numeric(x[3])/sum(df[tRNAs==x[1],3])}
        frac <- apply(df, 1, fraction)
        return(data.frame(df, optfact, frac, isopt))
      } )
}
