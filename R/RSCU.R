#' @title Relative Synonymous Codon Usage
#'
#' @description
#' \code{ComputeNEC} is a method for measuring synonymous codon usage bias in a gene by estimation of the “effective number of codons” Nc
#'
#' @details
#' Should compute the same RSCU as seqinr.
#'
#' @param cds Coding sequence in reading frame
#' @return Named (codons) numerical vector with relative synonymous codon usage for the 64 codons
#'
#' @author Roth, A.; Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \code{\link{readstats}}
#' @keywords CodonBias
#' @examples
#' RSCU('ATGTGGTACTCCGACTACGGAGGATAA')
#' RSCU(c2s(mylist(whatout=1)[[1]]))
#'
#' @import seqinr
#'
#' @export
#' @section Original code in Darwin:
#' \subsection{Relative synonymous codon usage.  AR (2007)}{\preformatted{
#' RSCU := proc(;d:string)
#'   if nargs>0 then
#'     cc:=CodonCount(d);
#'   else
#'     cc := CodonCount();
#'   fi;
#'   rscu := CreateArray(1..64);
#'   for i to 64 do
#'     s:=0;
#'     syn:=IntToCInt(CIntToInt(i));
#'     l:=length(syn);
#'     for j in syn do
#'       s:=s+cc[j];
#'     od;
#'     if s=0 then next fi;
#'     rscu[i]:=cc[i]/(s/l);
#'   od;
#'   rscu;
#' end: } }

RSCU <- function(cds) {
  if(!(checkCDS(cds))) {stop("non valid CDS)", call.=FALSE)}
  else {
    cc <- uco(s2c(cds), 0, "eff")
    rscu <- uco(s2c(cds), 0, "rscu")
    return(rscu)
  }
}
