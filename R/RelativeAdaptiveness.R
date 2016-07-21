#' @title Compute Relative Adaptiveness
#'
#' @description
#' \code{RelativeAdaptiveness} calculates RA for a List of Entry Numbers
#'
#' @details
#' Should compute the same RA as in Darwin
#'
#' @param entires list of positive integers

#' @return Named (codons) numerical vector with relative adaptiveness for the 64 codons
#'
#' @author Roth, A.; Friberg, M.; Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{ComputeCarboneRA}} \code{\link{statanacoseq}} \code{\link{SetupRA}}
#' @keywords CodonBias
#' @examples
#' RelativeAdaptiveness()
#'
#' @export
RelativeAdaptiveness <- function(entries=list(1)) {
  CodonCounts <- 1:64
#  for i in entries do
#dna := SearchTag('DNA', Entry(i));
#for j to length(dna) by 3 do
#cod := CodonToCInt(dna[j..j+2]);
#if cod=0 then next fi;   # to avoid XXX
#CodonCounts[cod] := CodonCounts[cod]+1;
#od;
#od;
  RA <- 1:64
  aa <- 1
#for aa to 20 do
#codons := IntToCInt(aa);
#counts := [seq(CodonCounts[i], i=codons)];
#freqs := counts / sum(counts);
#for i to length(codons) do
#cod := codons[i];
#RA[cod] := freqs[i] / max(freqs);
#od;
#od;
  #for i to length(RA) do       # set minimum RA value to 0.01
  # if RA[i] = 0 then
  #    RA[i] := 0.01 fi od;
  #for i in AToCInt('$') do     # set RA value of stop codons to 1
  #  RA[i] := 1; od;
  return(RA)
}
