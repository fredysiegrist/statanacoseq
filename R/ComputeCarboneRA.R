#' @title Compute Carbone Relative Adaptiveness
#'
#' @description
#' \code{ComputeCarboneRA} calculates RA
#'
#' @details
#' Should compute the same RA in Darwin
#'
#' @param t nonnegative
#' @param initfrac nonnegative
#' @param interfrac nonnegative
#' @param mode "reverse" for standard anything else for non-standard
#' @return Named (codons) numerical vector with relative synonymous codon usage for the 64 codons
#'
#' @author Roth, A.; Friberg, M.; Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{RelativeAdaptiveness}} \code{\link{statanacoseq}} \code{\link{SetupRA}}
#' @keywords CodonBias
#' @examples
#' ComputeCarboneRA()
#'
#' @export
#'
ComputeCarboneRA <- function(t=0.01, initfrac=1, iterfrac=0.5, mode='reverse') {
  RA <- NA
  #if not assigned(DB) then error('DB must be assigned') fi;
  x <- 1  # fraction of the sequences used to compute RA in this iteration
  #AllGenes := [seq(i, i=1..DB[TotEntries])]:
  #genes := Shuffle(AllGenes)[1..round(initfrac * DB[TotEntries])]:
  bestCorr <- 0
  #cai := CreateArray(1..DB[TotEntries]):
  #while length(genes) / DB[TotEntries] > t do
  #RA := RelativeAdaptiveness(genes);
  #for i to DB[TotEntries] do
  #dna:=SearchTag('DNA',Entry(i));
  #if SearchString('X', dna)<>-1 then next fi;
  #cai[i] := ComputeCAI(dna) od;
  #x := x * iterfrac;
  #res := transpose([AllGenes, cai]):
  #if mode='reverse' then
  #res := transpose(sort(res, res -> res[2])):
  #else
  #  res := transpose(sort(res, res -> -res[2])):
  #fi;
  #genes := res[1][1..round(x * DB[TotEntries])]:
  #od;
  return(RA)
}
