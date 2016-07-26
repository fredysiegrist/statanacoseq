#' @title Compute Carbone Relative Adaptiveness
#'
#' @description
#' \code{ComputeCarboneRA} calculates RA
#'
#' @details
#' Should compute the same RA in Darwin
#'
#' @param target nonnegative
#' @param initfrac nonnegative
#' @param interfrac nonnegative
#' @param reverse For standard anything else for non-standard
#' @param DB Database with working genes of type Entires()
#' @return Named (codons) numerical vector with relative synonymous codon usage for the 64 codons
#'
#' @author Roth, A.; Friberg, M.; Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{RelativeAdaptiveness}} \code{\link{statanacoseq}} \code{\link{SetupRA}}  \code{\link{Entries}}
#' @keywords CodonBias
#' @examples
#' ComputeCarboneRA(DB=mylist(whatout=1))
#'
#' @export
#' @section Original code in Darwin:
#' \subsection{Compute Carbone Relative Adaptiveness}{\preformatted{
#' ComputeCarboneRA := proc( ; t=0.01:nonnegative, initfrac=1:nonnegative, iterfrac=0.5:nonnegative, mode:string)
#'  global RA;
#'  if not assigned(DB) then error('DB must be assigned') fi;
#'  x := 1;  # fraction of the sequences used to compute RA in this iteration
#'  AllGenes := [seq(i, i=1..DB[TotEntries])]:
#'  genes := Shuffle(AllGenes)[1..round(initfrac * DB[TotEntries])]:
#'  bestCorr := 0;
#'  cai := CreateArray(1..DB[TotEntries]):
#'  while length(genes) / DB[TotEntries] > t do
#'    RA := RelativeAdaptiveness(genes);
#'    for i to DB[TotEntries] do
#'      dna:=SearchTag('DNA',Entry(i));
#'      if SearchString('X', dna)<>-1 then next fi;
#'      cai[i] := ComputeCAI(dna) od;
#'    x := x * iterfrac;
#'    res := transpose([AllGenes, cai]):
#'    if mode='reverse' then
#'      res := transpose(sort(res, res -> res[2])):
#'    else
#'      res := transpose(sort(res, res -> -res[2])):
#'    fi;
#'    genes := res[1][1..round(x * DB[TotEntries])]:
#'  od;
#' RA
#' end: } }
ComputeCarboneRA <- function(target=0.05, initfrac=1, iterfrac=0.5, reverse=FALSE, DB) {
  if (!(is.double(target) && length(target)==1 && target>(1/length(DB)) && target<1)) stop('target must be a positive number between 1/length(DB) and 1', call.=FALSE)
  RA <- NA
  if(missing(DB)) stop('DB must be assigned', call.=FALSE)
  x <- 1  # fraction of the sequences used to compute RA in this iteration
  AllGenes <- 1:length(DB)
  genes <- sample(AllGenes, round(initfrac * length(DB)))
  cai <- rep(NA, times=length(DB))
  while((length(genes) / length(DB)) > target) {
    suppressWarnings(RA <- RelativeAdaptiveness(DB[sapply(genes, function(j) attr(DB[[j]], "name"))]))
    for (i in 1:length(DB)) {
      dna <- c2s(DB[[i]])
      if (length(grep("n", dna)) == 1 | length(grep("x", dna)) == 1 ) break
      suppressWarnings(cai[i] <- ComputeCAI(dna, RA=RA))
    }
    x <- x * iterfrac
    res <- rbind(AllGenes, cai)
    res[1,order(mat[2,], decreasing = !reverse)]
    genes <- res[1, 1:round(x * length(DB))]
  }
  return(RA)
}
