#' @title Codon Adaptation Index Vector
#'
#' @description
#' \code{ComputeCAIVector} Computes Codon Adaptation Index
#'
#' @details
#' Should compute the same ComputeCAI or ComputeCAIVector as in Darwin and the cai of seqinr package.
#'
#' @param cds Coding sequence in reading frame
#' @return Vector of all 65 codon values
#'
#' @author Roth, A.; Friberg, M.; Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link[seqinr]{cai}} \code{\link{statanacoseq}} \code{\link{checkCDS}}
#' @keywords CodonBias
#' @examples
#' ComputeCAIVector('ATGTGGTACTCCGACTACGGAGGATAA', RA=SetupRA("yeast"))
#' ComputeCAIVector(c2s(mylist(whatout=1)[[1]]), RA=ComputeCarboneRA(DB=mylist(whatout=1)))
#'
#' @import seqinr
#'
#' @export
#' @section Original code in Darwin:
#' \subsection{Compute CAI, the Codon Adaptation Index (Sharp and Li 1987)}{\preformatted{
#' ComputeCAIVector := proc(e:Entry)
#'  if not assigned(RA) then
#'    error('Error in ComputeCAI: RA not assigned, use e.g. SetupRA(yeast);') fi;
#'  dna := SearchTag('DNA', e);
#'  wa := CreateArray(1..20);
#'  na := CreateArray(1..20);
#'  for j to length(dna) by 3 do
#'    cint := CodonToCInt(dna[j..j+2]);
#'    a := CIntToInt(cint);
#'    if a <= 20 then
#'      wa[a] := wa[a] + ln(RA[cint]);
#'      na[a] := na[a]+1;
#'    fi;
#'  od;
#'  res := CreateArray(1..21);
#'  for i to 20 do
#'    res[i] := If(na[i]=0, 'NA', exp(1/na[i] * wa[i])) od;
#'  res[21] := exp(1/sum(na) * sum(wa));
#'  res
#' end: } }

ComputeCAIVector <- function(cds, RA) {
  if(!(checkCDS(cds))) {stop("non valid CDS)", call.=FALSE)}
  else {
    if (missing(RA)) {
      stop('Error in ComputeCAI: RA not assigned, use e.g. RA=SetupRA("yeast")', call. = FALSE)
    }
    a.uni <- unique(a(substr(names(aa_ac), 1, 3)[c(-59, -60, -64, -65)]))
    dna <- toupper(cds)
    if(!(checkCDS(dna))) stop("non valid CDS)", call.=FALSE)
    wa <- rep(0, times=20)
    na <- rep(0, times=20)
    for (j in seq(from=1, to=nchar(dna), by=3)) {
      codon <- reversecomplement(substr(dna, j, j+2))
      cint <- which(substr(names(aa_ac), 5, 7) %in% codon)
      if (cint == 65) break # to avoid XXX
      a <- which(a.uni %in% a(substr(names(aa_ac), 1, 3)[cint]))
      if (a <= 20) {
        wa[a] <- wa[a] + log(RA[cint])
        na[a] <- na[a]+1
      }
    }
    res <- rep(NA, times=21)
    for (i in 1:20) {
      if(!(na[i]==0)) res[i] <- exp(1/na[i] * wa[i])
    }
    res[21] <- exp(1/sum(na) * sum(wa))
    return(res)
  }
}
