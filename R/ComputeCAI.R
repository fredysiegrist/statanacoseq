#' @title Codon Adaptation Index
#'
#' @description
#' \code{ComputeCAI} Computes Codon Adaptation Index
#'
#' @details
#' Should compute the same ComputeCAI or ComputeCAIVector as in Darwin and the cai of seqinr package.
#'
#' @param cds Coding sequence in reading frame
#' @param RA Relative Adaptiveness table to use
#' @param UseCodonProb Wheter to use Codon Probabilities (default = FALSE)
#' @return Numerical with CAI
#'
#' @author Roth, A.; Friberg, M.; Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link[seqinr]{cai}} \code{\link{statanacoseq}} \code{\link{checkCDS}}
#' @keywords CodonBias
#' @examples
#' ComputeCAI('ATGTGGTACTCCGACTACGGAGGATAA', RA=SetupRA("yeast"), UseCodonProb=TRUE)
#' ComputeCAI(mylist(whatout=1)[[1]], RA=ComputeCarboneRA(DB=mylist(whatout=1)))
#'
#' @import seqinr
#'
#' @export
#' @section Original code in Darwin:
#' \subsection{Compute CAI, the Codon Adaptation Index (Sharp and Li 1987)}{\preformatted{
#' # Markus Friberg and Alexander Roth (Dec 2005)
#' ComputeCAI := proc(DNA:{string, Entry})
#'  # check global variables and scan arguments
#'  if not assigned(RA) then
#'    error('Error in ComputeCAI: RA not assigned, use e.g. SetupRA(yeast);') fi;
#'  if type(DNA, Entry) then dna:=copy(SearchTag('DNA', DNA))
#'  else dna:=DNA fi;
#'  UseCodonProb := false;
#'  for i from 2 to nargs do
#'    if length(args[i]) = 2 and args[i, 1] = 'UseCodonProb' then
#'      UseCodonProb := args[i, 2]
#'    else
#'      error('Unknown argument ', args[i]);
#'    fi;
#'  od;
#'  # compute cai
#'  w := 0;
#'  n := length(dna)/3;
#'  for j to length(dna) by 3 do
#'    cint := CodonToCInt(dna[j..j+2]);
#'    codprob := If(UseCodonProb, CodonProb[cint], 1);
#'    if CIntToA(cint) <> '$' then    # don't consider stop codons
#'    w := w + ln(codprob * RA[cint]) fi;
#'  od;
#'  exp(1/n * w)
#' end: } }
ComputeCAI <- function(cds, RA, UseCodonProb=FALSE) {
  dna <- toupper(c2s(cds))
  if(!(checkCDS(dna))) {stop("non valid CDS)", call.=FALSE)}
  else {
    if (missing(RA)) {
      stop('Error in ComputeCAI: RA not assigned, use e.g. RA=SetupRA("yeast")', call. = FALSE)
    }
    cai <- NA
    w <- 0
    n <- nchar(dna)/3;
    for (j in seq(from=1, to=nchar(dna), by=3)) {
      codon <- reversecomplement(substr(dna, j, j+2))
      cint <- which(substr(names(aa_ac), 5, 7) %in% codon)
      if (UseCodonProb) {codprob <- CodonProb(cint)}
      else {codprob <- 1}
      if (!(cint %in% c(59,60,64,65)))  {  # don't consider stop codons
        w <- w + log(codprob * RA[cint])
      }
    }
    cai <- exp(1/n * w)
    return(unname(cai))
  }
}
