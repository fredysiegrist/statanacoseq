#' @title Effective number of codons
#'
#' @description
#' \code{ComputeNEC} is a method for measuring synonymous codon usage bias in a gene by estimation of the “effective number of codons” Nc
#'
#' @details
#' Based on early codon usage mesures proposed by Wright 1990 and Fuglsang 2004. Original implementation from Alexander Roth 2005-2007,
#' second of several indices for codon usage. Nc is the total number of different codons used in a sequence. The values for Nc range from 20 to 61.
#' Where in the first case only one codon is used per amino acid, and in the latter all posible ones are used.
#' Highly expressed genes should use fewer codons due to selection away from equal use. The here computed overall number of effective codons for a gene
#' is a sum of average homozygosities for different redundency classes.
#'
#' @param cds Coding sequence in reading frame
#' @return Numerical value for effective number of codons Nc (by the original definition)
#'
#' @author Roth, A.; Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \url{http://gtrnadb.ucsc.edu/GtRNAdb2/} \code{\link{readstats}}
#' @keywords CodonBias
#' @examples
#' ComputeNEC('ATGTGGTACTCCGACTACGGAGGATAA')
#' ComputeNEC(c2s(mylist(whatout=1)[[1]]))
#'
#' @import seqinr
#'
#' @export
#' @section Original code in Darwin:
#' \subsection{Effective number of codons* (Wright 1990, *Fuglsang 2004).  AR (April 2007)}{\preformatted{
#' ComputeNEC := proc(d:string)
#'  cod:=CreateArray(1..64);
#'  aa:=CreateArray(1..20);
#'  aviod:={op(AToCodon('$'))};
#'  count:=0;
#'  for i to length(d) by 3 do
#'    c := d[i..i+2];
#'    if not member(c,aviod) then
#'      ai:=CodonToInt(c);
#'      ci:=CodonToCInt(c);
#'      cod[ci]:=cod[ci]+1;
#'      aa[ai]:=aa[ai]+1;
#'      count:=count+1;
#'    fi;
#'  od;
#'
#'  Nc:=0;
#'  for i to 20 do
#'    Acods := IntToCInt(i);
#'    k := length(Acods);
#'    if k<2 then Nc := Nc + 1; next; fi;
#'    n := sum([seq(cod[Acods[x]], x=1..k)]);
#'    S := sum([seq((cod[Acods[x]]/n)^2, x=1..k)]);
#'    F := (n*S-1) / (n-1);
#'    Nc := Nc + 1/F;
#'  od;
#'  Nc;
#' end: } }
ComputeNEC <- function(cds) {
  if(!(checkCDS(cds))) {stop("non valid CDS)", call.=FALSE)}
  else {
    cds <- tolower(cds)
    cod <- rep(0, times=64)
    names(cod) <- sapply(as.character(Tef$codons), reversecomplement)
    cod <- cod[c(-59, -60, -64)]
    aa <- rep(0, times=20)
    names(aa) <- levels(Tef[,1])[c(-16,-18)]
    cods <- count(s2c(cds), word = 3, by=3)
    cod[names(cod)] <- as.vector(cods[tolower(names(cod))])
    AA <- table(translate(s2c(cds)))
    aa[names(aa)] <- AA[a(names(aa))]
    NC <- rep(1, times=length(aa))
    names(NC) <- names(aa)
    Acods <- sapply(names(NC), function(i) sapply(as.character(Tef[,2]), reversecomplement)[Tef[,1]==i])
    contributors <- names(aa)[sapply(Acods, length)>1&aa>1]
    Nc <- sapply(contributors, function(i) {
     n <- sum(cod[unlist(Acods[i])])
     S <- sum(sapply(cod[unlist(Acods[i])], function(x) (x/n)^2))
     F <- (n*S-1) / (n-1)
     return(1/F)
     }
    )
    NC[contributors] <- Nc[contributors]
    return(sum(NC))
  }
}
