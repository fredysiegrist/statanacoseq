#' @title Find Highly Expressed Genes
#'
#' @description
#' \code{FindHighlyExpressedGenes} Finds Highly Expressed Genes in Lost
#'
#' @details
#' Should work the same as FindHighlyExpressedGenes() in Darwin on DB
#'
#' @param n Integer
#' @param tag String either 'PROTEXPR' or 'MRNAEXPR'
#'
#' @return VOID
#'
#' @author Roth, A.; Friberg, M.; Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}} \code{\link{CodonProbabilities}}
#' @keywords CodonBias
#' @examples
#' FindHighlyExpressedGenes()
#'
#' @import seqinr
#'
#' @export
#' @section Original code in Darwin:
#' \subsection{Compute CAI, the Codon Adaptation Index (Sharp and Li 1987)}{\preformatted{
#' FindHighlyExpressedGenes := proc( ; n=100:integer, tag='PROTEXPR':string)
#'  # tags: 'PROTEXPR' 'MRNAEXPR'
#'  expr := CreateArray(1..DB[TotEntries]);
#'  for i to DB[TotEntries] do
#'    ex := sscanf(SearchTag(tag, Entry(i)), '%f');
#'    if ex <> [] then
#'      expr[i] := op(ex) fi;
#'  od;
#'  sorted := sort(expr);
#'  limit := sorted[length(sorted)-n+1];
#'  genes := [];
#'  for i to DB[TotEntries] do
#'    if expr[i] >= limit then
#'      genes := append(genes, i) fi od;
#'    genes
#' end: } }
FindHighlyExpressedGenes <- function(n=100, tag='PROTEXPR') {
  expr <- 1:length(mylist(whatout=1))
  if (n > length(mylist(whatout=1))) n <- length(mylist(whatout=1))
  for (i in mylist(whatout=1)) {
    # NEEDS IMPLEMENTATION OF RNAseq DATA!
#    ex <- sscanf(SearchTag(tag, Entry(i)), '%f');
#    if ex <> [] then
#  expr[i] := op(ex) fi;
  }
  sorted <- sort(expr)
  limit <- sorted[length(sorted)-n+1]
  genes <- NULL
  for (i in 1:length(mylist(whatout=1))) {
    if (expr[i] >= limit) {
      genes <- c(genes, substring(attr(mylist(whatout=1)[[i]], 'Annot'),2))
    }
  }
  return(genes)
}
