#' @title My List of Nucleotide-Sequences and Aminoacid-Sequences
#'
#' @description
#' \code{mylist} Read-in of Nucleotide- and Protein-Sequences from Fasta Files of Custom Genome (eg \emph{Eragostis tef})
#'
#' @details
#' Generate a list of myseq and myaa sequences from local fasta files for codon usage analysis. As default some transcripts and corresponding protein sequences from teff are loaded from the external data directory of the package installation.
#'
#' @param RNAfile Directory and File name pointing to fasta file for transcripts
#' @param AAfile Directory and File name pointing to fasta file for corresponding proteins
#' @param whatout Define the type of output; Default 0 means a list of mRNA and protein sequeces; with 1 the function returns only the list of RNA sequences and wiht 2 it returns the proein sequences
#'
#' @return List of sequence list for nucleotide and aminoacid sequences
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{seqinr}} \code{\link{statanacoseq}}
#' @keywords CodonBias
#' @examples
#' mylist()
#'
#' @import seqinr
#'
#' @export
mylist <- function(RNAfile=system.file("extdata",	"Etef.sample.transcript.fasta",	package	=	"statanacoseq"), AAfile=system.file("extdata",	"Etef.sample.protein.fasta",	package	=	"statanacoseq"), whatout=0) {
  require(seqinr)
  if (whatout==0) return( list(myseq <- read.fasta(file = RNAfile, seqtype="DNA"), myaa <- read.fasta(file = AAfile, seqtype="AA") ) )
  if (whatout==1) return( read.fasta(file = RNAfile, seqtype="DNA") )
  if (whatout==2) return( read.fasta(file = AAfile, seqtype="AA") )
}
