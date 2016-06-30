#' Optimal Codon Table for Eragostis tef Decided on tRNAscan-SE count of tRNA genes in sequenced genome
#'
#' A table with the number of tRNA genes in \emph{Eragostis tef} for optimal codon evaluation.
#'
#' @format A data frame with 64 rows (anticodons) and 6 descriptive variables:
#' \describe{
#'   \item{tRNAs}{Aminoacid coded by tRNA}
#'   \item{codons}{anticodon on tRNA gene}
#'   \item{reps}{Cove-confirmed tRNAs}
#'   \item{optfact}{# isotype (isoacceptor typ) / max(#) isotype for aminoacid}
#'   \item{frac}{# isotype / sum(#) isotype for aminoacid}
#'   \item{isopt}{boolean decision if it is (one of) the optimal codon - optfract == 1}
#' }
#' @source \url{http://gtrnadb.ucsc.edu/GtRNAdb2/}
#'
"Tef"
#' Viruses, Eukaryota, Archaea, Bacteria Species File Names on GtRNAdb2
#'
#' A Vector with .html File Names Containing the tRNAscan-SE v.2.0 Analysis of Complete Genomes from GtRNAdb
#'
#' @format A vector of stings (file names)
#' @source \url{http://gtrnadb.ucsc.edu/GtRNAdb2/}
#'
"veab"
#' Viruses, Eukaryota, Archaea, Bacteria Species Fasta File Names on GtRNAdb2
#'
#' A Vector with .fa File Names Containing the tRNA Gene List of Complete Genomes from GtRNAdb
#'
#' @format A vector of stings (file names)
#' @source \url{http://gtrnadb.ucsc.edu/GtRNAdb2/}
#'
"veabfa"
#' Viruses, Eukaryota, Archaea, Bacteria Species Names on GtRNAdb2
#'
#' A Vector with Species Names from GtRNAdb, edited for html conformatibility
#'
#' @format A vector of stings (partial file names)
#' @source \url{http://gtrnadb.ucsc.edu/GtRNAdb2/}
#'
"GtRNAdb2species"
#' Aminoacids and Anticodons
#'
#' Named Vector with Aminoacid and Anticodon Tokens with empty (0) counts
#'
#' @format Named (aminoacid and anticodon) vector of numerical (tRNA gene count) set to 0
#' @source \url{http://www.tef-research.org/}
#'
"aa_ac"
