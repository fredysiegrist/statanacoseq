#' @title Codon Probability table
#'
#' @description
#' \code{CodonProb} creates a table of codon probabilities
#'
#' @details
#' Generation of manually set tables from literature. Maybe enhanced to compute for the used database.
#'
#' @param pos Which Codon Probabilities should be printed out
#' @return Numerical vector with relative adaptiveness for the 64 codons
#'
#' @author Roth, A.; Friberg, M.; Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{statanacoseq}} \code{\link{SetupRA}}
#' @keywords CodonBias
#' @examples
#' CodonProb(1:64)
#'
#' @export

CodonProb <- function(pos=65) {
  CodonProb <- c(0.9910, 0.9750, 0.9793, 0.9691, 0.9318, 0.9268, 0.8389, 0.9636,
                 0.9622, 0.8830, 0.8633, 0.9223, 0.9179, 0.9598,      1, 0.9825,
                 0.9720, 0.8660, 0.8883, 0.9223, 0.9530, 0.8176, 0.7371, 0.9253,
                 0.5895, 0.5874, 0.4657, 0.8154, 0.9370, 0.7825, 0.8817, 0.9173,
                 0.9832, 0.9475, 0.9284, 0.9727, 0.9341, 0.9249, 0.8082, 0.9614,
                 0.8887, 0.8914, 0.8059, 0.9475, 0.9074, 0.9249, 0.9072, 0.9719,
                 0, 0.9460,      0, 0.9436, 0.9328, 0.9347, 0.8408, 0.9737,
                 0, 0.7542, 0.8870, 0.8534, 0.9722, 0.9748, 0.9819, 0.9703, 0)
  return(CodonProb[pos])
}
