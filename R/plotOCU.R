#' @title Plot of Optimality Score for Codons in a Single Protein
#'
#' @description
#' \code{plotOCU} Plots the Optimality Score for Codon in a Protein based on tRNA-Gene Abundancy in a Custom Genome (eg \emph{Eragostis tef})
#'
#' @details
#' Here a plot for the optimality of the used codons of a protein is plotted for its position in the protein sequence.
#' Two mooving averages are overlayed for better visuality.
#' In development stage!
#'
#' @param entry Index number in SeqFastadna list myaa
#' @param isfreq frequencies are used instead of boolean of optimaity
#'
#' @return void
#'
#' @author Siegrist, F. and Cannarozzi, G. M. \email{gina@@cannarozzi.com}
#' @seealso \code{\link{movingAverage}} \code{\link{statanacoseq}}
#' @keywords CodonBias
#' @examples
#' for (n in 1:length(mylist()[[2]]) try(plotOCU(n, TRUE))
#'
#' @export
#'

plotOCU <- function(entry=3, isfreq=FALSE) {
  areopt <- areopts(entry, !isfreq)
  x <- 1:length(areopt)
  y <- as.numeric(areopt)
  if (!exists("myaa")) myaa <- mylist()[[2]]
  title <- getName(myaa[[entry]])
  # Make same plots from before, with thicker lines
  plot(x, y, type="p", col=grey(.5), main=title, xlab='aa position on gene', ylab='moving average of "is optimal codon"')
  grid()
  y_lag <- filter(y, rep(1/20, 20), sides=1)
  lines(x, y_lag, col="red", lwd=4)         # Lagged average in red
  y_sym <- filter(y, rep(1/21,21), sides=2)
  lines(x, y_sym, col="blue", lwd=4)        # Symmetrical average in blue

  # Calculate lagged moving average with new method and overplot with green
  y_lag_na.rm <- movingAverage(y, 20)
  lines(x, y_lag_na.rm, col="green", lwd=2)

  # Calculate symmetrical moving average  with new method and overplot with green
  y_sym_na.rm <- movingAverage(y, 21, TRUE)
  lines(x, y_sym_na.rm, col="green", lwd=2)
}
