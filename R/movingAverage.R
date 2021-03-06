#' @title Calculating a Moving Average
#'
#' @description
#' A way to calculate a moving average.
#'
#' @details
#' Suppose your data is a noisy sine wave with some missing values. A different way to handle missing data is to simply ignore it, and not include it in the average. The function defined here will do that.
#'
#' @param x the vector
#' @param n the number of samples
#' @param centered if FALSE, then average current sample and previous (n-1) samples
#'                 if TRUE, then average symmetrically in past and future. (If n is even, use one more sample from future.)
#' @return sum divided by count
#' @author  Cookbook, R. \email{winston@@stdout.org}
#'
#' @seealso \url{http://www.cookbook-r.com/Manipulating_data/Calculating_a_moving_average/}
#'
#' @examples
#' # Make same plots from before, with thicker lines
#' plot(x, y, type="l", col=grey(.5))
#' grid()
#' y_lag <- filter(y, rep(1/20, 20), sides=1)
#' lines(x, y_lag, col="red", lwd=4)         # Lagged average in red
#' y_sym <- filter(y, rep(1/21,21), sides=2)
#' lines(x, y_sym, col="blue", lwd=4)        # Symmetrical average in blue
#' # Calculate lagged moving average with new method and overplot with green
#' y_lag_na.rm <- movingAverage(y, 20)
#' lines(x, y_lag_na.rm, col="green", lwd=2)
# Calculate symmetrical moving average  with new method and overplot with green
#' y_sym_na.rm <- movingAverage(y, 21, TRUE)
#' lines(x, y_sym_na.rm, col="green", lwd=2)
#'
#' @export

movingAverage <- function(x, n=1, centered=FALSE) {

  if (centered) {
    before <- floor  ((n-1)/2)
    after  <- ceiling((n-1)/2)
  } else {
    before <- n-1
    after  <- 0
  }

  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))

  # Add the centered data
  new <- x
  # Add to count list wherever there isn't a
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new

  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])

    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new

    i <- i+1
  }

  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))

    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new

    i <- i+1
  }

  # return sum divided by count
  s/count
}
