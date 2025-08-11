#' Age transformation function
#'
#' This function performs age transformation used in epigenetic clock calculations.
#' 
#' @param x Numeric vector of ages to transform
#' @param offset Numeric offset value (default: 0.06)
#' @param adult.age Numeric adult age threshold (default: 1.2)
#' @return Numeric vector of transformed ages
#' @export
#' @examples
#' ages <- c(0.5, 1.0, 2.0, 5.0)
#' transformed_ages <- trafo(ages)
trafo <- function(x, offset = 0.06, adult.age = 1.2) {
  y <- ifelse(x <= adult.age, 
              log(x + offset),
              x / (adult.age + offset) + log(adult.age + offset) - adult.age / (adult.age + offset))
  return(y)
}

#' Inverse age transformation function
#'
#' This function performs the inverse of the age transformation used in epigenetic clock calculations.
#' 
#' @param x Numeric vector of transformed ages to reverse
#' @param offset Numeric offset value (default: 0.06)
#' @param adult.age Numeric adult age threshold (default: 1.2)
#' @return Numeric vector of original ages
#' @export
#' @examples
#' transformed_ages <- trafo(c(0.5, 1.0, 2.0, 5.0))
#' original_ages <- anti.trafo(transformed_ages)
anti.trafo <- function(x, offset = 0.06, adult.age = 1.2) {
  ifelse(x <= log(adult.age + offset), 
         exp(x) - offset, 
         (adult.age + offset) * x - log(adult.age + offset) * (adult.age + offset) + adult.age)
}

#' Universal clock inverse transformation function
#'
#' Inverse transformation function used in universal clocks (F2 type).
#' 
#' @param y Numeric vector of predicted values to transform
#' @param y.maxAge Numeric maximum age for the species
#' @param y.gestation Numeric gestation time in years
#' @param const Numeric constant (default: 1)
#' @return Numeric vector of transformed ages
#' @export
#' @examples
#' # Example usage with mouse parameters
#' predicted_values <- c(-2, -1, 0, 1, 2)
#' ages <- F2_antitrans(predicted_values, y.maxAge = 4, y.gestation = 0.055)
F2_antitrans <- function(y, y.maxAge, y.gestation, const = 1) {
  x0 <- const * exp(-exp(-1 * y))
  x1 <- x0 * (y.maxAge + y.gestation)
  x <- x1 - y.gestation
  return(x)
}
