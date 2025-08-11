#' Log-linear transformation function for Clock 3
#'
#' This function performs the log-linear transformation used in Universal Clock 3.
#' 
#' @param age1 Numeric vector of ages to transform
#' @param m1 Numeric parameter m1 for transformation
#' @param m2 Numeric parameter m2 for transformation (default: same as m1)
#' @param c1 Numeric parameter c1 for transformation (default: 1)
#' @return Numeric vector of transformed values
#' @export
#' @examples
#' ages <- c(0.1, 0.5, 1.0, 2.0)
#' transformed <- F1_logli(ages, m1 = 0.5)
F1_logli <- function(age1, m1, m2 = m1, c1 = 1) {
  ifelse(age1 >= m1, 
         (age1 - m1) / m2, 
         c1 * log((age1 - m1) / m2 / c1 + 1))
}

#' Reverse transformation function for Clock 3
#'
#' This function performs the reverse transformation for Universal Clock 3 predictions.
#' 
#' @param y.pred Numeric vector of predicted values to reverse transform
#' @param m1 Numeric parameter m1 for transformation
#' @param m2 Numeric parameter m2 for transformation (default: same as m1)
#' @param c1 Numeric parameter c1 for transformation (default: 1)
#' @return Numeric vector of original scale values
#' @export
#' @examples
#' predictions <- c(-1, 0, 1, 2)
#' original <- F2_revtrsf(predictions, m1 = 0.5)
F2_revtrsf <- function(y.pred, m1, m2 = m1, c1 = 1) {
  ifelse(y.pred < 0, 
         (exp(y.pred / c1) - 1) * m2 * c1 + m1, 
         y.pred * m2 + m1)
}

#' Calculate parameters for log-linear transformation
#'
#' This function calculates the necessary parameters for the log-linear transformation
#' used in Universal Clock 3.
#' 
#' @param dat1 Data frame containing species information with columns: maxAgeCaesar, 
#'   GestationTimeInYears, averagedMaturity.yrs, Age
#' @param b1 Numeric parameter b1 (default: 1)
#' @param max_tage Numeric maximum transformed age (default: 4)
#' @param c1 Numeric parameter c1 (default: 5)
#' @param c2 Numeric parameter c2 (default: 0.38)
#' @param c0 Numeric parameter c0 (default: 0)
#' @return Data frame with additional columns for transformation parameters
#' @export
#' @importFrom dplyr mutate
#' @examples
#' \dontrun{
#' # Requires species data with appropriate columns
#' species_data <- data.frame(
#'   maxAgeCaesar = 4, 
#'   GestationTimeInYears = 0.055,
#'   averagedMaturity.yrs = 0.4,
#'   Age = c(0.1, 0.5, 1.0, 2.0)
#' )
#' transformed_data <- F3_loglifn(species_data)
#' }
F3_loglifn <- function(dat1, b1 = 1, max_tage = 4, c1 = 5, c2 = 0.38, c0 = 0) {
  n <- nrow(dat1)
  
  age1 <- (dat1$maxAgeCaesar + dat1$GestationTimeInYears) / 
          (dat1$averagedMaturity.yrs + dat1$GestationTimeInYears)
  
  a1 <- age1 / (1 + max_tage)
  dat1$a1_Logli <- a1 # x/m1 in manuscript
  
  a2 <- (dat1$GestationTimeInYears + c0) / (dat1$averagedMaturity.yrs)
  dat1$a_Logli <- a_Logli <- c1 * a2^c2
  # m=5*(G/ASM)^0.38 from regression analysis/formula(7)
  
  x <- dat1$Age + dat1$GestationTimeInYears
  t2 <- dat1$averagedMaturity.yrs * b1 + dat1$GestationTimeInYears
  x2 <- x / t2 #### log(x/t2)
  y <- F1_logli(x2, a_Logli, a_Logli)
  
  dat1$LogliAge <- y
  return(dat1)
}
