#' The screening sensitivity
#'
#' \deqn{sensitivity = \frac{exp(\beta_0+\beta_1 diameter)}{1+exp(\beta_0+\beta_1 diameter)}}
#'
#'
#' @param d the tumor diameter
#' @param beta0 intercept parameter of the logistic function
#' @param beta1 slope parameter associated with tumor diameter
#' @param d0 the diameter threshold in mm from which a tumor is detectable by screening.
#'
#' @return the screening sensitivity
#' @export
Sensitivity <- function(d, beta0, beta1, d0){
  return(ifelse(d >= d0, 1 / (1 + exp(-(beta0 + beta1 * d))), 0))
}
