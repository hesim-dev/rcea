#' Present value
#'
#' Compute the present value of a numeric quantiy `z` given a discount rate
#' `dr` and times `t`.
#'
#' @param z A numeric quantity.
#' @param dr The discount rate.
#' @param t A vector of times to compute the present value.
#' 
#' @examples 
#' pv(1000, dr = .03, t = 0:4)
#' @export
pv <- function(z, dr, t) {
  z/(1 + dr)^t
}