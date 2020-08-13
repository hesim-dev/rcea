#' Simulate a Markov chain
#'
#' This is a simple function that will simulate a time-homogeneous Markov chain;
#' that is, a Markov chain with constant transition probabilities.
#'
#' @param x0 The state vector at model cycle 0 (i.e., the initial state vector).
#' @param p The transition probability matrix.
#' @param n_cycles The number of model cycles. Default is 85 as in the
#' example from the tutorials. 
#' 
#' @examples 
#' x0 <- c(1, 0, 0, 0)
#' p <- matrix(
#'   c(0.848, 0.15,  0,     0.002,
#'     0.5,   0.395, 0.105, 0.006,
#'     0,     0,     0.98,  0,02,
#'     0,     0,     0,     1),
#'   byrow = TRUE,
#'   nrow = 4, ncol = 4
#' )
#' sim_markov_chain(x0 = x0, p = p, n_cycles = 2)
#' @export
sim_markov_chain <- function(x0, p, n_cycles = 85){
  x <- matrix(NA, ncol = length(x0), nrow = n_cycles) # Initialize Markov trace
  x <- rbind(x0, x) # Markov trace at cycle 0 is initial state vector
  colnames(x) <- colnames(p) # Columns are the health states
  rownames(x) <- 0:n_cycles # Rows are the model cycles
  for (t in 1:n_cycles){ # Simulating state vectors at each cycle with for loop
    x[t + 1, ] <- x[t, ] %*% p
  }
  return(x)
}