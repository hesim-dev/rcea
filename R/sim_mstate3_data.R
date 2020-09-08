#' Simulate example 3-state data
#'
#' Simulate example multi-state data for a 3-state model 
#' using Weibull distributions for each transition. Patients are 
#' right censored after 15 years of followup. 
#'
#' @param n The number of patients to simulate.
#' @param seed An integer to pass to `set.seed()`. If `NULL`, 
#' then no seed is set. 
#' @return Returns a `data.table` with multiple rows for patient. That is, 
#' for a given patient there are rows for each possible transition from a given
#' health state. The returned columns are as follows:
#' \describe{
#' \item{from}{The health state ID transitioned from.}
#' \item{to}{he health state ID transitioned to.}
#' \item{female}{`1` if female and `0` if male.}
#' \item{age}{The age of the patient in years.}
#' \item{patient_id}{The patient ID.}
#' \item{final}{An indicator equal to 1 if a patient is in their final health 
#'  state during the simulation and 0 otherwise.}
#' \item{time_start}{The time at the start of the interval.}
#' \item{time_stop}{The time and the end of the interval.}
#' \item{time}{The duration of the interval: `time_stop` - `time_start`.}
#' \item{status}{`1` if `time` is an observed transition and `0` if there was 
#' right censoring.}
#' \item{transition_id}{The transition ID representing a unique transition from
#' `from` to `to`.}
#' \item{strategy_id}{The treatment strategy ID.}
#' \item{strategy_name}{The name of the treatment strategy. Either `"New` 
#' or `"SOC`.}
#' }
#' 
#' @examples 
#' sim_mstate3_data(n = 3, seed = 101)
#' @import hesim 
#' @import data.table
#' @export
sim_mstate3_data <- function(n = 2000, seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  # Data  
  data <- data.table(
    intercept = 1,
    strategy_id = 1,
    patient_id = 1:n,
    female = rbinom(n, 1, .5),
    new = rbinom(n, 1, .5),
    age = rnorm(n, mean = 60, sd = 5.5)
  )
  attr(data, "id_vars") <- c("strategy_id", "patient_id")
  
  # Transition matrix
  tmat <- rbind(
    c(NA, 1,  2),
    c(NA,  NA, 3),
    c(NA, NA, NA)
  )
  trans_dt <- create_trans_dt(tmat)
  
  # Parameters for each transition
  get_scale <- function(shape, mean) {
    return(mean/(gamma(1 + 1/shape)))
  }
  
  matrixv <- function(v) {
    x <- matrix(v); colnames(x) <- "intercept"
    return(x)
  }
  
  params_wei <- function(shape, mean, 
                         beta_new = log(.6), 
                         beta_female = log(1.4)){
    log_shape <- matrixv(log(shape))
    scale = get_scale(shape, mean)
    beta_intercept <- log(scale) - beta_new
    scale_coefs <-  matrix(c(beta_intercept, beta_new, beta_female), 
                           ncol = 3)
    colnames(scale_coefs) <- c("intercept", "new", "female")
    params_surv(coefs = list(shape = log_shape,
                             scale = scale_coefs),
                dist = "weibull")
  }
    
  mstate_params <- params_surv_list(
    
    # 1. H -> S
    params_wei(shape = 2, mean = 1/.16),
    
    # 2. H -> D
    params_wei(shape = 3, mean = 10),
    
    # 3. S -> D
    params_wei(shape = 3.5, mean = 1/.12)
  )
  
  # Create multi-state model
  mstatemod <- create_IndivCtstmTrans(mstate_params, 
                                     input_data = data,
                                     trans_mat = tmat,
                                     clock = "reset",
                                     start_age = data$age)
  
  # Simulate data
  ## Observed "latent" transitions
  sim <- mstatemod$sim_disease(max_age = 100)
  sim[, c("sample", "grp_id", "strategy_id") := NULL]
  sim <- cbind(
    data[match(sim$patient_id, data$patient_id)][, patient_id := NULL],
    sim
  )
  sim[, ":=" (intercept = NULL, strategy_id = NULL, status = 1, added = 0)]
  
  ## Add all possible states for each transition
  ### Observed 1->2 add 1->3
  sim_13 <- sim[from == 1 & to == 2]
  sim_13[, ":=" (to = 3, status = 0, final = 0,  added = 1)]
  sim <- rbind(sim, sim_13)
  
  ### Observed 1->3 add 1->2
  sim_12 <- sim[from == 1 & to == 3 & added == 0]
  sim_12[, ":=" (to = 2, status = 0, final = 0, added = 1)]
  sim <- rbind(sim, sim_12)
  
  ### Sort and clean
  sim <- merge(sim, trans_dt, by = c("from", "to")) # Add transition ID
  setorderv(sim, c("patient_id", "from", "to"))
  sim[, added := NULL]
  
  ## Add right censoring
  sim[, status := ifelse(time_stop < 15, status, 0)]
  sim[, time_stop := pmin(time_stop, 15)]
  sim <- sim[time_start <= 15]

  ## Final data cleaning
  sim[, strategy_id := ifelse(new == 0, 1, 2)]
  sim[, strategy_name := factor(strategy_id, 
                                levels = c(1, 2),
                                labels = c("SOC", "New"))]
  sim[, new := NULL]
  
  # Return
  sim[, time := time_stop - time_start]
  return(sim[, ])
}