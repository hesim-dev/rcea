#' Simulate progression-free and overall survival
#'
#' Simulate a  multi-state dataset for a 3-state model with  
#' `sim_mstate3_date()` and compute progression-free survival (PFS)
#' and overall survival (OS).
#'
#' @inheritParams sim_mstate3_data
#' 
#' @return Returns a `data.table` with one rows for each patient. 
#' The returned columns are as follows:
#' \describe{
#' \item{patient_id}{The patient ID.}
#' \item{strategy_id}{The treatment strategy ID.}
#' \item{strategy_name}{The name of the treatment strategy. Either `"New` 
#' or `"SOC`.}
#' \item{female}{`1` if female and `0` if male.}
#' \item{age}{The age of the patient in years.}
#' \item{pfs_time}{Time to progression, death, or right censoring.}
#' \item{pfs_status}{`1` if a progression or death event occured and `0` 
#' if there was right censoring.}
#' \item{os_time}{Time to death or right censoring.}
#' \item{os_status}{`1` if a patient died and `0` if there was right censoring.}
#' }
#' 
#' @examples 
#' sim_pfs_os_data(n = 3, seed = 101)
#' @import hesim 
#' @import data.table
#' @export
sim_pfs_os_data <- function(n = 2000, seed = NULL){
  sim <- sim_mstate3_data(n = n, seed = seed)
  
  # Reshape wider with one row per patient
  sim <- dcast(sim, 
               patient_id + strategy_id + strategy_name + 
                 female + age ~ transition_id,
               value.var = c("status", "time_stop"))
  
  # Compute PFS endpoint
  sim[, pfs_time := pmin(time_stop_1, time_stop_2)]
  sim[, pfs_status := pmax(status_1, status_2)]
  
  # Compute OS endpoint
  sim[, os_status := fcase(
    status_2 == 1 | status_3 == 1, 1L,
    status_2 == 0 & status_3 == 0, 0L,
    status_2 == 0 & is.na(status_3), 0L
    )
  ]
  
  sim[, os_time := fcase(
    status_2 == 1, time_stop_2, # Observed death 1 -> 3
    status_2 == 0 & is.na(status_3), time_stop_2, # Right censored 1 -> 3
    status_2 == 0 & status_3 == 1, time_stop_3, # Observed death 2 -> 3
    status_2 == 0 & status_3 == 0,  time_stop_3 # Right censored 2 -> 3
    )
  ]
  
  # Clean and return
  sim[, c("status_1", "status_2", "status_3", 
          "time_stop_1", "time_stop_2", "time_stop_3")
      := NULL]
  return(sim[, ])
  
}