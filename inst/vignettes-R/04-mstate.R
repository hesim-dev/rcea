## ---- Overview ---------------------------------------------------------------
## @knitr R-setup
library("rcea")
library("hesim")
library("data.table")
library("ggplot2")
library("flexsurv")

set.seed(101) # Make random number generation reproducible

## ---- Model setup ------------------------------------------------------------
## @knitr tmat
tmat <- rbind(
  c(NA, 1, 2),
  c(NA, NA, 3),
  c(NA, NA, NA)
)
colnames(tmat) <- rownames(tmat) <- c("Stable", "Progression", "Dead")
print(tmat)

## @knitr hesim_data
n_patients <- 1000
patients <- data.table(
  patient_id = 1:n_patients,
  age = rnorm(n_patients, mean = 45, sd = 7),
  female = rbinom(n_patients, size = 1, prob = .51)
)

states <- data.table(
  state_id = c(1, 2),
  state_name = c("Stable", "Progression") # Non-death health states
)
n_states <- nrow(states)

strategies <- data.frame(
  strategy_id = 1:2,
  strategy_name = c("SOC", "New")
)
n_strategies <- nrow(strategies)

hesim_dat <- hesim_data(
  strategies = strategies,
  patients = patients,
  states = states
)
print(hesim_dat)

## ---- Parameter estimation ---------------------------------------------------
## @knitr sim_mstate3_data
data <- rcea::sim_mstate3_data(n = 2000)
data[patient_id %in% c(1, 2)]

## @knitr fit-mstate
n_trans <- max(tmat, na.rm = TRUE) # Number of transitions
wei_fits <- vector(length = n_trans, mode = "list")
for (i in 1:length(wei_fits)){
  wei_fits[[i]] <- flexsurvreg(
    Surv(time, status) ~ strategy_name + female,
    data = data,
    subset = (transition_id == i) ,
    dist = "weibullPH")
}
wei_fits <- flexsurvreg_list(wei_fits)

## @knitr utility_tbl
utility_tbl <- stateval_tbl(
  data.table(state_id = states$state_id,
             mean = c(.8, .6),
             se = c(0.02, .05)
  ),
  dist = "beta",
  hesim_data = hesim_dat)
print(utility_tbl)

## @knitr medcost_tbl
medcost_tbl <- stateval_tbl(
  data.table(state_id = states$state_id,
             mean = c(2000, 9500),
             se = c(2000, 9500)
  ),
  dist = "gamma",
  hesim_data = hesim_dat)
print(medcost_tbl)

## @knitr drugcost_tbl
n_times <- 2
drugcost_tbl <- stateval_tbl(
  data.table(
    strategy_id = rep(strategies$strategy_id, each = n_states * n_times),
    state_id = rep(rep(states$state_id, each = n_strategies), n_times),
    time_start = rep(c(0, 3/12), n_states * n_strategies),
    est = c(rep(2000, 4), # Costs are always the same with SOC
            12000, 12000, 12000, 10000 # Costs with new drop after 3 months in progression state
    )  
  ),
  dist = "fixed",
  hesim_data = hesim_dat)
print(drugcost_tbl)

## ---- Simulation -------------------------------------------------------------
## ---- Construct model --- ##
## @knitr psa-iterations
n_samples <- 500

## @knitr transition-model-data
transmod_data <- expand(hesim_dat,
                        by = c("strategies", "patients"))
head(transmod_data)

## @knitr transition-model
transmod <- create_IndivCtstmTrans(wei_fits, transmod_data,
                                   trans_mat = tmat, n = n_samples,
                                   clock = "reset",
                                   start_age = patients$age)

## @knitr utility-cost-models
# Utility
utilitymod <- create_StateVals(utility_tbl, n = n_samples)
# Costs
drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples,
                                time_reset = TRUE)
medcostmod <- create_StateVals(medcost_tbl, n = n_samples)
costmods <- list(Drug = drugcostmod,
                 Medical = medcostmod)

## @knitr economic-model
econmod <- IndivCtstm$new(trans_model = transmod,
                          utility_model = utilitymod,
                          cost_models = costmods)

## ---- Simulate outcomes --- ##
## @knitr sim_disease
econmod$sim_disease(max_age = 100)
head(econmod$disprog_)

## @knitr sim_stateprobs
econmod$sim_stateprobs(t = seq(0, 30 , 1/12))

## @knitr sim_qalys
econmod$sim_qalys(dr = c(0,.03))
head(econmod$qalys_)

## @knitr sim_costs
econmod$sim_costs(dr = 0.03)
head(econmod$costs_)

## ---- Cost-effectiveness analysis --------------------------------------------
## @knitr icer
ce_sim <- econmod$summarize()
cea_pw_out <- cea_pw(ce_sim, comparator = 1, 
                       dr_qalys = .03, dr_costs = .03,
                       k = seq(0, 25000, 500))
icer_tbl(cea_pw_out, colnames = strategies$strategy_name)

## @knitr save-ce_sim
saveRDS(ce_sim, "ce_sim.rds")