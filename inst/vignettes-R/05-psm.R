## ---- Overview ---------------------------------------------------------------
## @knitr R-setup
library("rcea")
library("hesim")
library("data.table")
library("ggplot2")
library("flexsurv")
theme_set(theme_bw())

set.seed(101) # Make PSA reproducible

## ---- Model setup ------------------------------------------------------------
## @knitr hesim_data
patients <- data.table( 
  patient_id = 1:4,
  patient_wt = rep(1/4, 4), # Each patient has same weight
  age = c(45, 45, 65, 65),
  female = c(0, 1, 0, 1)
)

states <- data.table(
  state_id = c(1, 2),
  state_name = c("Stable", "Progression") # Non-death health states
)

strategies <- data.frame(
  strategy_id = 1:2,
  strategy_name = c("SOC", "New")
)

hesim_dat <- hesim_data(
  patients = patients,
  strategies = strategies,
  states = states
)
print(hesim_dat)

## @knitr labels
labs <- get_labels(hesim_dat)

## ---- Parameter estimation ---------------------------------------------------
## @knitr pfs_os_data
onc3 <- hesim::onc3[strategy_name != "New 1"]
onc3[, strategy_name := droplevels(strategy_name)]
levels(onc3$strategy_name) <- c("SOC", "New")
surv_est_data <- as_pfs_os(onc3, patient_vars = c("patient_id", "female", "age",
                                                  "strategy_name"))
surv_est_data[patient_id %in% c(1, 2)]

## @knitr fit-survival-models
fit_pfs_wei <- flexsurvreg(
  Surv(pfs_time, pfs_status) ~ strategy_name + female,
  data = surv_est_data,
  dist = "weibull")

fit_os_wei <- flexsurvreg(
  Surv(os_time, os_status) ~ strategy_name + female,
  data = surv_est_data,
  dist = "weibull")

psfit_wei <- flexsurvreg_list(fit_pfs_wei, fit_os_wei)

## @knitr utility-medcost-tables
# Utility
utility_tbl <- stateval_tbl(
  data.table(state_id = states$state_id,
             mean = c(.8, S = .6),
             se = c(0.02, .05)
  ),
  dist = "beta"
)
print(utility_tbl)

# Medical costs
medcost_tbl <- stateval_tbl(
  data.table(state_id = states$state_id,
             mean = c(2000, 9500),
             se = c(2000, 9500)
  ),
  dist = "gamma"
)
print(medcost_tbl)

## @knitr drugcost_tbl
drugcost_tbl <- stateval_tbl(
  data.table(strategy_id = strategies$strategy_id,
             est = c(2000, 12000)),
  dist = "fixed"
)
print(drugcost_tbl)

## ---- Simulation -------------------------------------------------------------
## ---- Construct model --- ##
## @knitr psa-iterations
n_samples <- 100

## @knitr survival-models-data
# Input data for survival models
survmods_data <- expand(hesim_dat, by = c("strategies", "patients"))
print(survmods_data)

## @knitr survival-models
survmods <- create_PsmCurves(psfit_wei, 
                             input_data = survmods_data, 
                             n = n_samples,
                             uncertainty = "bootstrap", 
                             est_data = surv_est_data)

## @knitr utility-cost-models 
# Utility
utilitymod <- create_StateVals(utility_tbl, n = n_samples,
                               hesim_data = hesim_dat)

# Costs
drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples,
                                hesim_data = hesim_dat)
medcostmod <- create_StateVals(medcost_tbl, n = n_samples,
                               hesim_data = hesim_dat)
costmods <- list(Drug = drugcostmod,
                 Medical = medcostmod)

## @knitr economic-model
econmod <- Psm$new(survival_models = survmods,
                   utility_model = utilitymod,
                   cost_models = costmods)

## ---- Simulate outcomes --- ##
## @knitr sim_survival
times <- seq(0, 50, by = .1)
econmod$sim_survival(t = times)
autoplot(econmod$survival_, ci = TRUE,
         labels = labs)

## @knitr sim_stateprobs
econmod$sim_stateprobs()
econmod$stateprobs_[sample == 1 & state_id == 2 & t == 12]

## @knitr sim-costs-qalys
econmod$sim_costs(dr = .03)
econmod$sim_qalys(dr = .03)

## ---- Cost-effectiveness analysis --------------------------------------------
## @knitr icer
ce_sim <- econmod$summarize()
cea_pw_out <- cea_pw(ce_sim, comparator = 1, dr_qalys = .03, dr_costs = .03)
format(icer(cea_pw_out, labels = labs))