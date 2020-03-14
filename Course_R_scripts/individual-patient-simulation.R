# Required R-packages
# install.packages("devtools")
# devtools::install_github("hesim-dev/hesim")
# install.packages("data.table")
# install.packages("ggplot2")
library("hesim")
library("data.table")
library("ggplot2")

rm(list=ls())


# Define the population, treatment strategies, and model structure
tmat <- rbind(c(NA, 1, 2),
              c(3, NA, 4),
              c(NA, NA, NA))
colnames(tmat) <- rownames(tmat) <- c("Healthy", "Sick", "Dead")
print(tmat)

strategies <- data.table(strategy_id = c(1, 2))
n_patients <- 1000
patients <- data.table(patient_id = 1:n_patients,
                       age = rnorm(n_patients, mean = 45, sd = 7),
                       female = rbinom(n_patients, size = 1, prob = .51))
states <- data.table(state_id = c(1, 2),
                     state_name = c("Healthy", "Sick")) # Non-death health states
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients, 
                        states = states)



# Parameter estimation
## Multi-state model
library("flexsurv")
n_trans <- max(tmat, na.rm = TRUE) # Number of transitions
wei_fits_cr <- vector(length = n_trans, mode = "list") 
for (i in 1:length(wei_fits_cr)){
  wei_fits_cr[[i]] <- flexsurv::flexsurvreg(Surv(years, status) ~ factor(strategy_id), 
                                            data = mstate3_exdata$transitions, 
                                            subset = (trans == i) , 
                                            dist = "weibull") 
}
wei_fits_cr <- flexsurvreg_list(wei_fits_cr)

wei_fits_cf <- vector(length = n_trans, mode = "list") 
for (i in 1:length(wei_fits_cf)){
  wei_fits_cf[[i]] <- flexsurv::flexsurvreg(Surv(Tstart, Tstop, status) ~ factor(strategy_id), 
                                            data = mstate3_exdata$transitions, 
                                            subset = (trans == i) , 
                                            dist = "weibull") 
}
wei_fits_cf <- flexsurvreg_list(wei_fits_cf)

## Utility and costs
utility_tbl <- stateval_tbl(data.table(state_id = states$state_id,
                                       mean = mstate3_exdata$utility$mean,
                                       se = mstate3_exdata$utility$se),
                            dist = "beta",
                            hesim_data = hesim_dat)
head(utility_tbl)

drugcost_tbl <- stateval_tbl(data.table(strategy_id = strategies$strategy_id,
                                       est = mstate3_exdata$costs$drugs$costs),
                            dist = "fixed",
                            hesim_data = hesim_dat) 
medcost_tbl <- stateval_tbl(data.table(state_id = states$state_id,
                                       mean = mstate3_exdata$costs$medical$mean,
                                       se = mstate3_exdata$costs$medical$se),
                            dist = "gamma",
                            hesim_data = hesim_dat)  





# Simulation
## Constructing the economic model
n_samples <- 1000

### Disease model
transmod_data <- expand(hesim_dat, 
                        by = c("strategies", "patients"))
head(transmod_data)

transmod_cr <- create_IndivCtstmTrans(wei_fits_cr, transmod_data,
                                      trans_mat = tmat, n = n_samples,
                                      clock = "reset",
                                      start_age = patients$age)
transmod_cf <- create_IndivCtstmTrans(wei_fits_cf, transmod_data,
                                      trans_mat = tmat, n = n_samples,
                                      clock = "forward",
                                      start_age = patients$age)

### Utility
utilitymod <- create_StateVals(utility_tbl, n = n_samples)


### Costs
drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples)
medcostmod <- create_StateVals(medcost_tbl, n = n_samples)
costmods <- list(Drug = drugcostmod,
                 Medical = medcostmod)


### Combining the disease progression, cost, and utility models
econmod_cr <- IndivCtstm$new(trans_model = transmod_cr,
                             utility_model = utilitymod,
                             cost_models = costmods)
econmod_cf <- IndivCtstm$new(trans_model = transmod_cf,
                             utility_model = utilitymod,
                             cost_models = costmods)






## Simulating outcomes
### Disease progression

#### "Clock reset"
econmod_cr$sim_disease()
head(econmod_cr$disprog_)

#### "Clock forward"
econmod_cf$sim_disease()

econmod_cr$sim_stateprobs(t = seq(0, 20 , 1/12)) 

#### Short function to create state name variable to data.tabale
add_state_name <- function(x){
  x[, state_name := factor(state_id,
                           levels = 1:nrow(tmat),
                           labels = colnames(tmat))] 
}

#### Short function to create state probability "dataset" for plotting
summarize_stprobs <- function(stateprobs){
  x <- stateprobs[, .(prob_mean = mean(prob)),
                  by = c("strategy_id", "state_id", "t")]
  add_state_name(x)
}

#### Plot of state probabilities
stprobs_cr <- summarize_stprobs(econmod_cr$stateprobs_)
ggplot(stprobs_cr, aes(x = t, y = prob_mean, col = factor(strategy_id))) +
  geom_line() + facet_wrap(~state_name) + 
  xlab("Years") + ylab("Probability in health state") +
  scale_color_discrete(name = "Strategy") +
  theme(legend.position = "bottom") +
  theme_bw()

econmod_cf$sim_stateprobs(t = seq(0, 20 , 1/12)) 
stprobs_cf <- summarize_stprobs(econmod_cf$stateprobs_)

#### Compare "clock forward" and "clock reset" cases
stprobs <- rbind(data.table(stprobs_cf, clock = "Forward"),
                 data.table(stprobs_cr, clock = "Reset"))
ggplot(stprobs[strategy_id == 1], 
       aes(x = t, y = prob_mean, col = clock)) +
  geom_line() + facet_wrap(~state_name) + 
  xlab("Years") + ylab("Probability in health state") +
  scale_color_discrete(name = "Clock") +
  theme(legend.position = "bottom") +
  theme_bw()



### QALYs
econmod_cr$sim_qalys(dr = c(0,.03))
head(econmod_cr$qalys_)

qalys_summary <- econmod_cr$qalys_[, .(mean = mean(qalys)),
                                    by = c("strategy_id", "state_id", "dr")]
add_state_name(qalys_summary)
ggplot(qalys_summary[dr == .03],
       aes(x = factor(strategy_id), y = mean, fill = state_name)) + 
  geom_bar(stat = "identity") +
  scale_fill_discrete(name = "") +
  xlab("Strategy") + ylab("Mean QALYs") +
  theme_bw()



### Costs
econmod_cr$sim_costs(dr = 0.03)
head(econmod_cr$costs_)

library("scales")
costs_summary <- econmod_cr$costs_[dr == .03 , .(mean = mean(costs)),
                                   by = c("strategy_id", "category")]
ggplot(costs_summary,
       aes(x = factor(strategy_id), y = mean, fill = category)) + 
  geom_bar(stat = "identity") +
  scale_fill_discrete(name = "Category") +
  scale_y_continuous(label = scales::dollar_format()) +
  xlab("Strategy") + ylab("Mean costs") +
  theme_bw()




# Decision analysis
ce_sim <- econmod_cr$summarize()
icea_out <- icea(ce_sim, dr_qalys = .03, dr_costs = .03)
icea_pw_out <- icea_pw(ce_sim, comparator = 1, dr_qalys = .03, dr_costs = .03)

