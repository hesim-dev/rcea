#install.packages("hesim")
#install.packages("data.table")
library("data.table")
library("hesim")
library("ggplot2")
library("scales")


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

## PSA output 
print(icea_pw_out$delta)

## cost-effectiveness plane
ylim <- max(icea_pw_out$delta[, ic]) * 1.1
xlim <- ceiling(max(icea_pw_out$delta[, ie]) * 1.1)
ggplot2::ggplot(icea_pw_out$delta, aes(x = ie, y = ic)) + 
  geom_jitter(size = .5)  + 
  xlab("Incremental QALYs") + ylab("Incremental cost") +
  scale_y_continuous(limits = c(-ylim, ylim)) +
  scale_x_continuous(limits = c(-xlim, xlim), breaks = seq(-6, 6, 1)) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Strategy") +
  geom_abline(slope = 150000, linetype = "dashed") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

## cost-effectiveness acceptability curves
theme_set(theme_minimal())
ggplot2::ggplot(icea_out$mce, aes(x = k, y = prob, col = factor(strategy_id))) +
  geom_line()  + 
  xlab("Willingness to pay") +
  ylab("Probability most cost-effective") +
  scale_x_continuous(breaks = seq(0, 300000, 100000), label = scales::dollar) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Strategy")

## Cost-effectiveness acceptability frontier
ggplot2::ggplot(icea_out$mce[best == 1], aes(x = k, y = prob, col = factor(strategy_id))) +
  geom_line() + 
  xlab("Willingness to pay") +
  ylab("Probability most cost-effective") +
  scale_x_continuous(breaks = seq(0, 300000, 100000), label = scales::dollar) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Strategy")

## Value of perfect information
ggplot2::ggplot(icea_out$evpi, aes(x = k, y = evpi)) +
  geom_line() + 
  xlab("Willingness to pay") +
  ylab("Expected value of perfect information") +
  scale_x_continuous(breaks = seq(0, 300000, 100000), label = scales::dollar) +
  scale_y_continuous(label = scales::dollar) +
  theme(legend.position = "bottom") 
