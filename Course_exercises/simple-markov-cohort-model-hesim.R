#install.packages("hesim")
#install.packages("data.table")
library(hesim)
library(data.table)





# Define the model
## Data: two treatment strategies and one representative patient.
strategies <- data.table(strategy_id = 1:2,
                         strategy_name = c("Monotherapy", "Combination therapy"))
patients <- data.table(patient_id = 1)
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients)
print(hesim_dat)

data <- expand(hesim_dat, by = c("strategies", "patients"))
head(data)

## Parameters
trans_mono <- matrix(c(1251, 350, 116, 17,
                       0, 731, 512, 15,
                       0, 0, 1312, 437,
                       0, 0, 0, 469),
                      ncol = 4, nrow = 4, byrow = TRUE)
colnames(trans_mono) <- rownames(trans_mono) <-  c("A", "B", "C", "D")
print(trans_mono)

params <- list(
  alpha_mono = trans_mono, 
  lrr_mean = log(.509), 
  lrr_lower <- log(.365),
  lrr_upper = log(.710),
  c_dmed_mean = c(A = 1701, B = 1774, C = 6948),
  c_cmed_mean = c(A = 1055, B = 1278, C = 2059),
  c_zido = 2278,
  c_lam = 2086.50,
  u = 1
)

## Random number generation
rng_def <- define_rng({
  lrr_se <- (lrr_upper - lrr_lower)/(2 * qnorm(.975)) # Local object 
                                                      # not returned
  list( # Parameters to return
    p_mono = dirichlet_rng(alpha_mono),
    rr_comb = lognormal_rng(lrr_mean, lrr_se),
    c_zido = c_zido,
    c_lam = c_lam,
    c_dmed = gamma_rng(mean = c_dmed_mean, sd = c_dmed_mean),
    c_cmed = gamma_rng(mean = c_cmed_mean, sd = c_cmed_mean),
    u = u
  )
}, n = 1000)

## Transformed parameters
tparams_def <- define_tparams({
  ## The treatment effect (relative risk) is transformed so that it varies by 
  ## strategies and only applies for the first 2 years (Monotherapy is 
  ## the reference strategy). Time intervals are closed on the left
  ## and open on the right so we use strict equality.  
  rr <- ifelse(strategy_name == "Monotherapy" | time > 2, 1, rr_comb)
  
  list(
    tpmatrix = tpmatrix(
      C, p_mono$A_B * rr, p_mono$A_C * rr, p_mono$A_D * rr,
      0, C, p_mono$B_C * rr, p_mono$B_D * rr,
      0, 0, C, p_mono$C_D * rr,
      0, 0, 0, 1
    ),
    utility = u,
    costs = list(
        drug = ifelse(strategy_name == "Monotherapy" | time > 2,
                    c_zido, c_zido + c_lam),
        community_medical = c_cmed,
        direct_medical = c_dmed
    )
  )
}, times = c(2, Inf))

## The model
mod_def <- define_model(tparams_def = tparams_def, 
                        rng_def = rng_def, 
                        params = params)







# Simulation
econmod <- create_CohortDtstm(mod_def, data)

## Health state probabilities
econmod$sim_stateprobs(n_cycles = 20)

### Plot
library(ggplot2)
theme_set(theme_bw())

stateprob_summary <- econmod$stateprobs_[, .(prob_mean = mean(prob),
                                              prob_lower = quantile(prob, .025),
                                              prob_upper = quantile(prob, .975)),
                                          by = c("strategy_id", "state_id", "t")]
stateprob_summary[, strategy_name := factor(strategy_id,
                                          labels = strategies$strategy_name)]
ggplot(stateprob_summary, aes(x = t, y = prob_mean)) +
  geom_line(aes(col = strategy_name)) +
  geom_ribbon(aes(x = t, ymin = prob_lower, ymax = prob_upper,
                  fill = strategy_name), alpha = .3) +
  facet_wrap(~factor(state_id, labels = LETTERS[1:4])) +
  xlab("Year") + ylab("Probability") +
  scale_fill_discrete("Strategy") + scale_color_discrete("Strategy")

## Costs and QALYs
econmod$sim_qalys(dr = 0, integrate_method = "riemann_right")
econmod$sim_costs(dr = 0.06, integrate_method = "riemann_right")






# Decision analysis
ce_sim <- econmod$summarize()
wtp <- seq(0, 25000, 500)
icea_pw_out <- icea_pw(ce_sim, comparator = 1, dr_qalys = 0, dr_costs = .06,
                       k = wtp)
icer_tbl(icea_pw_out)

ggplot(icea_pw_out$ceac, 
       aes(x = k, y = prob, 
           col = factor(strategy_id, labels = strategies$strategy_name[-1]))) +
  geom_line()  + xlab("Willingness to pay") +
  ylab("Probability cost-effective") +
  scale_x_continuous(breaks = seq(0, max(wtp), 5000), label = scales::dollar) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Strategy")

