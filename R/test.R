devtools::install_github("hesim-dev/hesim")
library(hesim)
library(data.table)



strategies <- data.table(strategy_id = 1:2,
                         strategy_name = c("Monotherapy", "Combination therapy"))
patients <- data.table(patient_id = 1)
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients)
print(hesim_dat)


### disease model
trans_mono <- matrix(c(1251, 350, 116, 17,
                       0, 731, 512, 15,
                       0, 0, 1312, 437,
                       0, 0, 0, 469),
                     ncol = 4, nrow = 4, byrow = TRUE)

params <- list(
  alpha_mono = trans_mono, 
  lrr_mean = log(.509), 
  lrr_lower <- log(.365),
  lrr_upper = log(.710)
)

rng_def <- define_rng({
  lrr_se <- (lrr_upper - lrr_lower)/(2 * qnorm(.975)) # Local object 
  # not returned
  list( # Parameters to return
    p_mono = dirichlet_rng(alpha_mono),
    rr_comb = lognormal_rng(lrr_mean, lrr_se)
  )
}, n = 1000)

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
    )
  )
  
}, times = c(2, Inf))


mod_def <- define_model(tparams_def = tparams_def, 
                        rng_def = rng_def, 
                        params = params)

####
print(eval_model(mod_def, data))

# utility
utility_tbl <- stateval_tbl(data.table(state_id = states$state_id,
                                    mean = c(0.9, 0.8, 0.6),
                                    se = c(0.1, 0.1, 0.1)),
                         dist = "beta",
                         hesim_data = hesim_dat)
###

# cost
cost_tbl <- stateval_tbl(data.table(state_id = states$state_id,
                                       mean = c(1701, 1774, 6948),
                                       se = c(170, 177, 694)),
                            dist = "gamma",
                            hesim_data = hesim_dat)
####

data <- expand(hesim_dat, 
                        by = c("strategies", "patients"))

print(data)




#model_inputs <- eval_model(mod_def, data)
#tparams <- tparams_transprobs(model_inputs)
#transmod <- CohortDtstmTrans$new(params = tparams) 
transmod <- create_CohortDtstmTrans(tparams_def, data, n=1000)

utilitymod <- create_StateVals(utility_tbl, n = 1000)
costmod <- create_StateVals(cost_tbl, n = 1000)

econmod <- CohortDtstm$new(trans_model = transmod,
                           utility_model = utilitymod,
                           cost_models = costmods)




