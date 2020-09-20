## ---- Overview ---------------------------------------------------------------
## @knitr R-packages
library("hesim")

## ---- Model setup ------------------------------------------------------------
## @knitr hesim-data
strategies <- data.frame(
  strategy_id = 1:2,
  strategy_name = c("SOC", "New")
)
patients <- data.frame(
  patient_id = 1,
  age = 25
)
hesim_dat <- hesim_data(
  strategies = strategies,
  patients = patients
)
print(hesim_dat)

## ---- Model parameters -------------------------------------------------------
## @knitr transitions
transitions_soc <- matrix(
  c(848, 150, 0,   2,
    450, 355, 95,  5,
    0,   0,   784, 16,
    0,   0,   0,   23),
  nrow = 4, byrow = TRUE)
state_names <- c("H", "S1", "S2", "D")
colnames(transitions_soc) <- rownames(transitions_soc) <- tolower(state_names)

## @knitr all-parameters
params <- list(
  alpha_soc = transitions_soc,
  lrr_mean = log(.8),
  lrr_lower = log(.71),
  lrr_upper = log(.9),
  c_medical = c(H = 2000, S1 = 4000, S2 = 15000),
  c_soc = 2000,
  c_new = 12000,
  u_mean = c(H = 1, S1 = .075, S2 = 0.5),
  u_se = c(H = 0, S1 = 0.03, S2 = 0.05)
)

## @knitr define_rng
rng_def <- define_rng({
  lrr_se <- (lrr_upper - lrr_lower)/(2 * qnorm(.975)) # Local object 
  # not returned
  list( # Parameters to return
    p_soc = dirichlet_rng(alpha_soc),
    rr_new = lognormal_rng(lrr_mean, lrr_se),
    c_medical = gamma_rng(mean = c_medical, sd = c_medical),
    c_soc = c_soc,
    c_new = c_new,
    u = beta_rng(mean = u_mean, sd = u_se)
  )
}, n = 1000)

## @knitr expanded-data
input_data <- hesim::expand(hesim_dat, by = c("strategies", "patients"))
head(input_data)

## @knitr define_tparams
tparams_def <- define_tparams({
  ## The treatment effect (relative risk) is transformed so that it varies by 
  ## strategies (SOC is the reference strategy)
  rr <- ifelse(strategy_name == "SOC", 1, rr_new)
  
  list(
    tpmatrix = tpmatrix(
      C,          p_soc$h_s1 * rr, p_soc$h_s2 * rr,  p_soc$h_d * rr,
      p_soc$s1_h, C,               p_soc$s1_s2 * rr, p_soc$s1_d * rr,
      p_soc$s2_h, p_soc$s2_s1,     C,                p_soc$s2_d * rr,
      0,          0,               0,                1
    ),
    utility = u,
    costs = list(
      treatment = ifelse(strategy_name == "SOC", c_soc, c_new),
      medical = c_medical
    )
  )
})

## ---- Simulation -------------------------------------------------------------
## ---- Construct model --- ##
## @knitr define_model
mod_def <- define_model(tparams_def = tparams_def, 
                        rng_def = rng_def, 
                        params = params)

## @knitr initialize-model
econmod <- create_CohortDtstm(mod_def, input_data)

## ---- Simulate outcomes --- ##
## @knitr sim_stateprobs
econmod$sim_stateprobs(n_cycles = 85)

## @knitr sim_qalys
econmod$sim_qalys(dr = 0.03, integrate_method = "riemann_left")
head(econmod$qalys_)

## @knitr sim_costs
econmod$sim_costs(dr = 0.03, integrate_method = "riemann_left")

## ---- Cost-effectiveness analysis --------------------------------------------
## @knitr cea
ce_sim <- econmod$summarize()
cea_pw_out <- cea_pw(ce_sim, comparator = 1, 
                     dr_qalys = 0.03, dr_costs = 0.03,
                     k = seq(0, 25000, 500))

## @knitr icer
icer_tbl(cea_pw_out, colnames = strategies$strategy_name) 

## ---- Exercises --------------------------------------------------------------