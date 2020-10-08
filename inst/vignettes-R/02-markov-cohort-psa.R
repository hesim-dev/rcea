## ---- Overview ---------------------------------------------------------------
## @knitr R-packages
library("rcea")
library("hesim")
library("data.table")
library("magrittr")
library("ggplot2")

## ---- Model parameters -------------------------------------------------------
## @knitr tpmatrix
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
  c_medical = c(H = 2000, S1 = 4000, S2 = 15000, D = 0),
  c_soc = 2000,
  c_new = 12000,
  u_mean = c(H = 1, S1 = .75, S2 = 0.5, D = 0),
  u_se = c(H = 0, S1 = 0.03, S2 = 0.05, D = 0.0)
)

## ---- Simulation -------------------------------------------------------------
## @knitr sample-parameters
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
params_rng <- eval_rng(rng_def, params = params)
attr(params_rng, "n") <- rng_def$n
names(params_rng)
head(as.matrix(params_rng$p_soc))

## @knitr input-data
data <- data.frame(
  strategy = c("New", "SOC")
)
head(data)

## @knitr sim_stateprobs-fun
sim_stateprobs <- function(p0, rr, strategy, n_cycles){
  rr <- ifelse(strategy == "New", rr, 1)
  p <- tpmatrix(
    C,       p0$h_s1 * rr,  p0$h_s2 * rr,  p0$h_d * rr,
    p0$s1_h, C,             p0$s1_s2 * rr, p0$s1_d * rr,
    p0$s2_h, p0$s2_s1,      C,             p0$s2_d * rr,
    0,       0,             0,             1
  )
  x <- sim_markov_chain(x0 = c(1, 0, 0, 0),
                        p = matrix(as.matrix(p), ncol = 4, byrow = TRUE),
                        n_cycles = n_cycles)
  return(x)
}

## @knitr compute_qalys-fun
compute_qalys <- function(x, utility, dr = .03){
  n_cycles <- nrow(x) - 1
  pv(x %*% utility, dr, 0:n_cycles)
}

## @knitr compute_costs-fun
compute_costs <- function(x, costs_medical, costs_treat, dr = .03){
  n_cycles <- nrow(x) - 1
  costs_treat <- c(rep(costs_treat, 3), 0)
  costs <- cbind(
    pv(x %*% costs_medical, dr, 0:n_cycles),
    pv(x %*% costs_treat, dr, 0:n_cycles)
  )
  colnames(costs) <- c("dcost_med", "dcost_treat")
  return(costs)
}

## @knitr sim_model-fun
sim_model <- function(params_rng, data, n_cycles = 85, 
                      dr_qalys = .03, dr_costs = .03){
  # Initialize array of matrices
  n_samples <- attr(params_rng, "n")
  n_strategies <- nrow(data)
  out <- array(NA, dim = c(n_cycles + 1, 7, n_samples * n_strategies))
  dimnames(out) <- list(NULL, 
                        c("H", "S1", "S2", "D",
                          "dqalys", "dcosts_med", "dcosts_treat"), 
                        NULL)
  
  # Run the simulation
  i <- 1
  for (s in 1:n_samples){ # Start PSA loop
    for (k in 1:n_strategies) { # Start treatment strategy loop
      x <- sim_stateprobs(p0 = params_rng$p_soc[s, ],
                          rr = params_rng$rr_new[s],
                          strategy = data$strategy[k],
                          n_cycles = n_cycles)
      dqalys <- compute_qalys(x, utility = unlist(params_rng$u[s]), 
                              dr = dr_qalys)
      dcosts <- compute_costs(x, 
                              costs_medical = unlist(params_rng$c_medical[s]), 
                              costs_treat = ifelse(data$strategy[k] == "SOC", 
                                                   params_rng$c_soc,
                                                   params_rng$c_new),
                              dr = dr_costs)
      out[, , i] <- cbind(x, dqalys, dcosts)
      i <- i + 1
    } # End treatment strategy loop
  } # End PSA loop
  
  # Store metadata and return
  attr(out, "n_samples") <- n_samples
  attr(out, "strategies") <- data$strategy
  return(out)
}

## @knitr run_sim
sim_out <- sim_model(params_rng, data = data)
head(sim_out[, , 1])

## @knitr array-to-data-table
# rbind an array of matrices into a single matrix
rbind_array <- function(x){
  n_rows <- dim(x)[3] * dim(x)[1]
  x_mat <- matrix(c(aperm(x, perm = c(2, 1, 3))),
                  nrow = n_rows, byrow = TRUE)
  colnames(x_mat) <- dimnames(x)[[2]]
  return(x_mat)
}

# Convert the array into a long dataframe with ID columns
array_to_dt <- function(x){
  id_df <- expand.grid(cycle = 0:(dim(x)[1] - 1),
                       strategy = attr(x, "strategies"),
                       sample = 1:attr(x, "n_samples"))
  x_mat <- rbind_array(x)
  return(as.data.table(cbind(id_df, x_mat)))
}

sim_out <- array_to_dt(sim_out)
head(sim_out)

## ---- Cost-effectiveness analysis --------------------------------------------
## @knitr ce_output
ce_out <- sim_out[cycle != 0, 
                  .(dqalys = sum(dqalys),
                    dcosts = sum(dcosts_med) + sum(dcosts_treat)), 
                  by = c("sample", "strategy")]
ce_out

## @knitr ce_output_wider
ce_out_wider <- dcast(ce_out, sample ~ strategy, 
                      value.var = c("dqalys", "dcosts"))
ce_out_wider

## @knitr icer
ce_out_wider[, idcosts := dcosts_New - dcosts_SOC]
ce_out_wider[, idqalys := dqalys_New - dqalys_SOC]
ce_out_wider[, .(icer = mean(idcosts)/mean(idqalys))]

## @knitr ceplane
format_dollar <- function(x) {
  paste0("$", formatC(x, format = "d", big.mark = ","))
}

ylim <- max(ce_out_wider$idcosts) * 1.2
xlim <- max(ce_out_wider$idqalys) * 1.2
ggplot(ce_out_wider, 
       aes(x = idqalys, y = idcosts)) + 
  geom_jitter(size = .5)  + 
  xlab("Incremental QALYs") + 
  ylab("Incremental cost") +
  scale_y_continuous(limits = c(-ylim, ylim),
                     labels = format_dollar) +
  scale_x_continuous(limits = c(-xlim, xlim), breaks = seq(-6, 6, 2)) +
  geom_abline(slope = 100000, linetype = "dashed") +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  theme_bw()