## ---- Overview ---------------------------------------------------------------
## @knitr R-packages
library("rcea")
library("knitr")
library("kableExtra")
library("magrittr")
library("tibble")

## ---- Model parameters -------------------------------------------------------
## @knitr transition-probabilities
p_hd <- 0.002 # constant probability of dying when Healthy (all-cause mortality)
p_hs1 <- 0.15 # probability of becoming Sick when Healthy
p_s1h <- 0.5 # probability of becoming Healthy when Sick
p_s1s2 <- 0.105 # probability of becoming Sicker when Sick
p_s1d <- .006 # constant probability of dying when Sick
p_s2d <- .02 # constant probability of dying when Sicker

## @knitr transition-probability-complements
p_hh <- 1  - p_hs1 - p_hd
p_s1s1 <- 1 - p_s1h - p_s1s2 
p_s2s2 <- 1 - p_s2d

## @knitr tpmatrix
p_soc <- matrix(
  c(p_hh,  p_hs1,  0,      p_hd,
    p_s1h, p_s1s1, p_s1s2, p_s1d,
    0,     0,      p_s2s2, p_s2d,
    0,     0,      0,      1),
  byrow = TRUE,
  nrow = 4, ncol = 4
)
state_names <- c("H", "S1", "S2", "D")
colnames(p_soc) <- rownames(p_soc) <- state_names
print(p_soc)

## @knitr apply_rr
apply_rr <- function(p, rr = .8){
  p["H", "S1"] <- p["H", "S1"] * rr
  p["H", "S2"] <- p["H", "S2"] * rr
  p["H", "D"] <- p["H", "S2"] * rr
  p["H", "H"] <- 1 - sum(p["H", -1])
  
  p["S1", "S2"] <- p["S1", "S2"] * rr
  p["S1", "D"] <- p["S1", "D"] * rr
  p["S1", "S1"] <- 1 - sum(p["S1", -2])
  
  p["S2", "D"] <- p["S2", "D"] * rr
  p["S2", "S2"] <- 1 - sum(p["S2", -3])
  
  return(p)
}
p_new <- apply_rr(p_soc, rr = .8)

## @knitr utility-costs
utility <- c(1, .75, 0.5, 0)
costs_medical <- c(2000, 4000, 15000, 0)
costs_treat_soc <- c(rep(2000, 3), 0)
costs_treat_new <- c(rep(12000, 3), 0)

## ---- Simulation -------------------------------------------------------------
## @knitr x_init
x_init <- c(1, 0, 0, 0)

## @knitr sim-markov-chain
x_soc <- sim_markov_chain(x_init, p_soc)
x_new <- sim_markov_chain(x_init, p_new)

head(x_soc)
head(x_new)

## @knitr compute_qalys
compute_qalys <- function(x, utility, dr = .03){
  n_cycles <- nrow(x) - 1
  pv(x %*% utility, dr, 0:n_cycles)
}

qalys_soc <- x_soc %*% utility 
dqalys_soc <- compute_qalys(x_soc, utility = utility)
dqalys_new <- compute_qalys(x_new, utility = utility)

head(qalys_soc)
head(dqalys_soc)
head(dqalys_new)

## @knitr compute_costs
compute_costs <- function(x, costs_medical, costs_treat, dr = .03){
  n_cycles <- nrow(x) - 1
  costs <- cbind(
    pv(x %*% costs_medical, dr, 0:n_cycles),
    pv(x %*% costs_treat, dr, 0:n_cycles)
  )
  colnames(costs) <- c("medical", "treatment")
  return(costs)
}

dcosts_soc <- compute_costs(x_soc, costs_medical, costs_treat_soc)
dcosts_new <- compute_costs(x_new, costs_medical, costs_treat_new)

head(dcosts_soc)
head(dcosts_new)

## ---- Cost-effectiveness analysis --------------------------------------------
## @knitr icer
(sum(dcosts_new) - sum(dcosts_soc))/(sum(dqalys_new) - sum(dqalys_soc))

## @knitr icer_tbl
format_costs <- function(x) formatC(x, format = "d", big.mark = ",")
format_qalys <- function(x) formatC(x, format = "f", digits = 2)
make_icer_tbl <- function(costs0, costs1, qalys0, qalys1){
  # Computations
  total_costs0 <- sum(costs0)
  total_costs1 <- sum(costs1)
  total_qalys0 <- sum(qalys0)
  total_qalys1 <- sum(qalys1)
  incr_total_costs <- total_costs1 - total_costs0
  inc_total_qalys <- total_qalys1 - total_qalys0
  icer <- incr_total_costs/inc_total_qalys
  
  # Make table
  tibble(
    `Costs` = c(total_costs0, total_costs1) %>%
      format_costs(), 
    `QALYs` = c(total_qalys0, total_qalys1) %>%
      format_qalys(),
    `Incremental costs` = c("--", incr_total_costs %>% 
                             format_costs()),
    `Incremental QALYs` = c("--", inc_total_qalys %>% 
                             format_qalys()),
    `ICER` = c("--", icer %>% format_costs())
  ) %>%
    kable() %>%
    kable_styling() %>%
    footnote(general = "Costs and QALYs are discounted at 3% per annum.",
             footnote_as_chunk = TRUE)
}
make_icer_tbl(costs0 = dcosts_soc, costs1 = dcosts_new,
              qalys0 = dqalys_soc, qalys1 = dqalys_new)

## ---- Exercises --------------------------------------------------------------