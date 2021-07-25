## ---- Overview ---------------------------------------------------------------
## @knitr R-setup
library("hesim")
library("ggplot2")
library("magrittr")
theme_set(theme_minimal()) # Set ggplot2 theme

## ---- Application ------------------------------------------------------------
## @knitr conduct-cea
# First load the output saved during the "Semi-Markov 
# Multi-state Model" tutorial (i.e., from 04-mstate.Rmd)
load("04-mstate.Rdata")
head(ce_sim)

wtp <- seq(0, 250000, 500) # Willingness to pay per QALY
cea_pw_out <- cea_pw(ce_sim, comparator = 1, # Comparator is SOC (ID = 1)
                     dr_qalys = 0.03, dr_costs = 0.03,
                     wtp)
cea_out <- cea(ce_sim, 
               dr_qalys = 0.03, dr_costs = 0.03,
               k = wtp)

## ---- Incremental cost-effectiveness ratio -----------------------------------
## @knitr icer
labs <- get_labels(hesim_dat)
icer(cea_pw_out, wtp = 50000, labels = labs) %>%
  format()

## ---- Cost-effectiveness plane -----------------------------------------------
## @knitr ceplane-plot
plot_ceplane(cea_pw_out, k = 100000, labels = labs)

## ---- Cost-effectiveness acceptability curves --------------------------------
## @knitr ceac-simultaneous-plot
plot_ceac(cea_out, labels = labs)

## @knitr ceac-pairwise-plot
plot_ceac(cea_pw_out, labels = labs)

## ---- Cost-effectiveness acceptability frontier ------------------------------
## @knitr ceaf-plot
plot_ceaf(cea_out, labels = labs)

## ---- Value of perfect information -------------------------------------------
## @knitr evpi-plot
plot_evpi(cea_out)

## ---- Custom ggplot for cost-effectiveness plane -----------------------------
## @knitr helper-functions
strategy_factor <- function (x) { 
  factor(x, levels = 1:2, labels = c("SOC", "New"))
}

format_dollar <- function(x) {
  paste0("$", formatC(x, format = "d", big.mark = ","))
}

## @knitr custom-ceplane-plot
ylim <- max(cea_pw_out$delta[, ic]) * 1.1
xlim <- ceiling(max(cea_pw_out$delta[, ie]) * 1.1)
ggplot(cea_pw_out$delta, 
       aes(x = ie, y = ic, col = strategy_factor(strategy_id))) + 
  geom_jitter(size = .5)  + 
  xlab("Incremental QALYs") + 
  ylab("Incremental cost") +
  scale_y_continuous(limits = c(-ylim, ylim),
                     labels = format_dollar) +
  scale_x_continuous(limits = c(-xlim, xlim), breaks = seq(-6, 6, 2)) +
  theme(legend.position = "bottom") + 
  scale_colour_discrete(name = "Strategy") +
  geom_abline(slope = 100000, linetype = "dashed") +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0)

## ---- Custom ggplot for cost-effectiveness acceptability curve ---------------
## @knitr custom-ceac-simultaneous-plot
ggplot(cea_out$mce, 
       aes(x = k, y = prob, col = strategy_factor(strategy_id))) +
  geom_line() + 
  xlab("Willingness to pay") +
  ylab("Probability most cost-effective") +
  scale_x_continuous(breaks = seq(0, max(wtp), length.out = 6), 
                     label = format_dollar) +
  theme(legend.position = "bottom") + 
  scale_colour_discrete(name = "Strategy")

## @knitr custom-ceac-pairwise-plot
ggplot(cea_pw_out$ceac, 
       aes(x = k, y = prob, col = strategy_factor(strategy_id))) +
  geom_line()  + 
  xlab("Willingness to pay") +
  ylab("Probability most cost-effective") +
  scale_x_continuous(breaks =  seq(0, max(wtp), length.out = 6), 
                     label = format_dollar) +
  theme(legend.position = "bottom") +
  scale_colour_discrete(name = "Strategy")

# ---- Custom ggplot for cost-effectiveness acceptability frontier -------------
## @knitr custom-ceaf-plot
ggplot(cea_out$mce[best == 1], 
       aes(x = k, y = prob, col = strategy_factor(strategy_id))) +
  geom_line() + 
  xlab("Willingness to pay") +
  ylab("Probability most cost-effective") +
  scale_x_continuous(breaks = seq(0, max(wtp), length.out = 6), 
                     label = format_dollar) +
  theme(legend.position = "bottom") + 
  scale_colour_discrete(name = "Strategy")

# ---- Custom ggplot for value of perfect information --------------------------
## @knitr custom-evpi-plot
ggplot(cea_out$evpi, aes(x = k, y = evpi)) +
  geom_line()  + 
  xlab("Willingness to pay") +
  ylab("Expected value of perfect information") +
  scale_x_continuous(breaks = seq(0, max(wtp), length.out = 6), 
                     label = format_dollar) +
  scale_y_continuous(label = format_dollar) +
  theme(legend.position = "bottom") 