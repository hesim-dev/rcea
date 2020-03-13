# Required R-package 
#install.packages("expm",repos = "http://cran.us.r-project.org")
library(expm)



# Model settings
initial<-matrix(c(1,0,0,0),nrow=1)
initial

follow_up <- 20

discount_LY <- 0.00
discount_cost <- 0.06


# Model parameters
## Transition matrix monotherapy
transitions_mono <- matrix(c(0.721, 0.202, 0.067, 0.010,
                        0, 0.581, 0.407, 0.012,
                        0, 0, 0.750, 0.250,
                        0, 0, 0, 1), 
                        ncol = 4, nrow = 4, byrow = TRUE)
transitions_mono

## Relative treatment effect with combination therapy
RR<-0.509

## Transition matrix combination therapy
transitions_combi <- matrix(c((1-0.202*RR-0.067*RR-0.010*RR), 0.202*RR, 0.067*RR, 0.010*RR,
                                      0, (1-0.407*RR-0.012*RR), 0.407*RR, 0.012*RR,
                                      0, 0, (1-0.250*RR), 0.250*RR,
                                      0, 0, 0, 1),
                                      ncol = 4, nrow = 4, byrow = TRUE)
transitions_combi

## Lifeyears and costs by health state
LY <- matrix(c(1, 1, 1, 0),nrow=1)

drug_cost_mono <- matrix(c(2278, 2278, 2278, 0),nrow=1)
drug_cost_combi <- matrix(c(2278+2087, 2278+2087, 2278+2087, 0),nrow=1)
other_cost <- matrix(c(2756, 3052, 9007, 0),nrow=1)

LY
drug_cost_mono
drug_cost_combi
other_cost



# Basic example calculations
## Health state occupancy with monotherapy
initial %*% transitions_mono
initial %*% transitions_mono %*% transitions_mono %*% transitions_mono %*% transitions_mono %*% transitions_mono
initial %*% (transitions_mono %^% 5)

## Expected life years and costs with monotherapy
expected_LY_mono <- sum(
            initial %*% (transitions_mono %^% 1) * LY * 1/(1+discount_LY)^1,
            initial %*% (transitions_mono %^% 2) * LY * 1/(1+discount_LY)^2,
            initial %*% (transitions_mono %^% 3) * LY * 1/(1+discount_LY)^3,
            initial %*% (transitions_mono %^% 4) * LY * 1/(1+discount_LY)^4,
            initial %*% (transitions_mono %^% 5) * LY * 1/(1+discount_LY)^5 
            )
expected_LY_mono

expected_total_cost_mono <- sum(
            initial %*% (transitions_mono %^% 1) * (drug_cost_mono + other_cost) * 1/(1+discount_cost)^1,
            initial %*% (transitions_mono %^% 2) * (drug_cost_mono + other_cost) * 1/(1+discount_cost)^2,
            initial %*% (transitions_mono %^% 3) * (drug_cost_mono + other_cost) * 1/(1+discount_cost)^3,
            initial %*% (transitions_mono %^% 4) * (drug_cost_mono + other_cost) * 1/(1+discount_cost)^4,
            initial %*% (transitions_mono %^% 5) * (drug_cost_mono + other_cost) * 1/(1+discount_cost)^5  
            ) 
expected_total_cost_mono



# Model output
## Monotherapy
trace_mono<-matrix(rep(NA, 6 * follow_up), nrow=follow_up)  
colnames(trace_mono) <- c("State A", "State B", "State C", "Death", "Disc LY", "Disc cost") 

for(t in 1:follow_up){
trace_mono[t,1:4] <- initial %*% (transitions_mono %^%t)
trace_mono[t,5] <- sum(  trace_mono[t,1:4] * LY * 1/(1+discount_LY)^t  )
trace_mono[t,6] <- sum(  trace_mono[t,1:4] * (drug_cost_mono + other_cost) * 1/(1+discount_cost)^t  )
}
trace_mono

Expected_LY_mono <- sum(trace_mono[,5])
Expected_LY_mono
Expected_total_cost_mono <- sum(trace_mono[,6])
Expected_total_cost_mono


## Combination therapy
trace_combi<-matrix(rep(NA, 6 * follow_up), nrow=follow_up)
colnames(trace_combi) <- c("State A", "State B", "State C", "Death", "Disc LY", "Disc cost")

for(t in 1:2){
  trace_combi[t,1:4] <- initial %*% (transitions_combi %^%t)
  trace_combi[t,5] <- sum(  trace_combi[t,1:4] * LY * 1/(1+discount_LY)^t  )
  trace_combi[t,6] <- sum(  trace_combi[t,1:4] * (drug_cost_combi + other_cost) * 1/(1+discount_cost)^t  )
}
for(t in 3:follow_up){
  trace_combi[t,1:4] <- trace_combi[2,1:4] %*% (transitions_mono %^% (t-2))
  trace_combi[t,5] <- sum(  trace_combi[t,1:4] * LY * 1/(1+discount_LY)^t  )
  trace_combi[t,6] <- sum(  trace_combi[t,1:4] * (drug_cost_mono + other_cost) * 1/(1+discount_cost)^t  )
}
trace_combi

Expected_LY_combi <- sum(trace_combi[,5])
Expected_LY_combi
Expected_total_cost_combi <- sum(trace_combi[,6])
Expected_total_cost_combi

## Incremental Cost-Effectiveness Ratio
ICER  <- (Expected_total_cost_combi - Expected_total_cost_mono) / (Expected_LY_combi - Expected_LY_mono)
ICER

