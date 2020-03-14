# Required R-packages
# install.packages("devtools")
# devtools::install_github("hesim-dev/hesim")
# install.packages("data.table")
# install.packages("expm")
library(expm)
library(hesim)
library(data.table)

rm(list=ls())



# Model settings
initial<-matrix(c(1,0,0,0),nrow=1)

follow_up <- 20

discount_LY <- 0.00
discount_cost <- 0.06

PSA_samples  <- 1000


# Model parameters
## Transition matrix monotherapy
transitions_mono <- matrix(c(1251, 350, 116, 17,
                             0, 731, 512, 15,
                             0, 0, 1312, 437,
                             0, 0, 0, 469),
                             nrow = 4, byrow = TRUE)
p_mono = rdirichlet_mat(n=PSA_samples,transitions_mono)

print(p_mono[, , 1:10])

## Relative treatment effect with combination therapy
RR<-0.509
RR_low<-0.365
RR_high<-0.710

RR_combi <- exp( rnorm(PSA_samples, log(RR), (log(RR_high)-log(RR_low))/3.92 )   )

RR_combi[1:10]

## Transition matrix combination therapy
p_combi <- array(NA, dim = c(4, 4, PSA_samples))

for (i in 1:PSA_samples){ 
p_combi[,,i] <- matrix(c((1-p_mono[1,2,i]*RR_combi[i]-p_mono[1,3,i]*RR_combi[i]-p_mono[1,4,i]*RR_combi[i]), p_mono[1,2,i]*RR_combi[i], p_mono[1,3,i]*RR_combi[i], p_mono[1,4,i]*RR_combi[i],
                                      0, (1-p_mono[2,3,i]*RR_combi[i]-p_mono[2,4,i]*RR_combi[i]), p_mono[2,3,i]*RR_combi[i], p_mono[2,4,i]*RR_combi[i],
                                      0, 0, (1-p_mono[3,4,i]*RR_combi[i]), p_mono[3,4,i]*RR_combi[i],
                                      0, 0, 0, 1),
                                      ncol = 4, nrow = 4, byrow = TRUE)
}

print(p_combi[, , 1:10])

## Other parameters
LY <- matrix(c(1, 1, 1, 0),nrow=1)

drug_cost_mono <- matrix(c(2278, 2278, 2278, 0),nrow=1)
drug_cost_combi <- matrix(c(2278+2087, 2278+2087, 2278+2087, 0),nrow=1)
other_cost <- matrix(c(2756, 3052, 9007, 0),nrow=1)


# Model output
## Monotherapy
trace_mono <- array(NA, dim = c(follow_up, 6, PSA_samples))
colnames(trace_mono) <- c("State A", "State B", "State C", "Death", "Disc LY", "Disc cost")

Expected_LY_mono <- matrix(rep(NA, PSA_samples), nrow=PSA_samples)
Expected_total_cost_mono <- matrix(rep(NA, PSA_samples), nrow=PSA_samples)

for(i in 1:PSA_samples){
for(t in 1:follow_up){
trace_mono[t,1:4,i] <- initial %*% (p_mono[,,i] %^%t)
trace_mono[t,5,i] <- sum(  trace_mono[t,1:4,i] * LY * 1/(1+discount_LY)^t  )
trace_mono[t,6,i] <- sum(  trace_mono[t,1:4,i] * (drug_cost_mono + other_cost) * 1/(1+discount_cost)^t  )
}
Expected_LY_mono[i] <- sum(trace_mono[,5,i])
Expected_total_cost_mono[i] <- sum(trace_mono[,6,i])
}

trace_mono[,,1:10]
Expected_LY_mono[1:10]
Expected_total_cost_mono[1:10]


## Combination therapy
trace_combi <- array(NA, dim = c(follow_up, 6, PSA_samples))
colnames(trace_combi) <- c("State A", "State B", "State C", "Death", "Disc LY", "Disc cost")

Expected_LY_combi <- matrix(rep(NA, PSA_samples), nrow=PSA_samples)
Expected_total_cost_combi <- matrix(rep(NA, PSA_samples), nrow=PSA_samples)

for(i in 1:PSA_samples){
for(t in 1:2){
trace_combi[t,1:4,i] <- initial %*% (p_combi[,,i] %^%t)
trace_combi[t,5,i] <- sum(  trace_combi[t,1:4,i] * LY * 1/(1+discount_LY)^t  )
trace_combi[t,6,i] <- sum(  trace_combi[t,1:4,i] * (drug_cost_combi + other_cost) * 1/(1+discount_cost)^t  )
}
for(t in 3:follow_up){
trace_combi[t,1:4,i] <- trace_combi[2,1:4,i] %*% (p_mono[,,i] %^% (t-2))
trace_combi[t,5,i] <- sum(  trace_combi[t,1:4,i] * LY * 1/(1+discount_LY)^t  )
trace_combi[t,6,i] <- sum(  trace_combi[t,1:4,i] * (drug_cost_mono + other_cost) * 1/(1+discount_cost)^t  )
}
Expected_LY_combi[i] <- sum(trace_combi[,5,i])
Expected_total_cost_combi[i] <- sum(trace_combi[,6,i])
}

trace_combi[,,1:10]
Expected_LY_combi[1:10]
Expected_total_cost_combi[1:10]


## Incremental Cost-Effectiveness Ratio
ICER <- matrix(rep(NA, PSA_samples), nrow=PSA_samples)
for(i in 1:PSA_samples){
ICER[i]  <- (Expected_total_cost_combi[i] - Expected_total_cost_mono[i]) / (Expected_LY_combi[i] - Expected_LY_mono[i])
}

ICER[1:10]

## Overview PSA results
Results_by_PSA_sample <- cbind(Expected_LY_mono, Expected_LY_combi, Expected_total_cost_mono, Expected_total_cost_combi, ICER)
colnames(Results_by_PSA_sample) <- c("Disc LY mono", "Disc LY combi", "Disc cost mono", "Disc cost combi", "ICER")

Results_by_PSA_sample[1:10,]

