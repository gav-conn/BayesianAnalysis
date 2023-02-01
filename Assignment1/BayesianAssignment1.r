setwd("C:/Users/Gavin/Documents/Bayesian Analysis")

library(tidyr) 
library(rstan) 
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
library(ggplot2)
library(gridExtra)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(shinystan)
source('stan_utility.R')
SEED <- 151515 # set random seed for reproducability

# Create data list
d_bin <- list(N1 = 300, y1 = 145, N2 = 350, y2 = 144)

writeLines(readLines("binomHypertension1.stan"))
fit_bin <- stan(file = 'binomHypertension1.stan', data = d_bin, iter = 10000, seed = SEED)

monitor(fit_bin, probs = c(0.1, 0.5, 0.9))

draws <- as.data.frame(fit_bin)

mean(draws$theta1)
var(draws$theta1)

# Calculate 95% Credible Interval for Posterior Estimate of theta of New Data
CredInts = HDInterval::hdi(draws, 0.95)
CredIntNew = CredInts[,1]
CredIntNew

mcmc_hist(draws, pars = 'theta1') +
  geom_vline(xintercept = CredIntNew[1], linetype='dotted') + # plot the 95% credible interval
  geom_vline(xintercept = CredIntNew[2], linetype='dotted')

mean(draws$theta2)
var(draws$theta2)

# Calculate 95% Credible Interval for Posterior Estimate of theta of Old Data
CredIntOld = CredInts[,2]
CredIntOld

mcmc_hist(draws, pars = 'theta2') +
  geom_vline(xintercept = CredIntOld[1], linetype='dotted') + # plot the 95% credible interval
  geom_vline(xintercept = CredIntOld[2], linetype='dotted')

mean(draws$oddsratio)
var(draws$oddsratio)

# Calculate 95% Credible Interval for Posterior Estimate of Odds Ratio
CredIntOdds = CredInts[,3]
CredIntOdds

mcmc_hist(draws, pars = 'oddsratio') +
  geom_vline(xintercept = CredIntOdds[1], linetype='dotted') + # plot the 95% credible interval
  geom_vline(xintercept = CredIntOdds[2], linetype='dotted') + 
  geom_vline(xintercept = 1, linetype='longdash', col = '2') # plot x = 1 line (odds ratio of thetas equal)

mean(draws$difference)
var(draws$difference)

# Calculate 95% Credible Interval for Posterior Estimate of Difference between thetas
CredIntDiff = CredInts[,4]
CredIntDiff

mcmc_hist(draws, pars = 'difference') +
  geom_vline(xintercept = CredIntDiff[1], linetype='dotted') + # plot the 95% credible interval
  geom_vline(xintercept = CredIntDiff[2], linetype='dotted') + 
  geom_vline(xintercept = 0, linetype='longdash', col = '2') # plot x = 0 line (difference between thetas equals 0)

writeLines(readLines("binomHypertension2.stan"))
fit_bin2 <- stan(file = 'binomHypertension2.stan', data = d_bin, iter = 10000, seed = SEED)

monitor(fit_bin, probs = c(0.1, 0.5, 0.9))

draws2 <- as.data.frame(fit_bin2)

mean(draws2$oddsratio)
var(draws2$oddsratio)

# Calculate 95% Credible Interval for Posterior Estimate of Odds Ratio
CredInts2 = HDInterval::hdi(draws2, 0.95)
CredIntOdds2 = CredInts2[,3]
CredIntOdds2

mcmc_hist(draws2, pars = 'oddsratio') +
  geom_vline(xintercept = CredIntOdds2[1], linetype='dotted') + # plot 95% credible interval
  geom_vline(xintercept = CredIntOdds2[2], linetype='dotted') + 
  geom_vline(xintercept = 1, linetype='longdash', col = '2') # plot x = 1 line (odds ratio of thetas equal)

mean(draws2$difference)
var(draws2$difference)

# Calculate 95% Credible Interval for Posterior Estimate of Difference
CredIntDiff2 = CredInts2[,4]
CredIntDiff2

mcmc_hist(draws2, pars = 'difference') +
  geom_vline(xintercept = CredIntDiff2[1], linetype='dotted') + # plot the 95% credible interval
  geom_vline(xintercept = CredIntDiff2[2], linetype='dotted') + 
  geom_vline(xintercept = 0, linetype='longdash', col = '2') # plot x = 0 line (difference between thetas is 0)
