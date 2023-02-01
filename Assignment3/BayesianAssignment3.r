setwd("C:/Users/Gavin/Documents/Bayesian Analysis")

library(tidyr) 
library(rstan) 
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(shinystan)
source('stan_utility.R')
library(GGally)
library(gridExtra)
library(loo)

SEED <- 151515 # set random seed for reproducability

# Q1
data <- read.csv("ozone.csv")
str(data)
data <- data[,-1] # Remove X variable which is just an index of observation numbers

ggpairs(data) # look at plot of the data to check for patterns

N = nrow(data)
dat1 <- list(N = NROW(data),  K =2, 
             X = cbind(rep(1,N), data$solar), ozone = data$Y)
writeLines(readLines("Ozone.stan"))
fit_1 <- stan(file="Ozone.stan", data=dat1, iter=10000)
print(fit_1, pars=c("beta", "sigma"), probs = c(0.05, 0.5, 0.95)) # rhat = 1 for each variable, chain converged

post1 <- as.data.frame(fit_1)
cor(cbind(post1$`beta[1]`,post1$`beta[2]`))
stan_plot(fit_1, pars = c("beta"))
f_mu2 <- function(x) post1$`beta[1]` + post1$`beta[2]` * x
solar_new <- seq(-2, 2, 0.01) # plot showed data lying on interval (0,300)
mu2 <- sapply(solar_new, f_mu2)

# calculate 90% credible interval for regression mean
y_hdi = HDInterval::hdi(mu2, credMass=0.9)
hpdi_l = y_hdi[1,]
hpdi_u = y_hdi[2,]

p <- ggplot() 
# store plot of regression line as well as 90% credible interval for regression mean
p2 <- p +
  geom_point(data = as.data.frame(dat1),
             aes(dat1$X[,2], dat1$ozone), shape = 1, color = 'dodgerblue') +
  geom_ribbon(aes(solar_new, ymin = hpdi_l, ymax = hpdi_u),
              alpha = .1) +
  geom_abline(data = post1,
              aes(intercept = mean(`beta[1]`), slope = mean(`beta[2]`))) +
  labs(subtitle="HPDI Interval = 0.9") + xlab("Solar")+ylab("Ozone")

y_pi <- sapply(solar_new,
               function(x) rnorm(NROW(post1), post1$`beta[1]` + post1$`beta[2]` * x, post1$sigma)
)

# calculate 90% credible intervals for prediction intervals
y_phdi = HDInterval::hdi(y_pi, credMass=0.9)
pi_l = y_phdi[1,]
pi_u = y_phdi[2,]

# plot prediction intervals along with regression line & credible interval for linear mean
p2 + geom_ribbon(mapping = aes(solar_new, ymin=pi_l, ymax=pi_u), alpha = 0.05) +
  labs(subtitle = 'Prediction Intervals = 0.9')
# We see from the graph that the model seems to do a decent job at predicting low solar value observations
# But seems to perform poorly in predicting the higher valued solar observations, perhaps indicating a non-linear relationship between the 2 variables

dat2 <- list(N = NROW(data), K =3, X = cbind(rep(1,N), data$solar, data$wind), ozone = data$Y)
fit_2 <- stan(file="Ozone.stan", data=dat2, iter=10000)
print(fit_2, pars=c("beta", "sigma"), probs = c(0.05, 0.5, 0.95)) # rhat = 1 for each variable, chain converged

post2 <- as.data.frame(fit_2)
cor(cbind(post2$`beta[1]`,post2$`beta[2]`, post2$`beta[3]`))
stan_plot(fit_2, pars = c("beta"))

f_mu2_solar <- function(x) post2$`beta[1]` + post2$`beta[2]` * x + post2$`beta[3]` * mean(data$wind)
mu2_solar <- sapply(solar_new,f_mu2_solar)
# calculate 90% credible interval for regression mean
y_hdi_solar = HDInterval::hdi(mu2_solar, credMass=0.9)
hpdi_l_solar = y_hdi_solar[1,]
hpdi_u_solar = y_hdi_solar[2,]

# store plot of regression line as well as 90% credible interval for regression mean
p2 <- p +
  geom_point(data = as.data.frame(dat2),
             aes(dat2$X[,2], dat2$ozone), shape = 1, color = 'dodgerblue') +
  geom_ribbon(aes(solar_new, ymin = hpdi_l_solar, ymax = hpdi_u_solar),
              alpha = .1) +
  geom_abline(data = post2,
              aes(intercept = mean(`beta[1]`), slope = mean(`beta[2]`))) +
  labs(subtitle="HPDI Interval = 0.9") +  xlab("solar") + ylab("Ozone")

y_pi_solar <- sapply( solar_new,
                      function(x) rnorm(NROW(post2), post2$`beta[1]` + post2$`beta[2]` * x + post2$`beta[3]`*mean(data$wind), post2$sigma)
)

# calculate 90% credible intervals for prediction intervals
y_phdi_solar = HDInterval::hdi(y_pi_solar, credMass=0.9)
pi_l_solar = y_phdi_solar[1,]
pi_u_solar = y_phdi_solar[2,]

# plot prediction intervals along with regression line & credible interval for linear mean
p2_solar <- p2 + geom_ribbon(mapping = aes(solar_new, ymin=pi_l_solar, ymax=pi_u_solar), alpha = 0.05) +
  labs(subtitle = 'Prediction Intervals = 0.9')

f_mu2_wind <- function(x) post2$`beta[1]` + post2$`beta[2]`*mean(data$solar) + post2$`beta[3]` * x
wind_new <- seq(-2.2, 3.1, 0.01)
mu2_wind <- sapply(wind_new,f_mu2_wind)
# calculate 90% credible interval for regression mean
y_hdi_wind = HDInterval::hdi(mu2_wind, credMass=0.9)
hpdi_l_wind = y_hdi_wind[1,]
hpdi_u_wind = y_hdi_wind[2,]

y_pi_wind <- sapply( wind_new,
                     function(x) rnorm(NROW(post2), post2$`beta[1]` + post2$`beta[2]`*mean(data$wind) + post2$`beta[3]` * x, post2$sigma)
)

# calculate 90% credible intervals for prediction intervals
y_phdi_wind = HDInterval::hdi(y_pi_wind, credMass=0.9)
pi_l_wind = y_phdi_wind[1,]
pi_u_wind = y_phdi_wind[2,]

# store plot of regression line as well as 90% credible interval for regression mean
p2 <- p +
  geom_point(data = as.data.frame(dat2),
             aes(dat2$X[,3], dat2$ozone), shape = 1, color = 'dodgerblue') +
  geom_ribbon(aes(wind_new, ymin = hpdi_l_wind, ymax = hpdi_u_wind),
              alpha = .1) +
  geom_abline(data = post2,
              aes(intercept = mean(`beta[1]`), slope = mean(`beta[3]`))) +
  labs(subtitle="HPDI Interval = 0.9") +  xlab("wind") + ylab("Ozone")

# plot prediction intervals along with regression line & credible interval for linear mean
p2_wind <- p2 + geom_ribbon(mapping = aes(wind_new, ymin=pi_l_wind, ymax=pi_u_wind), alpha = 0.05) +
  labs(subtitle = 'Prediction Intervals = 0.9')

grid.arrange(p2_solar,p2_wind, nrow = 2)
p2_wind
dat3 <- list(N = NROW(data), K=4, X = cbind(rep(1,N), data$solar, data$wind, data$temp), ozone = data$Y)
fit_3 <- stan(file="Ozone.stan", data=dat3, iter=10000)
print(fit_3, pars=c("beta", "sigma"), probs =c(0.05, 0.5, 0.95)) # rhat = 1 for each variable, chain converged

post3 <- as.data.frame(fit_3)
cor(cbind(post3$`beta[1]`,post3$`beta[2]`, post3$`beta[3]`, post3$`beta[4]`))
stan_plot(fit_3, pars = c("beta"))

f_mu2_solar <- function(x) post3$`beta[1]` + post3$`beta[2]` * x +
  post3$`beta[3]` * mean(data$wind) + post3$`beta[4]` * mean(data$temp)
mu2_solar <- sapply(solar_new,f_mu2_solar)
# calculate 90% credible interval for regression mean
y_hdi_solar = HDInterval::hdi(mu2_solar, credMass=0.9)
hpdi_l_solar = y_hdi_solar[1,]
hpdi_u_solar = y_hdi_solar[2,]

# store plot of regression line as well as 90% credible interval for regression mean
p2 <- p +
  geom_point(data = as.data.frame(dat3),
             aes(dat3$X[,2], dat3$ozone), shape = 1, color = 'dodgerblue') +
  geom_ribbon(aes(solar_new, ymin = hpdi_l_solar, ymax = hpdi_u_solar),
              alpha = .1) +
  geom_abline(data = post3,
              aes(intercept = mean(`beta[1]`), slope = mean(`beta[2]`))) +
  labs(subtitle="HPDI Interval = 0.9")+  xlab("solar") + ylab("Ozone")

y_pi_solar <- sapply( solar_new,
                      function(x) rnorm(NROW(post3), post3$`beta[1]` + post3$`beta[2]` * x + post3$`beta[3]`*mean(data$wind) + post3$`beta[4]`*mean(data$temp), post3$sigma)
)

# calculate 90% credible intervals for prediction intervals
y_phdi_solar = HDInterval::hdi(y_pi_solar, credMass=0.9)
pi_l_solar = y_phdi_solar[1,]
pi_u_solar = y_phdi_solar[2,]

# plot prediction intervals along with regression line & credible interval for linear mean
p2_solar <- p2 + geom_ribbon(mapping = aes(solar_new, ymin=pi_l_solar, ymax=pi_u_solar), alpha = 0.05) +
  labs(subtitle = 'Prediction Intervals = 0.9')

f_mu2_wind <- function(x) post3$`beta[1]` + post3$`beta[2]`*mean(data$solar) + post3$`beta[3]` * x + post3$`beta[4]` * mean(data$temp)
mu2_wind <- sapply(wind_new,f_mu2_wind)
# calculate 90% credible interval for regression mean
y_hdi_wind = HDInterval::hdi(mu2_wind, credMass=0.9)
hpdi_l_wind = y_hdi_wind[1,]
hpdi_u_wind = y_hdi_wind[2,]

y_pi_wind <- sapply( wind_new,
                     function(x) rnorm(NROW(post3), post3$`beta[1]` + post3$`beta[2]`*mean(data$wind) +
                                         post3$`beta[3]` * x + post3$`beta[4]`*mean(data$temp), post3$sigma)
)

# calculate 90% credible intervals for prediction intervals
y_phdi_wind = HDInterval::hdi(y_pi_wind, credMass=0.9)
pi_l_wind = y_phdi_wind[1,]
pi_u_wind = y_phdi_wind[2,]

# store plot of regression line as well as 90% credible interval for regression mean
p2 <- p +
  geom_point(data = as.data.frame(dat3),
             aes(dat3$X[,3], dat3$ozone), shape = 1, color = 'dodgerblue') +
  geom_ribbon(aes(wind_new, ymin = hpdi_l_wind, ymax = hpdi_u_wind),
              alpha = .1) +
  geom_abline(data = post3,
              aes(intercept = mean(`beta[1]`), slope = mean(`beta[3]`))) +
  labs(subtitle="HPDI Interval = 0.9") + xlab("wind") + ylab("Ozone")

# plot prediction intervals along with regression line & credible interval for linear mean
p2_wind <- p2 + geom_ribbon(mapping = aes(wind_new, ymin=pi_l_wind, ymax=pi_u_wind), alpha = 0.05) +
  labs(subtitle = 'Prediction Intervals = 0.9')

f_mu2_temp <- function(x) post3$`beta[1]` + post3$`beta[2]`*mean(data$solar) + post3$`beta[3]` * mean(data$wind)+ post3$`beta[4]`*x 
temp_new <- seq(-2.4, 2.2, 0.01)
mu2_temp <- sapply(temp_new,f_mu2_temp)
# calculate 90% credible interval for regression mean
y_hdi_temp = HDInterval::hdi(mu2_temp, credMass=0.9)
hpdi_l_temp = y_hdi_temp[1,]
hpdi_u_temp = y_hdi_temp[2,]

y_pi_temp <- sapply( temp_new,
                     function(x) rnorm(NROW(post3), post3$`beta[1]` + post3$`beta[2]`*mean(data$wind) + post3$`beta[3]` *mean(data$wind)+post3$`beta[4]`*x, post3$sigma)
)

# calculate 90% credible intervals for prediction intervals
y_phdi_temp = HDInterval::hdi(y_pi_temp, credMass=0.9)
pi_l_temp = y_phdi_temp[1,]
pi_u_temp = y_phdi_temp[2,]

# store plot of regression line as well as 90% credible interval for regression mean
p2 <- p +
  geom_point(data = as.data.frame(dat3),
             aes(dat3$X[,4], dat3$ozone), shape = 1, color = 'dodgerblue') +
  geom_ribbon(aes(temp_new, ymin = hpdi_l_temp, ymax = hpdi_u_temp),
              alpha = .1) +
  geom_abline(data = post3,
              aes(intercept = mean(`beta[1]`), slope = mean(`beta[4]`))) +
  labs(subtitle="HPDI Interval = 0.9") + xlab("temperature") + ylab("Ozone")

# plot prediction intervals along with regression line & credible interval for linear mean
p2_temp <- p2 + geom_ribbon(mapping = aes(temp_new, ymin=pi_l_temp, ymax=pi_u_temp), alpha = 0.05) +
  labs(subtitle = 'Prediction Intervals = 0.9')

p2_solar
p2_wind
p2_temp

# Q2
log_lik1 <- extract_log_lik(fit_1)
waic1 <- waic(log_lik1)
waic1
loo1 <-loo(fit_1, pars="log_lik", cores = 2) # provides same output as above
loo1 # pareto k estimates are good meaning indicating estimates are reliable

log_lik2 <- extract_log_lik(fit_2)
waic2 <- waic(log_lik2)
waic2
loo2 <- loo(fit_2, pars="log_lik", cores = 2)
loo2 # pareto k estimates are good indicating estimates are reliable

log_lik3 <- extract_log_lik(fit_3)
waic3 <- waic(log_lik3)
waic3
loo3 <- loo(fit_3, pars="log_lik", cores = 2)
loo3 # pareto k estimates are ok indicating estimates are reliable

loo_compare(waic1, waic2, waic3)
loo_compare(loo1, loo2, loo3)
        
# Q3 
y_rep1 <- as.matrix(fit_1, pars=c("y_rep"))
y_rep2 <- as.matrix(fit_2, pars=c("y_rep"))
y_rep3 <- as.matrix(fit_3, pars=c("y_rep"))

#plot 5 random samples of the posterior predictive distribution vs observed distribution
i = sample(1:nrow(data), size = 5)
q1 <- ggplot() + geom_density(aes(x=y_rep1[i[1],])) + theme_classic()
q2 <- ggplot() + geom_density(aes(x=y_rep1[i[2],])) + theme_classic()
q3 <- ggplot() + geom_density(aes(x=y_rep1[i[3],])) + theme_classic()
q4 <- ggplot() + geom_density(aes(x=y_rep1[i[4],])) + theme_classic()
q5 <- ggplot() + geom_density(aes(x=y_rep1[i[5],])) + theme_classic()
q6 <- ggplot() + geom_density(aes(x=dat1$ozone), color="blue") + theme_classic()

grid.arrange(q1,q2,q3,q4,q5,q6, nrow = 2)

q1 <- ggplot() + geom_density(aes(x=y_rep2[i[1],])) + theme_classic()
q2 <- ggplot() + geom_density(aes(x=y_rep2[i[2],])) + theme_classic()
q3 <- ggplot() + geom_density(aes(x=y_rep2[i[3],])) + theme_classic()
q4 <- ggplot() + geom_density(aes(x=y_rep2[i[4],])) + theme_classic()
q5 <- ggplot() + geom_density(aes(x=y_rep2[i[5],])) + theme_classic()
q6 <- ggplot() + geom_density(aes(x=dat2$ozone), color="blue") + theme_classic()

grid.arrange(q1,q2,q3,q4,q5,q6, nrow = 2)

q1 <- ggplot() + geom_density(aes(x=y_rep3[i[1],])) + theme_classic()
q2 <- ggplot() + geom_density(aes(x=y_rep3[i[2],])) + theme_classic()
q3 <- ggplot() + geom_density(aes(x=y_rep3[i[3],])) + theme_classic()
q4 <- ggplot() + geom_density(aes(x=y_rep3[i[4],])) + theme_classic()
q5 <- ggplot() + geom_density(aes(x=y_rep3[i[5],])) + theme_classic()
q6 <- ggplot() + geom_density(aes(x=dat3$ozone), color="blue") + theme_classic()

grid.arrange(q1,q2,q3,q4,q5,q6, nrow = 2)

# Take random sample of 50 posterior distributions & compare to observed data for each model
ppc_dens_overlay(data$Y, y_rep1[sample(1:length(y_rep1[1,]), size = 50),])+theme_classic()
ppc_dens_overlay(data$Y, y_rep2[sample(1:length(y_rep2[1,]), size = 50),])+theme_classic()
ppc_dens_overlay(data$Y, y_rep3[sample(1:length(y_rep3[1,]), size = 50),])+theme_classic()

# Q4
dat4 <- list(N = NROW(data), K = 4, X = cbind(rep(1,N), data$solar, data$wind, data$temp), ozone = data$Y)
# Use lognormal model for ozone rather than normal
writeLines(readLines("Ozone2.stan"))
fit_4 <- stan(file="Ozone2.stan", data=dat4, iter=10000)
print(fit_4, pars=c("beta", "sigma"), probs = c(0.05, 0.5, 0.95)) # rhat = 1 for each variable, chain converged

post4 <- as.data.frame(fit_4)
cor(cbind(post4$`beta[1]`,post4$`beta[2]`, post4$`beta[3]`, post4$`beta[4]`, post4$`beta[5]`))
stan_plot(fit_4, pars = c("beta"))

y_rep4 <- as.matrix(fit_4, pars=c("y_rep"))
ppc_dens_overlay(data$Y, y_rep4[sample(1:length(y_rep4[1,]), size = 50),])+theme_classic()
# Fit seems much better than previous models

# include squared solar parameter within model, as linear fit does not fit data well
dat5 <- list(N = NROW(data), K = 5, X = cbind(rep(1,N), data$solar, data$wind, data$temp, data$solar^2), ozone = data$Y)
fit_5 <- stan(file="Ozone2.stan", data=dat5, iter=10000)
print(fit_5, pars=c("beta", "sigma"), probs = c(0.05, 0.5, 0.95)) # rhat = 1 for each variable, chain converged

post5 <- as.data.frame(fit_5)
cor(cbind(post5$`beta[1]`,post5$`beta[2]`, post5$`beta[3]`, post5$`beta[4]`, post5$`beta[5]`))
stan_plot(fit_5, pars = c("beta"))

y_rep5 <- as.matrix(fit_5, pars=c("y_rep"))
ppc_dens_overlay(data$Y, y_rep5[sample(1:length(y_rep5[1,]), size = 50),])+theme_classic()

log_lik4 <- extract_log_lik(fit_4)
waic4 <- waic(log_lik4)
loo4 <- loo(fit_4, pars="log_lik", cores = 2)
loo4

log_lik5 <- extract_log_lik(fit_5)
waic5 <- waic(log_lik5)
loo5 <- loo(fit_5, pars="log_lik", cores = 2)
loo5

loo_compare(waic1, waic2, waic3, waic4, waic5)
loo_compare(loo1, loo2, loo3, loo4, loo5)