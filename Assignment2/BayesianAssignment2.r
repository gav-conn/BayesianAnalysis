setwd("C:/Users/Gavin/Documents/Bayesian Analysis")

library(tidyr) 
library(rstan) 
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(ggplot2)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(shinystan)
source('stan_utility.R')
SEED <- 151515 # set random seed for reproducability

dublinproperty <- read.csv("property.csv")

# Q1
str(dublinproperty) # Look at summary of data
area <- dublinproperty$area
price <- dublinproperty$price
plot(area, price) # look at plot of the data to check for patterns
# looks as if a regression line through data would cross y-axis at price of roughly 300
# may be a good centre for our prior on the alpha/intercept term of our model
d <- dublinproperty
  
dat <- list(N = NROW(dublinproperty), area = area, price = price)

writeLines(readLines("LinRegModel.stan"))
fit_1 <- stan(file="LinRegModel.stan", data=dat, iter=5000)

monitor(fit_1) # rhat = 1 for each variable, chain converged

post <- as.data.frame(fit_1)
cor(cbind(post$alpha,post$beta)) # Correlation between alpha and beta reveals strong correlation
f_mu <- function(x) post$alpha + post$beta * x
area_new <- seq(0, 300) # plot showed data lying on interval (0,300)
mu1 <- sapply(area_new, f_mu)

# calculate 90% credible interval for regression mean
y_hdi = HDInterval::hdi(mu1, credMass=0.9)
hpdi_l = y_hdi[1,]
hpdi_u = y_hdi[2,]

p <- ggplot() 
# store plot of regression line as well as 90% credible interval for regression mean
p2 <- p +
  geom_point(data = d,
             aes(area, price), shape = 1, color = 'dodgerblue') +
  geom_ribbon(aes(area_new, ymin = hpdi_l, ymax = hpdi_u),
              alpha = .1) +
  geom_abline(data = post,
              aes(intercept = mean(alpha), slope = mean(beta))) +
  labs(subtitle="HPDI Interval = 0.9")

y_pi <- sapply(area_new,
               function(x) rnorm(NROW(post), post$alpha + post$beta * x, post$sigma)
)

# calculate 90% credible intervals for prediction intervals
y_phdi = HDInterval::hdi(y_pi, credMass=0.9)
pi_l = y_phdi[1,]
pi_u = y_phdi[2,]

# plot prediction intervals along with regression line & credible interval for linear mean
p2 + geom_ribbon(mapping = aes(area_new, ymin=pi_l, ymax=pi_u), alpha = 0.05) +
  labs(subtitle = 'Prediction Intervals = 0.9')
# plot shows data far outside of 90% prediction interval at beginning/middle of data range 
# data does not follow regression line very closely

# Q2
logarea <- log(area)
d2 <- as.data.frame(matrix(data = c(logarea, price), nrow = 200, ncol = 2, byrow = FALSE))
dat2 <- list(N = NROW(dublinproperty), logarea = logarea, price = dublinproperty$price)

plot(logarea, price)
# harder to see where regression line would cross y-axis in this case
# prior with high variance centred around -1000 seems appropriate

writeLines(readLines("LogLinRegModel.stan"))
fit_2 <- stan(file="LogLinRegModel.stan", data=dat2, iter=5000)

monitor(fit_2) #rhat = 1 for all variables, meaning chain has converged

post <- as.data.frame(fit_2)
cor(cbind(post$alpha,post$beta)) # Correlation between alpha and beta reveals strong correlation between each. 

f_mu <- function(x) post$alpha + post$beta * x
area_new <- seq(3, 6) # plot seemed to have data lying on interval (3,6)
mu1 <- sapply(area_new, f_mu)

# calculate 90% credible interval for regression mean
y_hdi = HDInterval::hdi(mu1, credMass=0.9)
hpdi_l = y_hdi[1,]
hpdi_u = y_hdi[2,]

p <- ggplot() 
# store plot of regression line as well as 90% credible interval for regression mean
p2 <- p +
  geom_point(data = d2,
             aes(logarea, price), shape = 1, color = 'dodgerblue') +
  geom_ribbon(aes(area_new, ymin = hpdi_l, ymax = hpdi_u),
              alpha = .1) +
  geom_abline(data = post,
              aes(intercept = mean(alpha), slope = mean(beta))) +
  labs(subtitle="HPDI Interval = 0.9")

y_pi <- sapply(area_new,
               function(x) rnorm(NROW(post), post$alpha + post$beta * x, post$sigma)
)

# calculate 90% credible intervals for prediction intervals
y_phdi = HDInterval::hdi(y_pi, credMass=0.9)
pi_l = y_phdi[1,]
pi_u = y_phdi[2,]

# plot prediction intervals along with regression line & credible interval for linear mean
p2 + geom_ribbon(mapping = aes(area_new, ymin=pi_l, ymax=pi_u), alpha = 0.05) +
  labs(subtitle = 'Prediction Intervals = 0.9')

# Q3
# apply f_mu function to log of 75
mu_75 <- sapply(log(75), f_mu)
# calculate highest posterior density interval for the prediction
y_hdi_75 = HDInterval::hdi(mu_75, credMass=0.9)
# calculate mean of the predictions as point estimate of our prediction
mean(mu_75)
# show prediction interval
pi_75 = y_hdi_75[,1]
pi_75

# repeat same analysis for area of 175 m^2
mu_175 <- sapply(log(175), f_mu)
y_hdi_175 = HDInterval::hdi(mu_175, credMass=0.9)

pi_175 = y_hdi_175[,1]
mean(mu_175)
pi_175
