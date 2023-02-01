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
library("rethinking")

SEED <- 151515 # set random seed for reproducability

data <- read.csv("coffee.csv", header = TRUE)

# Initial analysis of data
str(data)
plot(data[,2], data[,1], xlab = "Cafe", ylab = 'Rating')

# Plot the density of ratings for each of the cafes
plot <- matrix(data = NA, nrow = nrow(data), ncol = 2)
plot$ratings <- data$ratings
plot$cafe <- as.character(data$cafe)

ggplot(as.data.frame(plot), aes(ratings, fill = cafe, colour = cafe)) +
  geom_density(alpha = 0.2) + theme_grey()+
  labs(x = 'rating', y = '', title = 'Distribution of ratings for each cafe') +
  scale_y_continuous(breaks = NULL) + xlim(0, 100)

# store the ratings associated with each cafe individually
ratings1 = data$ratings[which(data[,2]==1)]
ratings2 = data$ratings[which(data[,2]==2)]
ratings3 = data$ratings[which(data[,2]==3)]
ratings4 = data$ratings[which(data[,2]==4)]
ratings5 = data$ratings[which(data[,2]==5)]

writeLines(readLines("coffee.stan"))

cafe1 <- list(rating = ratings1, N = length(which(data[,2]==1)), cafe = data$cafe[which(data[,2]==1)])
fit_cafe1 <- stan("coffee.stan",  data=cafe1, iter=20000, refresh = 0)

print(fit_cafe1, probs = c(0.10, 0.5, 0.9))

cafe_sep_1 <- as.data.frame(fit_cafe1, pars=c("mu"))
p1 <- stack(as.data.frame(cafe_sep_1))
p1$ind = "cafe1"

cafe2 <- list(rating = ratings2, N = length(which(data[,2]==2)), cafe = data$cafe[which(data[,2]==2)])
fit_cafe2 <- stan("coffee.stan",  data=cafe2, iter=20000, refresh = 0)

print(fit_cafe2, probs = c(0.10, 0.5, 0.9))

cafe_sep_2 <- as.data.frame(fit_cafe2, pars=c("mu"))
p2 <- stack(as.data.frame(cafe_sep_2))
p2$ind = "cafe2"

cafe3 <- list(rating = ratings3, N = length(which(data[,2]==3)), cafe = data$cafe[which(data[,2]==3)])
fit_cafe3 <- stan("coffee.stan",  data=cafe3, iter=20000, refresh = 0)

print(fit_cafe3, probs = c(0.10, 0.5, 0.9))

cafe_sep_3 <- as.data.frame(fit_cafe3, pars=c("mu"))
p3 <- stack(as.data.frame(cafe_sep_3))
p3$ind = "cafe3"

cafe4 <- list(rating = ratings4, N = length(which(data[,2]==4)), cafe = data$cafe[which(data[,2]==4)])
fit_cafe4 <- stan("coffee.stan",  data=cafe4, iter=20000, refresh = 0)

print(fit_cafe4, probs = c(0.10, 0.5, 0.9))

cafe_sep_4 <- as.data.frame(fit_cafe4, pars=c("mu"))
p4 <- stack(as.data.frame(cafe_sep_4))
p4$ind = "cafe4"

cafe5 <- list(rating = ratings5, N = length(which(data[,2]==5)), cafe = data$cafe[which(data[,2]==5)])
fit_cafe5 <- stan("coffee.stan",  data=cafe5, iter=20000, refresh = 0)

print(fit_cafe5, probs = c(0.10, 0.5, 0.9))

cafe_sep_5 <- as.data.frame(fit_cafe5, pars=c("mu"))
p5 <- stack(as.data.frame(cafe_sep_5))
p5$ind = "cafe5"

# Plot posterior distributions of each of the individual coffee shops
posterior_sep <- rbind(p1,p2, p3, p4, p5)
ggplot(posterior_sep, aes(values, fill = ind, colour = ind)) +
  geom_density(alpha = 0.2) + theme_grey()+
  labs(x = 'rating', y = '', title = 'Seperate model - Distributions of posterior means') +
  scale_y_continuous(breaks = NULL) + xlim(0, 100)
  
y_rep1 <- as.matrix(fit_cafe1, pars=c("y_rep"))
ppc_dens_overlay(data[which(data[,2]==1),1], y_rep1[1:40,])+theme_classic()
y_rep2 <- as.matrix(fit_cafe2, pars=c("y_rep"))
ppc_dens_overlay(data[which(data[,2]==2),1], y_rep2[1:40,])+theme_classic()
y_rep3 <- as.matrix(fit_cafe3, pars=c("y_rep"))
ppc_dens_overlay(data[which(data[,2]==3),1], y_rep3[1:40,])+theme_classic()
y_rep4 <- as.matrix(fit_cafe4, pars=c("y_rep"))
ppc_dens_overlay(data[which(data[,2]==4),1], y_rep4[1:40,])+theme_classic()
y_rep5 <- as.matrix(fit_cafe5, pars=c("y_rep"))
ppc_dens_overlay(data[which(data[,2]==5),1], y_rep5[1:40,])+theme_classic()

# Pooled model
cafe_pooled <- list(rating = data$ratings, N = nrow(data), cafe = data$cafe)
fit_pooled <-stan("coffee.stan",  data=cafe_pooled, iter=5000, refresh = 0)

print(fit_pooled, probs = c(0.10, 0.5, 0.9))

pooled <- as.data.frame(fit_pooled, pars=c("mu"))
p6 <- stack(as.data.frame(pooled))

ggplot(p6, aes(values, fill = ind, colour = ind)) +
  geom_density(alpha = 0.2) + theme_grey()+
  labs(x = 'rating', y = '', title = 'Pooled Model - Distributions of posterior means') +
  scale_y_continuous(breaks = NULL) + xlim(0, 100)

y_rep_pooled <- as.matrix(fit_pooled, pars=c("y_rep"))
ppc_dens_overlay(data[,1], y_rep_pooled[1:40,])+theme_classic()

# Hierarchical Model
writeLines(readLines("cafe_hierarchical.stan"))

cafe_hier <- list(ratings1 = ratings1, ratings2 = ratings2, ratings3 = ratings3, ratings4 = ratings4, ratings5 = ratings5,
                  N1 = length(ratings1), N2 = length(ratings2), N3 = length(ratings3), N4 = length(ratings4), N5 = length(ratings5), cafe = data$cafe)
fit_hier <-stan("cafe_hierarchical.stan",  data=cafe_hier, iter=20000, refresh = 0)

print(fit_hier, probs = c(0.10, 0.5, 0.9))

post_hier <- as.data.frame(fit_hier, pars=c("mu"))
posterior_hier <- stack(as.data.frame(post_hier))

ggplot(posterior_hier, aes(values, fill = ind, colour = ind)) +
  geom_density(alpha = 0.2) + theme_grey()+
  labs(x = 'rating', y = '', title = 'Hierarchical model - Distributions of posterior means') +
  scale_y_continuous(breaks = NULL) + xlim(0, 100)

y_rep_hier1 <- as.matrix(fit_hier, pars=c("y_rep1"))
ppc_dens_overlay(data[which(data[,2]==1),1], y_rep_hier1[1:40,])+theme_classic()
y_rep_hier2 <- as.matrix(fit_hier, pars=c("y_rep2"))
ppc_dens_overlay(data[which(data[,2]==2),1], y_rep_hier2[1:40,])+theme_classic()
y_rep_hier3 <- as.matrix(fit_hier, pars=c("y_rep3"))
ppc_dens_overlay(data[which(data[,2]==3),1], y_rep_hier3[1:40,])+theme_classic()
y_rep_hier4 <- as.matrix(fit_hier, pars=c("y_rep4"))
ppc_dens_overlay(data[which(data[,2]==4),1], y_rep_hier4[1:40,])+theme_classic()
y_rep_hier5 <- as.matrix(fit_hier, pars=c("y_rep5"))
ppc_dens_overlay(data[which(data[,2]==5),1], y_rep_hier5[1:40,])+theme_classic()

# examine the effect using a hierarchical model has on posterior distributions
mean(data$ratings)
var(data$ratings)/nrow(data)

c(mean(ratings1), mean(ratings2), mean(ratings3), mean(ratings4), mean(ratings5))
c(var(ratings1)/length(ratings1), var(ratings2)/length(ratings2), var(ratings3)/length(ratings3), var(ratings4)/length(ratings4), var(ratings5)/length(ratings5))

post_sep <- cbind(cafe_sep_1, cafe_sep_2, cafe_sep_3, cafe_sep_4, cafe_sep_5)
colnames(post_sep) <- c("mu1", "mu2", "mu3", "mu4", "mu5")
# examine posterior means for the average rating for each cafe under each model
colMeans(post_sep)
mean(pooled$mu)
colMeans(post_hier)
# compare these figures to the overall mean
abs(colMeans(post_sep)-mean(data$ratings))
abs(mean(pooled$mu)-mean(data$ratings))
abs(colMeans(post_hier)-mean(data$ratings))
# examine the variance of estimates under each model
apply(post_sep, 2, var)
var(pooled$mu)
apply(post_hier, 2, var)