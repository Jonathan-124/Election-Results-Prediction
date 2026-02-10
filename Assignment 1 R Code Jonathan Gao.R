# 1a)
data <- read.delim("Mconc.txt",header = F)
data <- as.numeric(data[,1])
hist(data, main="Observed Glass", xlab="M Concentration")

#1b)
weights <- c(0.15, 0.3, 0.4, 0.15)
means <- c(32.4, 35.2, 38.9, 42.1)
vars <- c(9, 1, 1, 9)
prior_mixture <- function(x) {sum(weights*dnorm(x,means,sqrt(vars)))}
prior_mixture <- Vectorize(prior_mixture)

# Plotting the mixture normal
curve(prior_mixture(x), from=min(data), to=max(data))

#Plotting histogram with mixture
hist(data, main="Observed Glass Frequency vs. Mixture Prior", xlab="M Concentration", freq=FALSE, ylim = c(0, 0.18))
x_vals <- seq(min(data), max(data), length = 1000)
lines(x_vals, prior_mixture(x_vals), col = 2, lwd = 2)

#1c) 
# Finding the prior predictive function
prior_predictive_vars <- vars + c(0.09, 0.09, 0.09, 0.09)
prior_predictive <- function(x) {sum(weights*dnorm(x, means, sqrt(prior_predictive_vars)))}
prior_predictive <- Vectorize(prior_predictive)
lines(x_vals, prior_predictive(x_vals), col = 4, lwd = 2)
curve(prior_predictive, from=20, to=60)

#1d)
tau_0 <- 1/(0.3^2)
taus <- 1/vars

# Function takes 4 inputs - y, the observed data, weights, means, and standard deviations
# Calculates the posterior weights
posterior_weights <- function(y, weights, means, vars){
  
  # Calculating the numerator - reweighing the four densities in the mixture
  new_weights <- weights*dnorm(y, means, sqrt(vars + c(0.09, 0.09, 0.09, 0.09)))
  
  # Rescaling to one
  new_weights <- new_weights / sum(new_weights)
  return (new_weights)
}

# Calculates the posterior means
posterior_means <- function (y, means, taus){
  new_means <- (y*tau_0 + taus*means)/(tau_0+taus)
  return(new_means)
}

# Calculates the posterior vars
posterior_vars <- function(y, taus){
  new_vars <- 1/(tau_0 + taus)
  return(new_vars)
}

posterior_weights <- posterior_weights(33.7, weights, means, vars)
posterior_means <- posterior_means(33.7, means, taus)
posterior_vars <- posterior_vars(33.7, vars)

posterior_mixture <- function(x) {sum(posterior_weights*dnorm(x,posterior_means,sqrt(posterior_vars)))}
posterior_mixture <- Vectorize(posterior_mixture)
curve(posterior_mixture, from=32.7, to=34.9, main= "Posterior Mixture", ylab = "Density")
abline(v=33.7, col =2, lwd = 2)

#1e)
posterior_predictive_vars <- posterior_vars + c(0.09, 0.09, 0.09, 0.09)
posterior_predictive <- function(x) {sum(posterior_weights*dnorm(x, posterior_means, sqrt(posterior_predictive_vars)))}
posterior_predictive <- Vectorize(posterior_predictive)
# lines(x_vals, prior_predictive(x_vals), col = 4, lwd = 2)
curve(posterior_predictive, from=32.7, to=34.9, main= "Posterior Predictive", ylab = "Density")

#1f)
posterior_predictive(33.9)
prior_predictive(33.9)

likelihood_ratio <- posterior_predictive(33.9) / prior_predictive(33.9)

#2a)
data_2 <- read.table("bioProcess.txt", header=TRUE)
plot(data_2, main="Bio Process Scatterplot", pch=20)

#2b)
# Vectors of the x and y variables
bio_x <- data_2[,1]
bio_y <- data_2[,2]

# Computes log-likelihood using bio_x data given input params beta and gamma 
log_bio_likelihood <- function(beta, gamma){
  sum(log(dgamma(bio_y, 25, 25/(beta * bio_x + gamma * log(bio_x + 1)))))
}

# Vectorizes the likelihood functions so vectors of beta/gamma can be inputted
log_bio_likelihood <- Vectorize(log_bio_likelihood)

#2c)

# Calculates posterior given beta and gamma values
bio_posterior <- function(beta, gamma){
  # exp(log_bio_likelihood(beta, gamma)) / (beta * gamma)
  exp(log_bio_likelihood(beta, gamma) - log(beta) - log(gamma))
}

# Vectorizes the posterior
bio_posterior <- Vectorize(bio_posterior)

# Setting up bounds for beta and gamma values and number of discretizations
beta_0 <- 0.3
beta_1 <- 1.3
gamma_0 <- 4
gamma_1 <- 7
discretization <- 100

# Computes vectors of beta, gamma using the given bounds
beta_range <- seq(beta_0, beta_1, length.out = discretization)
gamma_range <- seq(gamma_0, gamma_1, length.out = discretization)

# Outer product of beta gamma pairs, and applies bio_posterior to get a matrix
# of posterior values corresponding to the coordinates on the grid
bio_posterior_grid <- outer(beta_range, gamma_range, bio_posterior)

# Improves numerics? From lecture 3 code
bio_posterior_grid <- bio_posterior_grid - max(bio_posterior_grid) 
  
# Re-scaling the posterior to sum to 1
bio_posterior_grid <- bio_posterior_grid/sum(bio_posterior_grid) 
  
# Heat map
image(beta_range, gamma_range, bio_posterior_grid, main="Posterior Distribution", xlab="Beta", ylab="Gamma")

#2e)
pred <- function(beta, gamma){
  1-pgamma(20, 25, 25/(10*beta + log(11)*gamma))
}
sum(outer(beta_range,gamma_range,Vectorize(pred))*bio_posterior_grid)

#2f)

# New log posterior - theta is a vector c(beta, gamma, alpha)
log_bio_posterior_new <- function(theta){
  beta <- theta[1]
  gamma <- theta[2]
  alpha <- theta[3]
  log_likelihood <- sum(log(dgamma(bio_y, alpha, alpha/(beta * bio_x + gamma * log(bio_x + 1)))))
  logpost <- log_likelihood - log(beta) - log(gamma) - log(alpha)
  return(logpost)
}
log_bio_posterior_new

mode_guess <- c(0.7, 5.5, 25)

library(LearnBayes)
new_result <- laplace(log_bio_posterior_new, mode_guess)
print(new_result)

