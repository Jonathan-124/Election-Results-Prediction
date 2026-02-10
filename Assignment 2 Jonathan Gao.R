#1a) 
data <- read.table("election25.txt",header=TRUE)
q <- data$ResultPreviousElection
y <- data[!is.na(data$exitPollResults), ]$exitPollResults
p <- data$finalResults

# Previous results graph
prev_proportion <- exp(q)/(1+exp(q))
hist(prev_proportion, main="Previous Results", xlab="Proportion of Votes won by party A")
num_prev_won <- sum(prev_proportion > 1/2)

# Exit poll results graph
hist(y/100, main="Exit Poll Results", xlab="Proportion of Votes won by party A")
exit_polls_won <- sum(y > 50)

# Final results graph
final_proportion <- exp(p)/(1+exp(p))
hist(prev_proportion, main="Final Results", xlab="Proportion of Votes won by party A")
num_won_final <- sum(final_proportion > 1/2)

# Checking the mu/gamma values using post-hoc results 
ds <- data$finalResults - data$ResultPreviousElection
hist(ds, main="d distribution", xlab=("d=p-q"), freq = TRUE)
mu_check <- mean(ds)
gamma_check <- -log((sd(ds))^2)

#1d)

# Hamiltonean Monte Carlo Function with inputs data and time t
HMC <- function(data, t){
  
  # Filtering out the omegas and psis at time t
  omega_set <- data[data$timeReceived <= t, ]
  psi_set <- data[data$timeReceived > t & !is.na(data$exitPollResults), ]
  s <- nrow(omega_set)
  r <- nrow(psi_set)

  # Extracting the vectors for omegas, zs, ws
  omegas <- omega_set$finalResults - omega_set$ResultPreviousElection
  zs <- psi_set$exitPollResults
  ws <- psi_set$ResultPreviousElection
  
  U <- function(q){
    # q is a vector - c(mu, gamma, psi)
    # mu and gamma are numbers;
    # psi is a vector of length r
    mu <- q[1]
    gamma <- q[2]
    
    if (length(q) >2){
      # If there are still districts with exit polls but no final results
      psi <- q[3:length(q)]
    
      # Implementing the formula for the negative logarithm posterior density
      (exp(gamma)/2)*(sum((omegas-mu)^2)+sum((psi-mu)^2)) - gamma*(s+r)/2 - sum(zs*psi - 100*log(1+exp(psi+ws)))
    } else{
      
      # If there are no more districts with exit polls
      (exp(gamma)/2)*(sum((omegas-mu)^2)) - gamma*(s)/2 
    }
  }
  
  gradU <- function(q){
    # q is a vector - c(mu, gamma, psi)
    # mu and gamma are numbers;
    # psi is a vector of length r
    mu <- q[1]
    gamma <- q[2]
    
    if (length(q) >2){
      # If there are still districts with exit polls but no final results
      psi <- q[3:length(q)]
    
      # Finding the gradients 
      grad_1 <- -exp(gamma)*(sum(omegas - mu) + sum(psi - mu))
      grad_2 <- (exp(gamma)/2)*(sum((omegas - mu)^2) + sum((psi - mu)^2)) - (s+r)/2
      grad_psi <- exp(gamma)*(psi - mu) - zs + (100*exp(psi + ws)/(1+exp(psi+ws)))
      
      # Returns gradients
      c(grad_1, grad_2, grad_psi)
    } else{
      # If there are no more districts with exit polls
      grad_1 <- -exp(gamma)*sum(omegas - mu)
      grad_2 <- (exp(gamma)/2)*sum((omegas - mu)^2)  - (s+r)/2
      
      # Returns gradients
      c(grad_1, grad_2)
    }
  }
  
  # Leapfrog algorithm - code from Lecture 10 R code
  # currentp  the current value of p
  # currentq  the current value of q
  # gradU     the gradient function for the U function 
  # epsilon   the size of the steps of the leapfrog
  # L         the number of the steps for the leapfrog
  # Output: 
  # newp      the new value for p
  # newq      the new value for q
  leapfrog <- function(currentp, currentq, gradU, epsilon, L) {
    q <- currentq 
    p <- currentp - epsilon*gradU(q)/2 # Start with a step of half length
    for (i in 1:L) {
      q <- q + epsilon*p
      if (i!=L) p <- p - epsilon*gradU(q)
    }
    p <- p - epsilon*gradU(q)/2 # End with half step
    list(newp=p, newq=q)
  }
  
  # Hamiltonean MCMC - code from lecture 10 R code
  N     <- 500
  epsilon <- 0.005
  L     <- 50
  trace <- matrix(0, N, r+2)
  # Set q0: 
  # At q0, set mu_0=0, gamma_0=3, psi_1...psi_r=0 (r constituencies have exit polls but no final results)
  trace[1,] <- c(0, 3, rep(0, r))
  
  for (i in 2:N) {
    currentp <- rnorm(r+2, 0, 1)
    # To compute the leapfrog result: 
    res <- leapfrog(currentp, trace[i-1,], gradU, epsilon, L) 
    if (runif(1)<exp(- U(res$newq) - sum(res$newp^2)/2 + U(trace[i-1,]) + sum(currentp^2)/2)){
      trace[i,] <- res$newq
    }
    else
      trace[i,] <- trace[i-1,]
  }
  
  # Returns all iterations of HMC as a matrix
  return(trace)
}

# Plots trace plots when iterations of HMC algorithm are returned
x <- HMC(data, 1.5)
mu_mean <- mean(x[,1])
gamma_mean <- mean(x[,2])

plot(x[,1], type="l", main="Mu Trace Plot", xlab="Iteration", ylab="Mu")
plot(x[,2], type="l", main="Gamma Trace Plot", xlab="Iteration", ylab="Gamma")
plot(x[,3], type="l", main="Psi_1 Trace Plot", xlab="Iteration", ylab="Psi_1")
smoothScatter(x[,1], x[,2], main="Scatterplot for mu, gamma HMC", xlab="Mu", ylab="Gamma")

#1e)

# The first part of this function is basically the same as the previous example
# However instead of returning the simulated data as a matrix, this function uses
# Monte Carlo Integration to return the probability that party A will win given data at time t

A_wins <- function(data, t){
  
  # Filtering out the omegas and psis at time t
  omega_set <- data[data$timeReceived <= t, ]
  psi_set <- data[data$timeReceived > t & !is.na(data$exitPollResults), ]
  s <- nrow(omega_set)
  r <- nrow(psi_set)

  # Extracting the vectors for omegas, zs, ws
  omegas <- omega_set$finalResults - omega_set$ResultPreviousElection
  zs <- psi_set$exitPollResults
  ws <- psi_set$ResultPreviousElection
  
  U <- function(q){
    # q is a vector - c(mu, gamma, psi)
    # mu and gamma are numbers;
    # psi is a vector of length r
    mu <- q[1]
    gamma <- q[2]
    
    if (length(q) >2){
      # If there are still districts with exit polls but no final results
      psi <- q[3:length(q)]
      
      # Implementing the formula for the negative logarithm posterior density
      (exp(gamma)/2)*(sum((omegas-mu)^2)+sum((psi-mu)^2)) - gamma*(s+r)/2 - sum(zs*psi - 100*log(1+exp(psi+ws)))
    } else{
      
      # If there are no more districts with exit polls
      (exp(gamma))/2*(sum((omegas-mu)^2)) - gamma*s/2 
    }
  }
  
  gradU <- function(q){
    # q is a vector - c(mu, gamma, psi)
    # mu and gamma are numbers;
    # psi is a vector of length r
    mu <- q[1]
    gamma <- q[2]
    
    if (length(q) >2){
      # If there are still districts with exit polls but no final results
      psi <- q[3:length(q)]
      
      # Finding the gradients 
      grad_1 <- -exp(gamma)*(sum(omegas - mu) + sum(psi - mu))
      grad_2 <- (exp(gamma)/2)*(sum((omegas - mu)^2) + sum((psi - mu)^2)) - (s+r)/2
      grad_psi <- exp(gamma)*(psi - mu) - zs + 100*exp(psi + ws)/(1+exp(psi+ws))
      
      # Returns gradients
      c(grad_1, grad_2, grad_psi)
    } else{
      # If there are no more districts with exit polls
      grad_1 <- (-exp(gamma))*sum(omegas - mu)
      grad_2 <- (exp(gamma))/2*sum((omegas - mu)^2)  - s/2
      
      # Returns gradients
      c(grad_1, grad_2)
    }
  }
  
  # Leapfrog algorithm - code from Lecture 10 R code
  # currentp  the current value of p
  # currentq  the current value of q
  # gradU     the gradient function for the U function 
  # epsilon   the size of the steps of the leapfrog
  # L         the number of the steps for the leapfrog
  # Output: 
  # newp      the new value for p
  # newq      the new value for q
  leapfrog <- function(currentp, currentq, gradU, epsilon, L) {
    q <- currentq 
    p <- currentp - epsilon*gradU(q)/2 # Start with a step of half length
    for (i in 1:L) {
      q <- q + epsilon*p
      if (i!=L) p <- p - epsilon*gradU(q)
    }
    p <- p - epsilon*gradU(q)/2 # End with half step
    list(newp=p, newq=q)
  }
  
  # Hamiltonean MCMC - code from lecture 10 R code
  N     <- 500
  epsilon <- 0.005
  L     <- 50
  trace <- matrix(0, N, r+2)
  # Set q0: 
  # At q0, set mu_0=0, gamma_0=3, psi_1...psi_r=0 (r constituents have results)
  trace[1,] <- c(0, 3, rep(0, r))
  
  for (i in 2:N) {
    currentp <-  rnorm(r+2, 0, 1) 
    # To compute the leapfrog result: 
    res <- leapfrog(currentp, trace[i-1,], gradU, epsilon, L) 
    if (runif(1)<exp(- U(res$newq) - sum(res$newp^2)/2 + U(trace[i-1,]) + sum(currentp^2)/2)){
      trace[i,] <- res$newq
    }
    else
      trace[i,] <- trace[i-1,]
  }
  
  # 1e) new stuff below
  # Starting from the 2nd row of HMC simulation (as the initial conditions are prescribed)
  # At time t, Computes the probability that party A will win the election 
  # At each iteration, simulates results for each district using each iteration posterior
  # A log-odds ratio >0 indicates party A won in each district
  #.# For each of the omegas, we already know the results so no need to simulate
  # # For each of the thetas, we use the mu and gamma to simulate results in each district as no other info known
  # # For each of the psis, we use the iteration's posterior directly
  
  # Calculates number of seats already won at time t
  num_won <- nrow(data[data$timeReceived <= t & data$finalResults > 0,])

  # Calculates the Monte Carlo Integral of all psi variables
  # For each iteration, determines the number of psi districts that party A wins
  # i.e. has log-odds ratio >0
  # pred_psi_wins is thus a vector where each entry is an int with [0, r]
  
  pred_psi_wins <- rep(0, N-1) # Vector of length N-1 (for each iteration, counts number of predicted wins in r districts with only exit polls)

  if(r>0){
    iterated_psis <- trace[2:N, 3:(r+2)] 
    # Matrix of N-1 x r columns (N-1 iterations, r districts with only exit polls)

    iterated_psis <- iterated_psis + matrix(rep(ws, N-1), ncol=r, byrow=TRUE)

    # Add the ws (rx1) vector to each row (i.e. each psi iteration has the previous year's election result added to it)
    # This is from the formula p_i = d_i + q_i, where we have now estimated the p_i, log-odds of winning this election
    if (r > 1){
      pred_psi_wins <- rowSums(iterated_psis > 0)
      # Counts the number of predicted wins in districts with exit polls for each iteration
    }
  }

  # Calculates the Monte Carlo Integral of all theta variables
  # In each step, the simulated mu and gamma is used
  # For each iteration, determines the number of theta districts that party A wins
  # i.e. has log-odds ratio >0
  # pred_omega_wins is thus a vector where each entry is an int with [0, 231-r-s]
  
  pred_theta_wins <- rep(0, N-1) # Vector of length N-1 (for each iteration, counts number of wins in 231-r-s districts with no data)
  
  if (231-r-s > 0){
    theta_results <- matrix(0, N-1, 231-r-s)
    # Matrix of (N-1) rows x (231-r-s) columns (N-1 iterations, 231-r-s districts with no results)
    
    prev_theta_results <- data[data$timeReceived > t & is.na(data$exitPollResults), ]$ResultPreviousElection
    # Vector of length (231-r-s) representing log_odds results from the previous election

    for (i in 2:N) {
      iterated_mu <- trace[i, 1]
      iterated_gamma <- trace[i, 2]
      theta_results[i-1,] <- rnorm(231-r-s, iterated_mu, exp(-iterated_gamma))
      # Simulates theta for each district with no data whatsoever for N-1 steps
      
    }
    theta_results <- theta_results + matrix(rep(prev_theta_results, N-1), ncol=231-r-s, byrow=TRUE)
    # Add the previous theta results (N-r-s x1) vector to each row (i.e. each theta iteration has the previous year's election result added to it)
    # This is from the formula p_i = d_i + q_i, where we have now estimated the p_i, log-odds of winning this election
    
    if (r+s <230){
      pred_theta_wins <- rowSums(theta_results > 0)
      # Counts the number of wins in districts with no data for each iteration
    }
  }
  
  # Finds the total predicted number of districts won for each iteration
  districts_won <- pred_psi_wins + pred_theta_wins + num_won
  # Returns proportion of iterations that have party A winning majority of votes
  sum(districts_won > 231/2) / (N-1)
}

A_wins(data, 3.4)

# 1f)
# Lists the unique times that new results come in 
update_times <- sort(unique(data$timeReceived))

# initializes proportion of determined seats 
determined_seats <- rep(0, length(update_times))
determined_seats_for_A <- rep(0, length(update_times))
prob_A_wins <- rep(0, length(update_times))

# Loops over the unique update times
# Finds the number of determined seats at time t, and also the prob. that A wins
for (i in 1:length(update_times)){
  # probably can optimize this a bit more...
  # Counts the number of determined seats at time t 
  determined_seats[i] <- nrow(data[data$timeReceived <= update_times[i], ])
  determined_seats_for_A[i] <- nrow(data[data$timeReceived <= update_times[i] & data$finalResults > 0, ]) 
  prob_A_wins[i] <- A_wins(data, update_times[i])
}

determined_seats_for_A_proportion <- determined_seats_for_A / determined_seats
determined_seats_proportion <- determined_seats / 231


# Changing margins to plot outside of graph region
par(mar=c(10.1, 4.1, 3.1, 2.1), xpd=TRUE)
plot(update_times, determined_seats_for_A_proportion, type="l", col="black", xlab="Time", ylab="Proportion", main="Election Results", ylim=c(0, 1))
lines(update_times, prob_A_wins, type="l", col="red")
lines(update_times, rep(0.5, length(update_times)), lty=2, col="blue")
legend("bottom", inset=c(1,-1), legend = c("Probability Party A wins", "Proportion of Determined Seats won by Party A", "Majority Threshold"), col = c("red", "black", "blue"), lty = 1)

######################
# Some bonus graphs...
plot(update_times, determined_seats_proportion, type="l", col="red", xlab="Time", ylab="Proportion", main="Election Results", ylim=c(0, 1))
lines(update_times, determined_seats_for_A_proportion, type="l", col="black")
lines(update_times, rep(0.5, length(update_times)), lty=2, col="blue")
legend("bottom", inset=c(1,-1), legend = c("Proportion of Final Results Released", "Proportion of Determined Seats won by Party A", "Majority Threshold"), col = c("red", "black", "blue"), lty = 1)

