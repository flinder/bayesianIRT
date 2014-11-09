## M-H sampler for Rasch model
set.seed(36832)

# Simulated data
J <- 100 # Number of subjects
K <- 500 # Number of items
N <- J * K # Number of observations
theta.t <- rnorm(J, 0, 1) # Ability scores
beta.t <- c(0, rnorm(K-1, 0, 1))


# Draw data
ip <- function(beta, theta){
  1 / (1 + exp(beta - theta))
} 

ps <- sapply(beta.t, ip, theta = theta.t)
Y <- apply(ps, 2, rbinom, n = J, size = 1)

# logl
logl <- function(beta, theta) -log(1 + exp(beta - theta))

# Log Likelihood
loglik <- function(Y, theta, beta) sum(sapply(beta , logl, theta = theta))

# Log Prior
logprior <- function(theta, beta) {
  sum(dnorm(theta, 0, 1, log = T)) + sum(dnorm(beta, 0, 1, log = T))
}

# Proposals 
logproposal=function(theta.new, theta, beta.new, beta.2, psd){
  sum(dnorm(theta.new, theta, psd, log = T) + dnorm(beta.new, beta, psd, log = T))
}

sample.theta=function(theta, psd) sapply(theta, rnorm, n = 1,  sd = psd)
sample.beta=function(beta, psd) sapply(beta, rnorm, n = 1,  sd = psd)


# starting values
theta <- rep(1, J)
beta <- rep(1, K)

# Proposal standard deviation
psd <- 0.7

##
## MCMC loop
##
n.mcmc <- 1000
res <- matrix(NA, nrow = n.mcmc, ncol = (J + K))
colnames(res) <- c(paste0("theta[", c(1:J), "]"), paste0("beta[", c(1:K), "]"))
accept <- 0

for(k in 1:n.mcmc){
  ## print out iteration
  cat(k, "\r")
  ##
  ## sample theta
  ##
  
  theta.star <- sample.theta(theta, psd)
  beta.star <- sample.beta(beta, psd)
    
  ##
  ## calculate MH acceptance probability
  ##
  
  
  mh1 <- loglik(Y, theta.star, beta.star) + logprior(theta.star, beta.star) + 
    logproposal(theta, theta.star, beta, beta.star, psd)
  
  mh2 <- loglik(Y, theta, beta) + logprior(theta, beta) + 
    logproposal(theta.star, theta, beta.star, beta, psd)
  
  AP <- exp(mh1 - mh2)
  
  ##
  ## accept / reject using MH step
  ##
  
  u <- runif(1)
  if(u<AP){
    accept <- accept+1
    theta <- theta.star
    beta <- beta.star
  }
  
  ##
  ## save parameters
  ##
  
  res[k, 1:J] <- theta
  res[k, (J + 1):(J + K)] <- beta
}

## look at acceptance ratio (try to keep between .2 to .4
accept/n.mcmc

# Plot single chains
plot_chain <- function(chain) plot(c(1:length(chain)), chain, type = "l")
plot_chain(res[, 600])
