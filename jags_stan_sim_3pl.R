#################################################
# Comparison of JAGS and STAN for 3pl IRT models
#################################################
rm(list=ls())
library(rjags)
library(rstan)
library(ggplot2)

set.seed(3439109)

#setwd("C:/Users/samsung/Dropbox/bd_irt/bayesianIRT")
setwd("C:/Users/flinder/Dropbox/bd_irt/bayesianIRT")

# Load plot function
source("plot_pmeans.R")

### Generate some data according to:
# Y[j, k] ~ Bern(pi[j, k])
# pi[j, k] = gamma[k] + (1 - gamma[k])/(1 + exp(-alpha[k] * (theta[j] - beta[k])))
 

J <- 100 # Number of subjects
K <- 10 # Number of items
N <- J * K # Number of observations
theta <- rnorm(J, 0, 1) # Ability scores

## True item-parameter matrix
# 1 alpha:Discrimination # 2  beta:Difficulty # 3 gamma:Guessing 
tpar <- rbind(rlnorm(K, -0.3, 0.3), rnorm(K, 0, 1), rbeta(K, 3, 8))
rownames(tpar) <- c("alpha", "beta", "gamma")

# Draw data
ip <- function(tpar, theta){
  tpar[3] + (1 - tpar[3]) / (1 + exp( - tpar[1] * (theta - tpar[2])))
} 
ps <- apply(tpar, 2, ip, theta = theta)
Y <- apply(ps, 2, rbinom, n = J, size = 1)

## JAGS
data.jags <- list("J" = J, "K" = K, "Y" = Y, "a.gamma" = 5, "b.gamma" = 17, 
                  "v.alpha" = 2, "v.beta" = 2)

# Starting values
inits <- list(list("alpha" = rep(1, K), # For chain 1
                   "beta" = rnorm(K, 0, 1.5), 
                   "gamma" = rep(0.1, K), 
                   "theta" = rnorm(J)
                   )
              )
par.tosave <- c("alpha", "beta", "gamma", "theta")

f <- function(){
  jmod <- jags.model(file = "models/jags_3pl_irt.txt", data = data.jags, 
                     inits = inits, n.chains = 1, n.adapt = 1000)
  mcmcres <- coda.samples(model = jmod, variable.names = par.tosave, 
                          n.iter = 5000, thin = 1)
  list(jmod, mcmcres)
}
t.jags <- system.time(jags.res <- f())


# Plot results
jagspost1 <- jags.res[[2]][1][[1]]


plot_pmeans(jagspost1, K, J, tpar, theta, "jags_3pl")
ggsave('plots/jags_3pl_big.png')

## STAN
fileName <- "models/stan_3pl_irt.txt"
model <- readChar(fileName, file.info(fileName)$size)
Y.vec <- stack(as.data.frame(Y))[, 1]
data.stan <- list("J" = J, 
                  "K" = K, 
                  "N" = N,  
                  "jj" = rep(c(1:J), K), 
                  "kk" = rep(c(1:K), each = J), 
                  "Y" = Y.vec
                  )
t.stan <- system.time(
  stan.res <- stan(model_code = model, model_name = "stan_3pl", data = data.stan, 
                   iter = 5000, warmup = 1000, chains = 1, verbose = TRUE)
  )


stanpost <- do.call(cbind,stan.res@sim$samples[[1]][- (J + 3 * K + 1)])

# Plot posterior means
plot_pmeans(stanpost, K, J, tpar, "stan_3pl")
ggsave('plots/stan_3pl_big.png')

##########################
# With hierarchical priors
##########################

## JAGS
data.jags <- list("J" = J, "K" = K, "Y" = Y, "sigma_0" = diag(2), "mu_0" = c(0,0))
f <- function(){
  jmod <- jags.model(file = "models/jags_3pl_h_irt.txt", data = data.jags, 
                     inits = inits, n.chains = 2, n.adapt = 500)
  mcmcres <- coda.samples(model = jmod, variable.names = par.tosave, 
                          n.iter = 2000, thin = 1)
  list(jmod, mcmcres)
}
t.jags <- system.time(jags.res.h <- f())
#save(jags.res, file = "jags_res.RData")
#load("jags_res.RData")