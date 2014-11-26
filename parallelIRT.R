rm(list=ls())
#################################################
# Comparison of JAGS and STAN for 2pl IRT models
#################################################

library(rjags)
library(rstan)
library(ggplot2)
library(truncnorm)
library(doParallel)

set.seed(2377344)
set_cppo("fast")

setwd("C:/Users/samsung/Dropbox/bd_irt/bayesianIRT")
#setwd("C:/Users/flinder/Dropbox/bd_irt/bayesianIRT")

# Load plot function
source("plot_pmeans.R")

### Generate some data according to:
# Y[j, k] ~ Bern(pi[j, k])
# pi[j, k] = 1/(1 + exp(-alpha[k] * (theta[j] - beta[k])))

J <- 200 # Number of subjects
K <- 10 # Number of items
N <- J * K # Number of observations
theta <- rnorm(J, 0, 1) # Ability scores

## True item-parameter matrix
# 1 alpha:Discrimination # 2 beta:Difficulty 
tpar <- rbind(rlnorm(K, -1, 0.35), rnorm(K, 0, 1))
rownames(tpar) <- c("alpha", "beta")

# Draw data
ip <- function(tpar, theta){
  1 / (1 + exp( - tpar[1] * (theta - tpar[2])))
} 
ps <- apply(tpar, 2, ip, theta = theta)
Y <- apply(ps, 2, rbinom, n = J, size = 1)

### Parallel
# Function run mcmc on single partition of the data
fit <- function(data){
  data.jags <- list("J" = nrow(data), "K" = ncol(data), "Y" = data)
  par.tosave <- c("alpha", "beta", "theta")
  mod <- jags.model(file = "models/jags_2pl_irt.txt", data = data.jags, 
                    n.chains = 1, n.adapt = 1000)
  smpls <- coda.samples(model = mod, variable.names = par.tosave, 
                        n.iter = 5000, thin = 1)
  smpls
}

# Function to partition the data and run mcmc on each partition
# and combines them according to neiswanger2014asymptotically
parMCMC <- function(data, nPart){
  CORES <- detectCores()
  cl <- makeCluster(CORES)
  registerDoParallel(cl)
  dat <- suppressWarnings(cbind(data, rep(1:nPart, each = floor(nrow(data) / nPart))))
  pdat <- split(data.frame(dat[, -ncol(dat)]), f = dat[, ncol(dat)])
  reform <- function(df){
    df <- as.matrix(df)
    colnames(df) <- NULL
    df
  }
  pdat <- lapply(pdat, reform)
  mcmc.out <- foreach(i = pdat, .packages = c("rjags"), .export = "fit") %dopar% fit(i)
  stopCluster(cl)
  mcmc.out
}

t.par <- system.time(par.mcmc <- parMCMC(Y, 4))


# Sequential
## JAGS
data.jags <- list("J" = J, "K" = K, "Y" = Y)
par.tosave <- c("alpha", "beta", "theta")

f <- function(){
  jmod <- jags.model(file = "models/jags_2pl_irt.txt", data = data.jags, 
                     n.chains = 1, n.adapt = 1000)
  mcmcres <- coda.samples(model = jmod, variable.names = par.tosave, 
                          n.iter = 5000, thin = 1)
  list(jmod, mcmcres)
}
t.seq <- system.time(jags.res <- f())
