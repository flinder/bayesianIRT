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

setwd("C:/Users/samsung/Dropbox/bd_irt/bayesianIRT")
#setwd("C:/Users/flinder/Dropbox/bd_irt/bayesianIRT")

### Generate some data according to:
# Y[j, k] ~ Bern(pi[j, k])
# pi[j, k] = 1/(1 + exp(-alpha[k] * (theta[j] - beta[k])))

J <- 1000 # Number of subjects
K <- 50 # Number of items
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

### Parallel (for item parameters)
# Function run mcmc on single partition of the data
fit <- function(data, n.adapt, n.iter){
  data.jags <- list("J" = nrow(data), "K" = ncol(data), "Y" = data)
  par.tosave <- c("alpha", "beta", "theta")
  mod <- jags.model(file = "models/jags_2pl_irt.txt", data = data.jags, 
                    n.chains = 1, n.adapt = n.adapt)
  smpls <- coda.samples(model = mod, variable.names = par.tosave, 
                        n.iter = n.iter, thin = 1)
  smpls
}

# Function to partition the data and run mcmc on each partition
# and combines them according to neiswanger2014asymptotically
# for now just gives item parameters
parMCMC <- function(data, nPart, n.adapt = 1000, n.iter = 5000){
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
  
  mcmc.out <- foreach(i = pdat, .packages = c("rjags"), .export = "fit") %dopar% {
    fit(i, n.adapt = n.adapt, n.iter = n.iter)}
  
  posteriors <- lapply(mcmc.out, function(x) x[[1]])
  
  ## Posterior combination functions for item parameters
  # Parametric BvM
  bvm <- function(posteriors){
    posteriors <- lapply(posteriors, function(x) x[, !grepl("theta", colnames(x))])
    vcms <- foreach(i = posteriors) %dopar% var(i)
    ivcms <- foreach(i = vcms) %dopar% solve(i)
    var.c <- solve(Reduce("+", ivcms))
    means <- foreach(i = posteriors) %dopar% apply(i, 2, mean)
    mprod <- foreach(i = ivcms, j = means) %dopar% (i %*% j)
    mean.c <- var.c %*% Reduce("+", mprod)
    list("post.means" = mean.c, "post.variance" = var.c)
  }
  
  out <- list(combined = bvm(posteriors), separate = posteriors)
  stopCluster(cl)
  out
}

t.par <- system.time(par.mcmc <- parMCMC(Y, 4))


par.mcmc[[1]]

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


m1 <- par.mcmc$combined$post.means
a <- jags.res[[2]][[1]]
a <- a[, !grepl("theta", colnames(a))]
m2 <- apply(a, 2, mean)
pdat <- data.frame(m1, m2, grp = factor(rep(c("alpha", "beta"), each = 50)))
p <- ggplot(pdat, aes(m1, m2, col = grp))
p + geom_point()
