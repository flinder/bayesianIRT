###############################################
# Comparison of JAGS and STAN for 3pl IRT model
###############################################
rm(list=ls())
library(rjags)
library(rstan)

set.seed(23773)

setwd("C:/Users/samsung/Dropbox/bd_irt/bayesianIRT")

### Generate some data according to:
# Y[j, k] ~ Bern(pi[j, k])
# pi[j, k] = gamma[k] + (1 - gamma[k])/(1 + exp(-alpha[k] * (theta[j] - beta[k])))
 
# Number of subjects
J <- 1000
# Number of items
K <- 50
# Number of observations
N <- J * K

# Ability scores
theta <- rnorm(J, 0, 1) 

## Parameter matrix
# gamma:Guessing # alpha:Discrimination # beta:Difficulty 
sigma <- rbind(runif(K), rlnorm(K, 0.9, 0.5), rnorm(K, 0.2, 0.1))

# Draw data
ip <- function(sigma, theta){
  sigma[1] + (1 - sigma[1]) / (1 + exp( - sigma[2] * (theta - sigma[3])))
} 
ps <- apply(sigma, 2, ip, theta = theta)
Y <- apply(ps, 2, rbinom, n = J, size = 1)


## JAGS
data.jags <- list("J" = J, "K" = K, "Y" = Y)
inits <- list(list("alpha" = rep(1.1, K), 
                   "beta" = rnorm(K, 0, 1.5), 
                   "gamma" = rep(0.2, K), 
                   "theta" = rnorm(J)
                   )
              )
par.tosave <- c("alpha", "beta", "gamma", "theta")

f <- function(){
  jmod <- jags.model(file = "models/jags_3pl_irt.txt", data = data.jags, 
                     inits = inits, n.chains = 1, n.adapt = 1e3)
  mcmcres <- coda.samples(model = jmod, variable.names = par.tosave, 
                          n.iter = 5e3, thin = 1)
  list(jmod, mcmcres)
}
t.jags <- system.time(jags.res <- f())
post.m.jags <- apply(jags.res[[2]][1][[1]], 2, mean)
comp.jags <- cbind(c(sigma[2, ], sigma[3, ], sigma[1, ], theta), post.m.jags)
plot(comp.jags)

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
  res.stan <- stan(model_code = model, model_name = "stan_3pl", data = data.stan, 
                   iter = 5e3, warmup = 1e3, chains = 1, verbose = TRUE)
  )

post.m.stan <- sapply(fit@sim$samples[[1]], mean)
comp.stan <- cbind(c(sigma[2, ], sigma[3, ], sigma[1, ], theta), post.m.stan[-(J+3*K+1)])
plot(comp.stan)
