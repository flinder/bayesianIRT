#################################################
# Comparison of JAGS and STAN for 3pl IRT models
#################################################
rm(list=ls())
library(rjags)
library(rstan)
library(ggplot2)

set.seed(3439109)

setwd("C:/Users/samsung/Dropbox/bd_irt/bayesianIRT")
#setwd("C:/Users/flinder/Dropbox/bd_irt/bayesianIRT")

### Generate some data according to:
# Y[j, k] ~ Bern(pi[j, k])
# pi[j, k] = gamma[k] + (1 - gamma[k])/(1 + exp(-alpha[k] * (theta[j] - beta[k])))
 

J <- 300 # Number of subjects
K <- 20 # Number of items
N <- J * K # Number of observations
theta <- rnorm(J, 0, 1) # Ability scores

## True item-parameter matrix
# 1 gamma:Guessing # 2 alpha:Discrimination # 3 beta:Difficulty 
tpar <- rbind(rbeta(K, 3, 8), rlnorm(K, -0.3, 0.3), rnorm(K, 0.5, 0.7))

# Draw data
ip <- function(tpar, theta){
  tpar[1] + (1 - tpar[1]) / (1 + exp( - tpar[2] * (theta - tpar[3])))
} 
ps <- apply(tpar, 2, ip, theta = theta)
Y <- apply(ps, 2, rbinom, n = J, size = 1)

## JAGS
data.jags <- list("J" = J, "K" = K, "Y" = Y)

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
plot_pmeans <- function(mcmcres, K, J, tpar, title=''){
  grp = c(rep("alpha", K), rep("beta",K), rep("gamma", K), rep("theta", J))
  post.m <- apply(mcmcres, 2, mean)
  post.q <- apply(mcmcres, 2, quantile, c(0.025, 0.975))
  pdat = data.frame("est" = grp, 
                    "true" = c(tpar[2, ], tpar[3, ], tpar[1, ], theta), 
                    "pmean" = post.m, "lwr" = post.q[1, ], 
                    "upr" = post.q[2, ])
  p <- ggplot(pdat, aes(true, pmean))
  p <- p + geom_point()
  p <- p + geom_errorbar(aes(ymin = lwr, y = pmean, ymax = upr), width = 0, size = .5)
  p <- p + facet_wrap(~ est, ncol = 2, scales = "free")
  p <- p + labs(x = "True value", y = "Posterior Mean")
  p <- p + theme_bw()
  p <- p + ggtitle(title)
  p
}

jagspost1 <- jags.res[[2]][1][[1]]


plot_pmeans(jagspost1, K, J, tpar, "jags_3pl")
ggsave('plots/jags_3pl.pdf')

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
                   iter = 2000, warmup = 500, chains = 1, verbose = TRUE)
  )


stanpost <- do.call(cbind,stan.res@sim$samples[[1]][- (J + 3 * K + 1)])

# Plot posterior means
plot_pmeans(stanpost, K, J, tpar, "stan_3pl")
ggsave('plots/stan_3pl.pdf')

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
save(jags.res, file = "jags_res.RData")
#load("jags_res.RData")