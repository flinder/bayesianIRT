#################################################
# Comparison of JAGS and STAN for 2pl IRT models
#################################################
rm(list=ls())
library(rjags)
library(rstan)
library(ggplot2)

set.seed(2377343)

setwd("C:/Users/samsung/Dropbox/bd_irt/bayesianIRT")
#setwd("C:/Users/flinder/Dropbox/bd_irt/bayesianIRT")

### Generate some data according to:
# Y[j, k] ~ Bern(pi[j, k])
# pi[j, k] = 1/(1 + exp(-alpha[k] * (theta[j] - beta[k])))

J <- 300 # Number of subjects
K <- 20 # Number of items
N <- J * K # Number of observations
theta <- rnorm(J, 0, 1) # Ability scores

## True item-parameter matrix
# 1 alpha:Discrimination # 2 beta:Difficulty 
tpar <- rbind(rlnorm(K, -0.3, 0.3), rnorm(K, 0.5, 0.7))

# Draw data
ip <- function(tpar, theta){
  1 / (1 + exp( - tpar[1] * (theta - tpar[2])))
} 
ps <- apply(tpar, 2, ip, theta = theta)
Y <- apply(ps, 2, rbinom, n = J, size = 1)

## JAGS
data.jags <- list("J" = J, "K" = K, "Y" = Y)

# Starting values
inits <- list(list("alpha" = rep(1, K), # For chain 1
                   "beta" = rnorm(K, 0, 1.5), 
                   "theta" = rnorm(J)
                   )
              )
par.tosave <- c("alpha", "beta", "theta")

f <- function(){
  jmod <- jags.model(file = "models/jags_2pl_irt.txt", data = data.jags, 
                     inits = inits, n.chains = 1, n.adapt = 1000)
  mcmcres <- coda.samples(model = jmod, variable.names = par.tosave, 
                          n.iter = 5000, thin = 1)
  list(jmod, mcmcres)
}
t.jags <- system.time(jags.res <- f())

# Plot results
plot_pmeans <- function(mcmcres, K, J, tpar, title=''){
  grp = c(rep("alpha", K), rep("beta",K), rep("theta", J))
  post.m <- apply(mcmcres, 2, mean)
  post.q <- apply(mcmcres, 2, quantile, c(0.025, 0.975))
  pdat = data.frame("est" = grp, 
                    "true" = c(tpar[1, ], tpar[2, ], theta), 
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

plot_pmeans(jagspost1, K, J, tpar, '2pl jags') # From chain 1
ggsave('plots/jags_2pl.pdf')

## STAN
fileName <- "models/stan_2pl_irt.txt"
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
                   iter = 5000, warmup = 1000, chains = 1, verbose = TRUE)
  )

stanpost <- do.call(cbind,res.stan@sim$samples[[1]][- (J + 2 * K + 1)])

# Plot posterior means
plot_pmeans(stanpost, K, J, tpar, '2pl stan')
ggsave('plots/stan_2pl.pdf')
