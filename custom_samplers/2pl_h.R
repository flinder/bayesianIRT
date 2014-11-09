rm(list = ls())
library(MCMCpack)
library(MASS)
library(mvtnorm)
set.seed(23429437)
### Generate some data according to:
# Y[j, k] ~ Bern(pi[j, k])
# pi[j, k] = 1/(1 + exp(-alpha[k] * (theta[j] - beta[k])))

J <- 500 # Number of subjects
K <- 10 # Number of items
N <- J * K # Number of observations
theta <- rnorm(J, 0, 1) # Ability scores
theta.t <- theta
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

### Functions

## Likelihood
# log-Likelihood for one row (one subject)
loglikp <- function(y, theta, zeta) {
  p <- 1/(1 + exp(-exp(zeta[1, ]) * (theta - zeta[2, ])))
  p[p == 1] = 0.999
  p[p == 0] = 0.0001
  sum(y * log(p) + (1 - y) * log(1 - p))
}

# log-Likelihood for one column (item)
logliki <- function(y, theta, zetat) {
  p <- 1/(1 + exp(-exp(zetat[1]) * (theta - zetat[2])))
  p[p == 1] = 0.999
  p[p == 0] = 0.0001
  sum(y * log(p) + (1 - y) * log(1 - p))
}


## Priors
# Theta
logp.t <- function(theta, mu_t, sig_t) {
  dnorm(theta, mu_t, sqrt(sig_t), log = T)
} 

# Zeta
logp.z <- function(zeta, mu_z, sig_z) {
  dmvnorm(zeta, mu_z, sig_z, log = T)
}

## Samplers


# Hyperparameters
samp.mu_t <- function(n, n0 = 1, mu_0, theta, sig_t) {
  mu <- n0 * mu_0 / (n + n0) + n * mean(theta) / (n + n0)
  sigma <- sig_t / (n + n0)
  rnorm(1, mu, sqrt(sigma))
}

samp.sig_t <- function(n, n0 = 1, g1, g2, theta, mu_t) {
  s <- 1 / (n - 1) * sum((theta - mu_t)^2)
  sig_n <- g2 + (n - 1) * s / 2 + n * n0 * (mean(theta) - mu_t)^2 / 2*(n + n0)
  rinvgamma(1, (g1 + n / 2), sig_n)
}

samp.mu_z <- function(K0, K, mu_0, sig_z, zeta) {
  zeta.b <- apply(zeta, 1, mean)
  mu <- (K0 * mu_0) / (K0 + K) + (K0 * zeta.b) / (K0 + K)
  sigma <- sig_z / (K0 + K)
  rmvnorm(1, mu, sigma)
}

samp.sig_z <- function(K, K0 = 1, sig_0) {
  zeta.b <- apply(zeta, 1, mean)
  dm <- matrix(zeta.b - mu_0)
  S <- sum(diag(t(zeta - zeta.b) %*% (zeta - zeta.b)))/K
  sig.s <- sig_0 + K * S + K * K0 / (K + K0) * dm %*% t(dm)
  riwish(K + 1, sig.s)
}

## Initialize parameters
# Primary Parameters
theta <- rep(0.1, J)
zeta <- rbind(rep(0.1, K),rep(0, K))
# Hyper-parameters
mu_t <- 0
mu_z <- c(1, 0)
sig_t <- 1
sig_z <- diag(2)

# Hyper-prior parameters
mu_0 <- c(1, 0)
sig_0 <- diag(2)
g1 <- 1/2
g2 <- 1/2

# Tuning parameters
phi_t  <- 1
phi_z  <- matrix(c(.1, 0, 0, 1), nc = 2)


# Number of mcmc iterations
B <- 300

# Result storage
res <- matrix(NA, nr = B, nc = (J + 2 * K))
colnames(res) <- c(paste0("theta[", c(1:J), "]"),
                   paste0("alpha[", c(1:K), "]"),
                   paste0("beta[", c(1:K), "]")
)

# MCMC Loop
for(t in 1:B) {
  cat("Iteration: ", t, "\r")
  # Step 1: M-H step for thetas
  theta.s <- rnorm(J, theta, phi_t)
  u <- runif(J)
  lnew <- lold <- rep(NA, J)
  for(j in 1:J) {
    lnew[j] <- loglikp(theta.s[j], y = Y[j, ], zeta = zeta)
    lold[j] <- loglikp(theta[j], y = Y[j, ], zeta = zeta)
  } 
  pnew <- logp.t(theta.s, mu_t, sig_t)
  pold <- logp.t(theta, mu_t, sig_t)
  ap <- exp(lnew * pnew - lold * pold)
  accept <- u < ap
  theta[accept] <- theta.s[accept]
  

  # Step 2: M-H step for zeta
  zeta.s <- apply(zeta, 2, mvrnorm, n = 1, Sigma = phi_z)
  u <- runif(K)
  lnew <- lold <- pnew <- pold <- rep(NA, K)
  for(k in 1:K) {
    lnew[k] <- logliki(theta = theta, y = Y[, k], zeta = zeta.s[, k])
    pnew[k] <- logp.z(zeta.s[, k], mu_z, sig_z)
    lold[k] <- logliki(theta = theta, y = Y[, k], zeta = zeta[, k])
    pold[k] <- logp.z(zeta[, k], mu_z, sig_z)
  }
  ap <- exp(lnew * pnew - lold * pold)
  accept <- u < ap
  zeta[, accept] <- zeta.s[, accept]
  
  # Step 3: Gibbs step for hyperparameters
  mu_t <- samp.mu_t(J, 1, mu_0, theta, sig_t)
  sig_t <- samp.sig_t(J, 1, g1, g2, theta, mu_t)
  mu_z <- samp.mu_z(K0 = 1, K = K, mu_0, sig_z, zeta)
  sig_z <- samp.sig_z(K0 = 1, K = K, sig_0)
  
  # Store results
  res[t, 1:J] <- theta
  res[t, (J+1):ncol(res)] <- c(exp(zeta[1, ]), zeta[2, ])
}


plot(c(1:B), res[, 4], type = "l")

theta.t[4]

