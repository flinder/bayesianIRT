library(foreign)
library(MCMCpack)
library(ggplot2)

setwd("C:/Users/samsung/Dropbox/Fall2014/PLSC597D_MLE/project")
source("parallel_mcmc.R")

# ================================================================
# Simulated data
# ================================================================

# according to:
# Y[j, k] = 1 if z > 0, 0 otherwise
# z[j, k] = -beta[k] + alpha[k] * theta[j] + epsilon[j, k]

set.seed(79748)
J <- 5000 # Number of subjects
K <- 20  # Number of items
N <- J * K # Number of observations
theta <- rnorm(J, 0, 1) # Ability scores
alpha <- matrix(rep(rnorm(K, 0, 0.5), each = J) , ncol = K) # Difficulty
beta <- rlnorm(K, -1, 0.3) # Discrimination
epsilon <- matrix(rnorm(N), ncol = K)
Z <- -alpha + theta %*% t(beta) + epsilon # Latent score
Y <- ifelse(Z > 0, 1, 0)
rownames(Y) <- c(1:J)

# Function to create constraints list
get_const <- function(data, theta) {
  w_theta_s <- as.integer(rownames(data))
  theta <- theta[w_theta_s]
  names(theta) <- w_theta_s
  neg <- names(which.min(theta))
  pos <- names(which.max(theta))
  out <- list(min(theta), max(theta))
  names(out) <- c(paste0("V", neg), paste0("V", pos))
  return(out)
}

# ================================================================
# Fit Models
# ================================================================

n_burnin <- 1e3
n_samples <- 1e4

# Not paralelized benchmark model
# ----------------------------------------------------------------

t_full <- system.time(
  full_mod <- MCMCirt1d(Y, store.item = T, store.ability = F, burnin = n_burnin, 
                      mcmc = n_samples, verbose = F, 
                      theta.constraints = get_const(Y, theta))
)

# Parallelized with parametric combination of subposteriors
# ----------------------------------------------------------------

t_param <- system.time(
  prl_p_mod <- parallel_mcmc(data = Y, cores = 4, combine = "parametric", 
                             n_burnin = n_burnin, n_samples = n_samples, 
                             theta = theta)
  )


# Parallelized with non parametric combination of subposteriors (non parallel comb)
# ----------------------------------------------------------------

#t_nparam <- system.time(
#  prl_n_mod <- parallel_mcmc(data = Y, cores = 4, packages = "MCMCpack", 
#                             combine = "non parametric", fun = fit, 
#                             n_burnin = n_burnin, n_samples = n_samples)
#)


# Look at the results
# ----------------------------------------------------------------

# Point estimates 
truep <- rep(NA, 2 * K)
truep[seq(1, 2 * K, 2)] <- alpha[1, ]
truep[seq(2, 2 * K, 2)] <- beta

pdat <- data.frame(
  true = rep(truep, 2),
  p_est = c(apply(full_mod, 2, mean), prl_p_mod[[1]][[1]]),
  model = rep(c("Full", "Param"), each = K * 2),
  prm = rep(c("alpha", "beta"), K * 2)
  )

p <- ggplot(pdat, aes(true, p_est, color = prm))
p <- p + geom_point()
p <- p + facet_wrap( ~ model, scales = "fixed")
p


# Sub posteriors
sub_post <- prl_p_mod[[2]]
sub_means <- lapply(sub_post, function(x) apply(x, 2, mean))
pldat <- data.frame(est = unlist(sub_means), name = names(unlist(sub_means)),
                    true = rep(truep, 4),
                    grp = rep(c(1, 2, 3, 4), each = K),
                    type = rep(c("alpha", "beta"), K * 4))

p <- ggplot(pldat, aes(true, est, color = type))
p <- p + geom_point()
p <- p + facet_wrap( ~ grp, scales = "fixed")
p


# ================================================================
# Roll call data
# ================================================================
rc <- read.dta("hou112kh.dta")


# Recode: 1,2,3 yea 4,5,6 nay, 0,7,8,9: NA
rec <- function(x) {
  x[x %in% c(0, 7:9)] <- NA  
  x[x %in% c(1:3)] <- 1
  x[x %in% c(4:6)] <- 0
  return(x)
}
Y <- sapply(rc[, c(10:ncol(rc))], rec)
