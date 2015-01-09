# =========================================================================
# Simulation for parallel_mcmc
# Fridolin Linder

# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------

setwd("C:/Users/samsung/Dropbox/Fall2014/PLSC597D_MLE/project")
library(doParallel)
library(ggplot2)
source("parallel_mcmc.R")

# Function to create constraints list
get_const <- function(data, theta) {
  w_theta_s <- as.integer(rownames(data))
  theta <- theta[w_theta_s]
  names(theta) <- w_theta_s
  neg <- names(which.min(theta))
  pos <- names(which.max(theta))
  out <- list(min(theta), max(theta))
  #names(out) <- c(paste0("V", neg), paste0("V", pos))
  names(out) <- c(as.character(neg), as.character(pos))
  return(out)
}

# Function to simulate data
gen_dat <- function(theta, alpha, beta, J, K) {
  N <- J*K
  epsilon <- matrix(rnorm(N), ncol = K)
  Z <- -alpha + theta %*% t(beta) + epsilon # Latent score
  Y <- ifelse(Z > 0, 1, 0)
  rownames(Y) <- c(1:J)
  return(Y)
}


# Parameters
B <- 50 # Number of MC simulations
n_burnin <- 1e3
n_samples <- 5e3
cores <- 8
K <- 10
J <- 4000
set.seed(123)
theta <- rnorm(J, 0, 1) # Ability scores
alpha <- matrix(rep(rnorm(K, 0, 0.5), each = J) , ncol = K) # Difficulty
beta <- rlnorm(K, -1, 0.3) # Discrimination
p <- 0.1 # *100 percent of the data for subsample method

# -------------------------------------------------------------------------
# Run Simulations
# -------------------------------------------------------------------------

### Run full models

# Function
full_run <- function() {
  Y <- gen_dat(theta, alpha, beta, J, K)
  t_full <- system.time(
    full_mod <- MCMCirt1d(Y, store.item = T, store.ability = F, burnin = n_burnin, 
                          mcmc = n_samples, verbose = T, 
                          theta.constraints = get_const(Y, theta))
  )
  means <- apply(full_mod, 2, mean)
  vars <- var(full_mod)
  out <- list("post_mean" = means, "post_var" = vars, "t" = t_full[3])
  return(out)  
}

# Run
cl <- makeCluster(cores)
registerDoParallel(cl) 
full_mods = foreach(i = seq(1, B), .packages = "MCMCpack") %dopar% {
                      full_run()
                    }
stopCluster(cl)

### Run parallelized models

# Function
par_run <- function() {
  Y <- gen_dat(theta, alpha, beta, J, K)
  t_param <- system.time(
    prl_p_mod <- parallel_mcmc(data = Y, cores = cores, combine = "parametric", 
                               n_burnin = n_burnin, n_samples = n_samples, 
                               theta = theta)
  )
  sub_means <- lapply(prl_p_mod[[2]], function(x) apply(x, 2, mean))
  sub_vars <- lapply(prl_p_mod[[2]], var)
  
  return(list(full = prl_p_mod[[1]], sub = list(sub_means, sub_vars), t = t_param[3]))
}

# Run
p_mods <- foreach(i = seq(1, B), .packages = "MCMCpack") %do% {
  par_run()
}

### Run Barberas method (subsamples)

# Function
barb_run <- function(){
  Y <- gen_dat(theta, alpha, beta, J, K)
  # Take p% sub sample of the data
  Ys <- Y[sample(c(1:J), p*J), ]
  t_barb <- system.time(
    barb_mod <- MCMCirt1d(Ys, store.item = T, store.ability = F, burnin = n_burnin, 
                          mcmc = n_samples, verbose = T, 
                          theta.constraints = get_const(Ys, theta))
  )
  means <- apply(barb_mod, 2, mean)
  vars <- var(barb_mod)
  out <- list("post_mean" = means, "post_var" = vars, "t" = t_barb[3])
  return(out)  
}

# Run
cl <- makeCluster(cores)
registerDoParallel(cl) 
barb_mods = foreach(i = seq(1, B), .packages = "MCMCpack") %dopar% {
  barb_run()
}
stopCluster(cl)

#save.image("par_simulation.RData")
#load("par_simulation.RData")

# ============================================================
# Visualize Results
# ============================================================

# -----------------------------------------------------------------------------
# Estimate plot

# Vector of true item parameters (in order as in the models)
truep <- rep(NA, 2 * K)
truep[seq(1, 2 * K, 2)] <- alpha[1, ]
truep[seq(2, 2 * K, 2)] <- beta

# Extract full estimates
full_means <- sapply(full_mods, function(x) x[[1]])

# Extract parallelized combined estimates
par_means <- sapply(p_mods, function(x) x$full$post_means)

# Extract subsample estimates
sub_means <- sapply(barb_mods, function(x) x$post_mean)

pdat <- data.frame(true = rep(truep, length.out = 3000),
                   est = c(c(full_means), c(par_means), c(sub_means)),
                   type = rep(c("Full", "Parallelized", "Sub sample"), each = 1000),
                   par = rep(c("alpha", "beta"), 1500)
                   )

q <- ggplot(pdat, aes(x = true, y = est, color = par))
q <- q + geom_point(alpha = 0.3)
q <- q + facet_wrap(~ type, scales = "fixed")
q <- q + theme_bw()
q <- q + xlab("True Parameter") + ylab("Estimate")
ggsave(plot = q, filename = "figures/estimates.png", height = 5, width = 15)

# Revert beta signs one model in full that is not identified
# so mse is not overly inflated
full_means[grep("beta", rownames(full_means)), 45]  <- -1*full_means[grep("beta", rownames(full_means)), 45]

# -----------------------------------------------------------------------------
# Plot for average computation time 

# Extract computation times
t_full <- sapply(full_mods, function(x) x[[3]])
t_par <- sapply(p_mods, function(x) x[[3]])
t_sub <- sapply(barb_mods, function(x) x[[3]])

tdat <- data.frame(t = c(t_full, t_par, t_sub),
                   type = rep(c("Full", "Parallelized", "Sub sample"), each = B))

q <- ggplot(tdat, aes(x = type, y = t))
q <- q + geom_boxplot()
q <- q + theme_bw()
q <- q + xlab("Method") + ylab("Estimation time in s")
ggsave(plot = q, filename = "figures/ctime.png", height = 5, width = 5)

# -----------------------------------------------------------------------------
# Plot for mse

# Calculate mean squared errors
mse_full <- apply((full_means - truep)^2, 1, mean)
mse_par <- apply((par_means - truep)^2, 1, mean)
mse_sub <- apply((sub_means - truep)^2, 1, mean)

mdat <- data.frame(mse = c(mse_full, mse_par, mse_sub),
                   par = rep(paste0(c("alpha", "beta"), rep(seq(1:10), each = 2)), 3),
                   type = rep(c("Full", "Parallelized", "Sub sample"), each = 20))
q <- ggplot(mdat, aes(fill = type, y = mse, x = par))
q <- q + geom_bar(position=position_dodge(), stat="identity")
q <- q + theme_bw()
q <- q + xlab("Parameter") + ylab("Mean Squared Error")
q <- q + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(plot = q, filename = "figures/mse.png", height = 5, width = 10)
