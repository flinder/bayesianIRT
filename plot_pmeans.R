# Function to plot mcmc results against true parameter values
plot_pmeans <- function(mcmcres, K, J, tpar, theta, title = ""){
  post.m <- apply(mcmcres, 2, mean)
  post.q <- apply(mcmcres, 2, quantile, c(0.025, 0.975))
  
  pars <- gsub("\\[\\d+\\]", "", colnames(mcmcres))
  upars <- unique(pars)
  
  # Work with stan hierarchical parametrization
  if(any(upars == "delta")){
    as <- grep("^alpha", colnames(mcmcres))
    mcmcres1 <- mcmcres
    mcmcres1[, as] <- exp(mcmcres[, as])
    # Add mean to all thetas
    ts <- grep("theta", colnames(mcmcres1))
    mcmcres1[, ts]  <- mcmcres1[, ts] + mean(mcmcres1[, "delta"])
    post.m <- apply(mcmcres1, 2, mean)
    post.q <- apply(mcmcres1, 2, quantile, c(0.025, 0.975))
    post.m[ts] <- post.m[ts] + post.m["delta"]    
    post.m <- post.m[-grep("delta", names(post.m))]
    post.m <- post.m[-grep("sigma_", names(post.m))]
    post.m <- post.m[-grep("lp_", names(post.m))]
    post.q <- post.q[, -grep("delta", colnames(post.q))]
    post.q <- post.q[, -grep("sigma_", colnames(post.q))]
    post.q <- post.q[, -grep("lp_", colnames(post.q))]
    
    pars <- gsub("\\[\\d+\\]", "", names(post.m))
    upars <- unique(pars)
    # Bring everything in the same order
    post.m <- c(post.m[grep("alpha",names(post.m))],
                post.m[grep("beta",names(post.m))],
                post.m[grep("theta",names(post.m))])
    post.q <- cbind(post.q[, grep("alpha",colnames(post.q))],
                    post.q[, grep("beta",colnames(post.q))],
                    post.q[, grep("theta",colnames(post.q))])
  }
  
  # Get correlation in par name for facet title
  for(i in upars) {
    if(i == "theta") r <- round(cor(theta, post.m[grep(i, names(post.m))]), 2)
      else r <- round(cor(tpar[i, ], post.m[grep(i, names(post.m))]), 2)
    pars <- gsub(i, paste(i, "r =", as.character(r)), pars)
  }
  
  pdat <- data.frame("est" = pars, 
                    "true" = c(c(t(tpar)), theta), 
                    "pmean" = post.m, "lwr" = post.q[1, ], 
                    "upr" = post.q[2, ])
  
  p <- ggplot(pdat, aes(true, pmean))
  p <- p + geom_point()
  p <- p + geom_errorbar(aes(ymin = lwr, y = pmean, ymax = upr), 
                         width = 0, size = .5, alpha = 0.4)
  p <- p + facet_wrap(~ est, ncol = 2, scales = "free")
  p <- p + labs(x = "True value", y = "Posterior Mean")
  p <- p + theme_bw()
  p <- p + ggtitle(title)
  p
}
