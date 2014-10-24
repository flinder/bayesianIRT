# Function to plot mcmc results against true parameter values
plot_pmeans <- function(mcmcres, K, J, tpar, theta, title = ""){
  post.m <- apply(mcmcres, 2, mean)
  post.q <- apply(mcmcres, 2, quantile, c(0.025, 0.975))

  pars <- gsub("\\[\\d+\\]","",colnames(mcmcres))
  upars <- unique(pars)
  mod <- length(upars) #two pl or 3pl +1
  
  # Get correlation in par name for facet title
  for(i in 1:mod) {
    if(i == mod) r <- round(cor(theta, post.m[grep(upars[i], names(post.m))]), 2)
      else r <- round(cor(tpar[i, ], post.m[grep(upars[i], names(post.m))]), 2)
    pars <- gsub(upars[i], paste(upars[i], "r =", as.character(r)), pars)
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
