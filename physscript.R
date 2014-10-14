rm(list=ls())
library(R2WinBUGS)
library(coda)
library(rjags)
library(rstan)

setwd("C:/Users/samsung/Dropbox/bd_irt/irt")
#setwd("C:/Users/flinder/Dropbox/bd_irt/phys212")

df = as.matrix(read.table('data/physics212posttest.csv',sep=','))

M <- ncol(df) # Number of Items
N <- nrow(df) # Number of subjects

# Set the initial values to be used.
a1 <- rep(1.1,M)
b1 <- rnorm(M,0,1.5)
c1 <- rep(0.2,M)
t1  <- rnorm(N)
inits <- list(list(b = b1, a = a1,c = c1, theta= t1))


data <- list("N"=N,"M"=M,"r"=df)

# WinBUGS

system.time(
cbstsim <- bugs(data,inits,model.file="models/threepmmodel.bug"
              ,parameters.to.save=c("a", "b","c","theta"), n.chains=1,n.burnin=1e3
              ,n.iter=1e4,inits=inits,n.thin=1
              ,bugs.directory="c:/Program Files (x86)/WinBUGS14/",debug=T)
)

# JAGS
par.tosave <- c("a", "b","c","theta")


a <- function(){
  jmod <- jags.model(file='models/threepmmodel.bug',data=data,inits=inits,n.chains=1,n.adapt=1e3)
  mcmcres <- coda.samples(model=jmod,variable.names=par.tosave,n.iter=1e4,thin=10)
  list(jmod,mcmcres)
}
system.time(a())
#user  system elapsed 
#1790.46    0.09 1790.77 

# STAN

