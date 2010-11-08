# In this demo we illustrate the use of algorithm on a small example
# with data generated from the MPMM model with dyadic data (2 vertices per event).

library(mpmm)
library(ggplot2)

# Initialize model parameters for a small 2-dimensional example
set.seed(123)
N <- 100
dims <- c(N,N)
C <- 4
pi <- c(.4,.2,.1,.3)
pi <- pi/sum(pi)
phi <- list()
phi[[1]] <- rbind(c(rep(1,60),rep(0,40)),
                  c(rep(0,80),rep(1,10),rep(0,10)),
                  c(rep(0,50),rep(1,50)),
                  c(rep(0,10),rep(1,20),rep(0,70)))
phi[[2]] <- rbind(c(rep(1,60),rep(0,40)),
                  c(rep(0,10),rep(1,60),rep(0,30)),
                  c(rep(0,80),rep(1,20)),
                  c(rep(0,70),rep(1,20),rep(0,10)))
phi <- lapply(phi,function(x) rownormalize(x+.001))
params <- list(pi=pi,phi=phi,C=C)

# Simulate synthetic dataset
T <- 1000
df <- mpmm.generate(T,params,dims)

# True predictive distribution of each event under the model
mpmm.plot(list(params=params,dims=dims))

# Observed data
obs <- mult.dir(df$dataset,dims)[,-4] # get counts, not probabilities
plotmat(obs)

# Set the priors use MCMC for inference
priors <- list(pi=1,phi=list(.1,.1,.1))
C <- 4
fit <- mpmm.cgs(df$dataset,C,dims,priors,niter=100)

# Plot the predictive distribution for a given sample from the posterior
mpmm.plot(fit)

# Plot the parameter estimates of phi
mpmm.plotphi(fit$params)

