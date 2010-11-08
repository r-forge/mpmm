library(mpmm)
library(ggplot2)
data(enron)
dataset <- as.matrix(enron[,2:3])

# Initialize model parameters for a small 2-dimensional example
set.seed(123)
N <- max(dataset)
dims <- c(N,N)
C <- 10

# Observed data
obs <- mult.dir(dataset,dims)[,-4] # get counts, not probabilities
plotmat(obs)

# Set the priors use MCMC for inference
priors <- list(pi=1,phi=list(1/N,1/N))
fit <- mpmm.cgs(dataset,C,dims,priors,niter=100)

# Plot the class assignments for the observed events after sorting rows and columns
mpmm.plotassignments(dataset,fit)

