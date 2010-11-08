\name{mpmm.predict}
\alias{mpmm.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict probability of events given model parameters.
}
\description{
This function is used to compute the probability of given events under
the MPMM with the provided parameters.  These parameters could be
estimated from fitting the model with mpmm.cgs, for example.
  }
\usage{
mpmm.predict(dataset, params, dims)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
A Tx2 or Tx3 matrix of integers, each row corresponding to an event, each
integer identifying the participating vertices.
}
  \item{params}{
    A list with the following variables:
    \describe{
    \item{C}{number of latent classes}
    \item{pi}{vector of length C, denoting the probability of each latent class}
    \item{phi}{ a list with 3 elements (one for each dimension).  Each
    element is a C x N_d matrix denoting the probability of choosing a
    particular actor i given the event has latent class c (where N_d is
    the number of possible actors in dimension d).}
    }
  }
  \item{dims}{a vector containing the number of possible actors for each dimension}
}

\value{
\item{probs}{ a vector (with an element foreach row in the provided
dataset) containing the probability of that event under the model with
the provided parameters.
}
}
\references{
Christopher DuBois and Padhraic Smyth. Modeling Relational Events via
Latent Classes.  Proceedings of the 16th ACM SIGKDD, 2010.
}
\author{
Christopher DuBois (<email: duboisc@ics.uci.edu>)
}
\seealso{
mpmm.cgs
}
\examples{
data(enron)
dataset <- as.matrix(enron[,2:3])
N <- max(dataset)
dims <- c(N,N)
C <- 5
priors <- list(pi=1,phi=list(1/N,1/N))
fit <- mpmm.cgs(dataset,C,dims,priors,niter=50)

# Compute probability of edge (1,2) occurring under model
newdata <- matrix(c(1,2),1,2)
mpmm.predict(newdata,fit$params,dims)

# Compute loglikelihood of observed data
llk <- sum(log(mpmm.predict(dataset,fit$params,dims)))
}
