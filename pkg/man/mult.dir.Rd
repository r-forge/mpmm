\name{mult.dir}
\alias{mult.dir}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a multinomial with a Dirichlet prior to the observed data.
}
\description{
A simple baseline for counts of event data is a multinomial with a symmetric
Dirichlet prior.  This method returns estimated probabilities for all
observed events.  (All unobserved events should be assigned a probability equal
to the prior of the Dirichlet.) 
}
\usage{
mult.dir(dataset, dims, prior = 1/prod(dims))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
A Tx2 or Tx3 matrix of integers, each row corresponding to an event, each
integer identifying the participating vertices.
}
  \item{dims}{
Maximum value for each categorical random variable (e.g. for dyadic
data, this might be the number of nodes in the graph).
in graph).
}
  \item{prior}{
The hyperparameter for a symmetric Dirichlet prior.
}
}
\value{
  a data.frame with the unique events from the provided dataset and two
  additional columns:
  \item{count}{the number of times that event occurred}
  \item{prob}{the posterior probability of each of these events given the observed
  data}
}

\author{
Christopher DuBois (<email: duboisc@ics.uci.edu>)
}

\examples{
data(enron)
dataset <- as.matrix(enron[,2:3])
N <- max(dataset)
dims <- c(N,N)
head(mult.dir(dataset,dims))
}
