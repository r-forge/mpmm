\name{mpmm.cgs}
\alias{mpmm.cgs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Collapsed Gibbs sampler for the marginal product mixture model (MPMM)
}
\description{
Performs collapsed Gibbs sampling for the marginal-product mixture model
(DuBois 2010).  This function takes a dataset of observed relational
events, perform inference, and return point estimates of the
latent parameters using the state at the last iteration of Gibbs
sampling.

The model is appropriate for analyzing event data in 1, 2, or
3-mode graphs.  The main assumption is that each event belongs to a
latent class.  We use MCMC to explore the posterior distribution of
these assignments given observed data.  The resulting class assignments
can be useful for exploratory data analysis and predicting missing or future data.

At the moment, the number of latent classes must be chosen a priori,
though a DP version is in the works.

With the current implementation, each dimension is treated separately,
e.g. for a given latent class, there will be a distribution over likely
senders and a separate distribution over likely recipients, even if this
involves the same set of actors.
}
\usage{
mpmm.cgs(dataset, C, dims, priors, niter = 10, z.init = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
A Tx2 or Tx3 matrix of integers, each row corresponding to an event, each
integer identifying the participating vertices.
}
  \item{C}{
Number of desired latent classes.
}
  \item{dims}{
Maximum value for each categorical random variable.
}
  \item{priors}{
    \describe{
      \item{pi}{the value for the symmetric Dirichlet prior on pi
	(called alpha in the paper)}
      \item{phi}{a list where each element is the value for the
	symmetric Dirichlet prior for the corresponding phi vector}
    }
}
  \item{niter}{
Number of Gibbs iterations to perform.
}
  \item{z.init}{
vector of initial latent class assignments for all observations.
}
}

\value{
  \item{assignments}{Latent class assignments for the T observations
    from the last iteration of the Gibbs sampler.}
  \item{params}{
      \describe{
      \item{pi}{MAP estimate of the pi vector}
      \item{phi}{MAP estimate for the phi vectors}
    }
  }
  \item{counts}{A list where each element is a matrix counting the
    number of times values of the corresponding dimension were observed
    and assigned to a particular latent class. (These matrices are used
    for estimates of phi.)}
  \item{dims}{Maximum value for each categorical random variable.}
}
\references{
Christopher DuBois and Padhraic Smyth. Modeling Relational Events via
Latent Classes.  Proceedings of the 16th ACM SIGKDD, 2010.
}
\author{
Christopher DuBois (<email: duboisc@ics.uci.edu>)
}

\examples{
## See demo(synthetic) and demo(enron).
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
