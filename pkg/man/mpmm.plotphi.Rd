\name{mpmm.plotphi}
\alias{mpmm.plotphi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot model parameters.
}
\description{
This function is useful for exploring the parameter estimates from the model.
  }
\usage{
mpmm.plotphi(fit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{Output from running mpmm.cgs().}
}

\value{
  a plot faceted by each dimension (e.g. sender, receiver) showing MAP
  estimates for the probability of each individual being chosen for a given latent class.
}
\references{
Christopher DuBois and Padhraic Smyth. Modeling Relational Events via
Latent Classes.  Proceedings of the 16th ACM SIGKDD, 2010.
}
\author{
Christopher DuBois (<email: duboisc@ics.uci.edu>)
}

\examples{
## See demo(synthetic).
}
