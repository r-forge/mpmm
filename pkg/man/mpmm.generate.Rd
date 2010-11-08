\name{mpmm.generate}
\alias{mpmm.generate}
\title{
Generate data from the MPMM.
}
\description{
[describe generative model]
}
\usage{
mpmm.generate(T, params, dims)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{Number of events to simulate.}
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
  \item{dataset}{ a matrix with T rows and 3 columns, where the d'th column
  contains integers in [1,...,N_d]}
  \item{assignments}{a vector of the class assigned to each event}
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
## See demo(synthetic) and demo(enron).
}
