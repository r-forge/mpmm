\name{mpmm.plotassignments}
\alias{mpmm.plotassignments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the latent class assigned to each observed event.
}
\description{
This function is useful for exploring the latent classes assigned to the
observed events.
}
\usage{
mpmm.plotassignments(dataset, fit, sort = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
A Tx2 or Tx3 matrix of integers, each row corresponding to an event, each
integer identifying the participating vertices.
  }
  \item{fit}{Output from running mpmm.cgs.}
  \item{sort}{
logical indicating whether the rows and columns should be sorted by most
likely per class.
  }
}

\author{
Christopher DuBois (<email: duboisc@ics.uci.edu>)
}

\examples{
## See demo(enron)
}
