\name{enron}
\alias{enron}
\docType{data}
\title{
Enron email data.
}
\description{
  The timestamp, sender, recipient, and message id of the emails.
[NOT COMPLETED]
}
\usage{data(enron)}
\format{
  A data frame with 34182 observations on the following 4 variables.
  \describe{
    \item{\code{time}}{a numeric vector}
    \item{\code{sender}}{a numeric vector}
    \item{\code{receiver}}{a numeric vector}
    \item{\code{mid}}{a numeric vector}
  }
}
\details{
This data was extracted from the MySQL dump.
}
\source{
  http://www.isi.edu/~adibi/Enron/Enron.htm
%  http://www.cs.cmu.edu/~enron/
}
\examples{
data(enron)
}
\keyword{datasets}
