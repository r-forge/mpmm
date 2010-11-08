\name{plotmat}
\alias{plotmat}
\title{
Plot a matrix.
}
\description{
An alternative to image().
}
\usage{
plotmat(mat, labels = c("Sender", "Receiver"),limits = c(0,max(mat[, 3])))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
A dataframe with 3 columns: d1, d2, and value
}
  \item{labels}{
The labels to use for the x- and y-axis.
}
\item{limits}{2 element vector specifying the lower and upper limits of
  value that will be colored white and black respectively.
}
}
\author{
Christopher DuBois
}
