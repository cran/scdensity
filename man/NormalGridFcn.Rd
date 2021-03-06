% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{NormalGridFcn}
\alias{NormalGridFcn}
\title{Compute the values of n normal PDFs (or their derivatives at m grid points).}
\usage{
NormalGridFcn(g, r, mu, sd)
}
\arguments{
\item{g}{Locations at which to evaluate the normal densities.}

\item{r}{Derivative degree of the densities to evaluate (0, 1, 2, or 3).}

\item{mu}{Means/locations for the normal densities.}

\item{sd}{Standard deviations of the normal densities.}
}
\value{
An m-by-n matrix as described above.
}
\description{
Grid points are specified in g. Returns an m-by-n matrix.  The (i,j)th element of
the matrix is the rth derivative of N(mu(j),sd(j)) at g(i).
}
\keyword{internal}
