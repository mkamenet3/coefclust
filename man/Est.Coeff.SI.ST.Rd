% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simultaneous_SpTm.R
\name{Est.Coeff.SI.ST}
\alias{Est.Coeff.SI.ST}
\title{Est.Coeff.SI.ST}
\usage{
Est.Coeff.SI.ST(yList, XList, long, lat, Clusters)
}
\arguments{
\item{yList}{The input data (as a list of vectors).}

\item{XList}{The input data (as a list of matrices).}

\item{long}{longitude.}

\item{lat}{latitude.}

\item{Clusters}{Inherited output from \code{Find.Clusters.SI.ST()}.}
}
\value{
List of coefficients
}
\description{
Estimate coefficients via simultaneous detection.
}
