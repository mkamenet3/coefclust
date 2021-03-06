% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twostage_SpTm.R
\name{Find.Clusters.TS.ST1}
\alias{Find.Clusters.TS.ST1}
\title{Find.Clusters.TS.ST1}
\usage{
Find.Clusters.TS.ST1(yList, XList, long, lat, MR, M, overlap, alpha)
}
\arguments{
\item{yList}{The input data (as a list of vectors).}

\item{XList}{The input data (as a list of matrices).}

\item{long}{longitude.}

\item{lat}{latitude.}

\item{MR}{Maximum radius.}

\item{M}{Number of simulations.}

\item{overlap}{Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters.}

\item{alpha}{Significance level}
}
\value{
List of clusters, coefficients, and indicator of cluster membership.
}
\description{
Find multiple (overlapping) cylindrical spatio-temporal clusters sequentially in the slopes via two-stage detection (in the 1st stage).
}
