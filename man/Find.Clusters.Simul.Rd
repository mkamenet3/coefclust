% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simultaneous.R
\name{Find.Clusters.Simul}
\alias{Find.Clusters.Simul}
\title{Find.Clusters.Simul}
\usage{
Find.Clusters.Simul(y, X, long, lat, MR, M, overlap, alpha)
}
\arguments{
\item{y}{The input data(as a vector)}

\item{X}{The input data(as a matrix)}

\item{long}{longitude}

\item{lat}{latitude}

\item{MR}{Maximum radius}

\item{M}{number of simulations}

\item{overlap}{Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters}

\item{alpha}{significance level}
}
\value{
list of cluster, coefficient, and indicator of cluster membership.
}
\description{
Find multiple clusters sequentially via simulataneous detection.
Find and test the cluster in the simple linear regression for given potential centroids via the simultaneous detection.
}
\examples{
data("SE_FakeData")
y <- SE_FakeData$y
X <- cbind(rep(1,length(y)), SE_FakeData$x)
long <- SE_FakeData$long; lat <- SE_FakeData$lat
MR <- 300; M <- 2
Clusters_sim <- Find.Clusters.Simul(y, X, long, lat, MR, M, overlap=TRUE, alpha=0.05)
}
