###################################################################################################
# Last Update: 2016/04/29
#
# Find the Most Likely Cluster (MLC) with Longitude/Latitude coordinates
# Junho Lee, 2014/10/19 ~
#
### List.C.Data(DMatrix,MR):
###           construct a list of all possible clusters where r <= MR for each center
###
### cum.sum.scalar(yc,xc): compute cumulative sum of scalars
###
### MLC.CDL.ID.SL(y,X,cdataL,ID,overlap)
###           find the Most Likely Cluster in the Simple Linear Regression 
###           with a pre-defined cdataL which is from List.C.Data(DMatrix,MR)
###           for given potential centroids (ID).
###
### MLC.CDL2.ID.SL(y, x, cdataL, sum_x1, sum_x2, cs.xxL, ID, overlap)
###           find the Most Likely Cluster in the Simple Linear Regression 
###           a pre-defined cdataL which is from List.C.Data(DMatrix,MR),
###           a pre-computed sum_x1 which is sum(x_i)
###           a pre-computed sum_x2 which is sum(x_i^2)
###           a pre-computed cs.xxL which is from cum.sum.scalar1(xc)
###           for given potential centroids (ID).
###
### Fstat.CDL.ID.SL1(y,X,cdataL,ID,overlap)
###           find the Most Likely Cluster in the Slope in the Simple Linear Regression
###           with a pre-defined cdataL which is from List.C.Data(DMatrix,MR)
###           for given potential centroids (ID).
###
### Fstat.CDL.ID.SL2(y,X,cdataL,ID,overlap)
###           find the Most Likely Cluster in the INTercept only; the same slope but different intercept
###           in the Simple Linear Regression
###           with a pre-defined cdataL which is from List.C.Data(DMatrix,MR)
###           for given potential centroids (ID).
#
###################################################################################################

# if (!require(geosphere)) {
#   install.packages("geosphere")
#   library(geosphere)
# }
#library(Rcpp)
#sourceCpp("cDist2.cpp")

###################################################################################################
### Added 2014/11/06
### List.C.Data(DMatrix,MR):
###           construct a list of all possible clusters where r <= MR for each center
#### DMatrix: distance matrix
#### MR: The Max Radius

#'@title List.C.Data 
#'@description  Constructs a list of all possible clusters where r <= MR for each center.
#'@param DMatrix Distance matrix.
#'@param MR Maximum radius.
#'@return List of all possible clusters.
#'@export
List.C.Data <- function(DMatrix,MR) {
  cdata <- list()
  N <- dim(DMatrix)[1]
  for (k in 1:N) {
    dat.k <- data.frame(center=k, id=c(1:N), r=DMatrix[,k])
    dat.k1 <- subset(dat.k, r <= MR)
    dat.k1 <- dat.k1[order(dat.k1$r, dat.k1$id),]
    cdata[[k]] <- dat.k1
  }
  return(cdata)
}


###################################################################################################
### Added 2015/07/11
### Updated 2015/07/12
### cum.sum.scalar(yc,xc): compute cumulative sum of
###                      
#### yc: y vector
#### xc: x vector

#'@title cum.sum.scalar 
#'@description Compute the cumulative sum of scalars.
#'@param yc y vector.
#'@param xc x vector.
#'@return List.
cum.sum.scalar <- function(yc, xc) {
  n <- length(yc)
  y1 <- rep(NA,n); y2 <- rep(NA,n)
  x1 <- rep(NA,n); x2 <- rep(NA,n)
  xy <- rep(NA,n)
  
  y1[1] <- yc[1];      y2[1] <- yc[1]^2
  x1[1] <- xc[1];      x2[1] <- xc[1]^2
  xy[1] <- xc[1]*yc[1]
  
  for (i in 2:n) {
    y1[i] <- y1[i-1] + yc[i];     y2[i] <- y2[i-1] + yc[i]^2
    x1[i] <- x1[i-1] + xc[i];     x2[i] <- x2[i-1] + xc[i]^2
    xy[i] <- xy[i-1] + xc[i]*yc[i]
  }
  return(list(y1 = y1, y2 = y2, x1 = x1, x2 = x2, xy = xy, nc = 1:n))
}

#'@title cum.sum.scalar1 
#'@description Compute the cumulative sum by x vector.
#'@param xc x vector.
#'@return List.
cum.sum.scalar1 <- function(xc) {
  n <- length(xc)  
  x1 <- rep(NA,n); x2 <- rep(NA,n)
  
  x1[1] <- xc[1];      x2[1] <- xc[1]^2
  
  for (i in 2:n) {
    x1[i] <- x1[i-1] + xc[i];     x2[i] <- x2[i-1] + xc[i]^2
  }
  return(list(x1 = x1, x2 = x2, nc = 1:n))
}

#'@title cum.sum.scalar2 
#'@description Compute the cumulative sum.
#'@param yc y vector.
#'@param xc x vector.
#'@return List.
cum.sum.scalar2 <- function(yc, xc) {
  n <- length(yc)
  y1 <- rep(NA,n); y2 <- rep(NA,n)
  xy <- rep(NA,n)
  
  y1[1] <- yc[1];      y2[1] <- yc[1]^2
  xy[1] <- xc[1]*yc[1]
  
  for (i in 2:n) {
    y1[i] <- y1[i-1] + yc[i];     y2[i] <- y2[i-1] + yc[i]^2
    xy[i] <- xy[i-1] + xc[i]*yc[i]
  }
  return(list(y1 = y1, y2 = y2, xy = xy))
}

###################################################################################################
### Added 2015/07/11
### Updated 2015/11/11
### MLC.CDL.ID.SL(y,X,cdataL,ID,overlap): 
###      find the Most Likely Cluster in the Simple Linear Regressionfor given potential centroids
#### y: The input data(in Vector Format)
#### X: The input data(in Matrix Format)
#### cdataL: Pre-defined cdata list which is from List.C.Data(DMatrix,MR)
#### ID: Indices for potential centroids
#### overlap: Boolean which is TRUE for overlapping clusters / FALSE for non-overlapping clusters

#' @title MLC.CDL.ID.SL 
#' @description Find the most likely cluster in the simple linear regression for given centroids.
#' @param y The input data (in vector format).
#' @param X The input data (in matrix format).
#' @param cdataL Pre-defined cdata list which is from \code{List.C.Data(DMatrix,MR)}
#' @param ID Indices for potential centroids.
#' @param overlap Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters.
#' @return Most likely cluster (defined by center and radius) and \eqn{\hat{\sigma}^2}.
#' @export

MLC.CDL.ID.SL <- function(y, X, cdataL, ID, overlap) {
  mlc <- NULL
  N <- dim(X)[1]; p <- dim(X)[2]
  
  # compute sum(y_i), sum(y_i^2), sum(x_i), sum(x_i^2), sum(x_i*y_i)
  # Updated 2015/07/11
  x <- X[,p]
  sum_y1 <- sum(y);  sum_y2 <- sum(y^2)
  sum_x1 <- sum(x);  sum_x2 <- sum(x^2)
  sum_xy <- sum(x*y)
  
  wholeID <- 1:N                             # Update: 2015/06/22
  nonID <- wholeID[!(wholeID %in% ID)]       # Update: 2015/06/22
  
  for (k in ID) {
    dat.k <- cdataL[[k]]
    dat.k$dr <- round(dat.k$r - c(dat.k$r[-1],0),10) # Update: 2015/08/05
    # For non-overlapping assumpsion
    # Update: 2015/06/22
    if (!overlap && length(which(dat.k$id %in% nonID)) > 0) {
      # remove all cluster which overlaps the previous clusters  
      # also remove all cluster which is greater than overlapped one
      dat.k <- dat.k[-(min(which(dat.k$id %in% nonID)):length(dat.k$id)),]
    }
    if (!is.null(dim(dat.k)) && dim(dat.k)[1] > p) {  # Updated 2015/07/11
      yc <- y[dat.k$id]; xc <- x[dat.k$id]
      # compute cumulative sum
      cs <- cum.sum.scalar(yc, xc)
      # compute N*sigSq hat
      dat.k$Nsig2.hat<- sum_y2 - (sum_y1-cs$y1)^2/(N-cs$nc) - (cs$y1)^2/(cs$nc) -
        ((sum_xy-cs$xy) - (sum_x1-cs$x1)*(sum_y1-cs$y1)/(N-cs$nc))^2/((sum_x2-cs$x2) - (sum_x1-cs$x1)^2/(N-cs$nc)) -
        ((cs$xy) - (cs$x1)*(cs$y1)/(cs$nc))^2/((cs$x2) - (cs$x1)^2/(cs$nc))
      # Remove -Inf which is in case of all x's are the same
      dat.k2 <- as.matrix(dat.k[dat.k$dr!=0 & dat.k$Nsig2.hat > 0,])  # Update: 2015/11/11
      dim1 <- dim(dat.k2)[1]
      dim2 <- dim(dat.k2)[2]
      # Pick the smallest one. First p-1 is NA or NaN
      kpick <- order(dat.k2[,dim2])[1]
      kpick1 <- c(k, dat.k2[kpick,dim2-2], dat.k2[kpick,dim2]/N)
      mlc <- rbind(mlc,kpick1)
    }      
  }
  colnames(mlc) <- c("center","radius","sig2.hat")
  d2 <- dim(mlc)[2]
  # Pick the smallest one
  Mlc <- mlc[order(mlc[,d2])[1],]
  return(list(Mlc = Mlc, mlc = mlc))
}

###################################################################################################
### Added 2015/07/12
### Updated 2015/11/11
### MLC.CDL2.ID.SL(y, x, cdataL, sum_x1, sum_x2, cs.xxL, ID, overlap):
###      find the Most Likely Cluster in the Simple Linear Regressionfor given potential centroids
#### y: The input data(in Vector Format)
#### x: The input data(in Vector Format)
#### cdataL: Pre-defined cdata list which is from List.C.Data(DMatrix,MR)
#### sum_x1: Pre-computed sum(x_i)
#### sum_x2: Pre-computed sum(x_i^2)
#### cs.xxL: Pre-computed cs.xx list which is from cum.sum.scalar1(xc)
#### ID: Indices for potential centroids
#### overlap: Boolean which is TRUE for overlapping clusters / FALSE for non-overlapping clusters

#' @title MLC.CDL2.ID.SL
#' @description Find the most likely cluster in the simple linear regression for given potential centroids.
#' @param y The input data (in vector format).
#' @param x The input data (in vector format).
#' @param cdataL Pre-defined cdata list which is from \code{List.C.Data(DMatrix,MR)}.
#' @param sum_x1 Pre-computed \eqn{\sum x_i}.
#' @param sum_x2 Pre-computed \eqn{\sum x_i^2}.
#' @param cs.xxL Pre-computed \code{cs.xx} list which is from \code{cum.sum.scalar1(xc)}.
#' @param ID Indices for potential centroids.
#' @param overlap Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters.
#' @return Most likely cluster (defined by center and radius) and \eqn{\hat{\sigma}^2}.
#' @export

MLC.CDL2.ID.SL <- function(y, x, cdataL, sum_x1, sum_x2, cs.xxL, ID, overlap) {
  mlc <- NULL
  N <- length(y); p <- 2
  
  # compute sum(y_i), sum(y_i^2), sum(x_i*y_i)
  # Updated 2015/07/12
  sum_y1 <- sum(y);  sum_y2 <- sum(y^2)
  sum_xy <- sum(x*y)
  
  wholeID <- 1:N                             # Update: 2015/06/22
  nonID <- wholeID[!(wholeID %in% ID)]       # Update: 2015/06/22
  
  for (k in ID) {
    dat.k <- cdataL[[k]]
    dat.k$dr <- round(dat.k$r - c(dat.k$r[-1],0),10) # Update: 2015/08/05
    # For non-overlapping assumpsion
    # Update: 2015/06/22
    if (!overlap && length(which(dat.k$id %in% nonID)) > 0) {
      # remove all cluster which overlaps the previous clusters  
      # also remove all cluster which is greater than overlapped one
      dat.k <- dat.k[-(min(which(dat.k$id %in% nonID)):length(dat.k$id)),]
    }
    
    if (!is.null(dim(dat.k)) && dim(dat.k)[1] > p) {  # Update: 2015/06/22
      #dat.k$dr <- dat.k$r - c(dat.k$r[-1],0)      
      yc <- y[dat.k$id]; xc <- x[dat.k$id]
      
      # compute cumulative sum
      tmp<- cs.xxL[[k]]
      cs <- cum.sum.scalar2(yc, xc)
      x1 <- tmp$x1[1:length(yc)]
      x2 <- tmp$x2[1:length(yc)]
      nc <- tmp$nc[1:length(yc)]
      
      # compute N*sigSq hat
      dat.k$Nsig2.hat<- sum_y2 - (sum_y1-cs$y1)^2/(N-nc) - (cs$y1)^2/(nc) -
        ((sum_xy-cs$xy) - (sum_x1-x1)*(sum_y1-cs$y1)/(N-nc))^2/((sum_x2-x2) - (sum_x1-x1)^2/(N-nc)) -
        ((cs$xy) - (x1)*(cs$y1)/(nc))^2/((x2) - (x1)^2/(nc))
      # Remove -Inf which is in case of all x's are the same
      dat.k2 <- as.matrix(dat.k[dat.k$dr!=0 & dat.k$Nsig2.hat > 0,])  # Update: 2015/11/11
      
      dim1 <- dim(dat.k2)[1]
      dim2 <- dim(dat.k2)[2]
      
      # Pick the smallest one. First p-1 is NA or NaN
      kpick <- order(dat.k2[,dim2])[1]
      kpick1 <- c(k, dat.k2[kpick,dim2-2], dat.k2[kpick,dim2]/N)
      mlc <- rbind(mlc,kpick1)
    }      
  }
  colnames(mlc) <- c("center","radius","sig2.hat")
  d2 <- dim(mlc)[2]
  # Pick the smallest one  
  Mlc <- mlc[order(mlc[,d2])[1],]
  return(list(Mlc = Mlc, mlc = mlc))
}

###################################################################################################
### Added 2015/07/12
### Updated 2015/11/11
### Fstat.CDL.ID.SL1(y,X,cdataL,ID,overlap):
###      find the Most Likely Cluster in the Simple Linear Regressionfor given potential centroids
###      in the Slope: different slope and different intercept
#### y: The input data(in Vector Format)
#### X: The input data(in Matrix Format)
#### cdataL: Pre-defined cdata list which is from List.C.Data(DMatrix,MR)
#### ID: Indices for potential centroids
#### overlap: Boolean which is TRUE for overlapping clusters / FALSE for non-overlapping clusters



#' @title Fstat.CDL.ID.SL1 
#' @description Find the most likely cluster in the simple linear regression for given potential centroids in the slope: different slope and different intercept.
#' @param y The input data (in vector format).
#' @param X The input data (in matrix format).
#' @param cdataL Pre-defined cdata list which is from \code{List.C.Data(DMatrix,MR)}.
#' @param ID Indices for potential centroids.
#' @param overlap Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters.
#' @return Most likely cluster (defined by center and radius), \eqn{\hat{\sigma}^2}, and maximum F-statistic.
#' @export
Fstat.CDL.ID.SL1 <- function(y, X, cdataL, ID, overlap) {
  mlc <- NULL
  N <- dim(X)[1]; p <- dim(X)[2]
  
  # compute sum(y_i), sum(y_i^2), sum(x_i), sum(x_i^2), sum(x_i*y_i)
  x <- X[,p]
  sum_y1 <- sum(y);  sum_y2 <- sum(y^2)
  sum_x1 <- sum(x);  sum_x2 <- sum(x^2)
  sum_xy <- sum(x*y)
  
  wholeID <- 1:N
  nonID <- wholeID[!(wholeID %in% ID)]
  
  for (k in ID) {
    dat.k <- cdataL[[k]]
    dat.k$dr <- round(dat.k$r - c(dat.k$r[-1],0),10) # Update: 2015/08/05
    # For non-overlapping assumpsion
    if (!overlap && length(which(dat.k$id %in% nonID)) > 0) {
      # remove all cluster which overlaps the previous clusters  
      # also remove all cluster which is greater than overlapped one
      dat.k <- dat.k[-(min(which(dat.k$id %in% nonID)):length(dat.k$id)),]
    }
    
    if (!is.null(dim(dat.k)) && dim(dat.k)[1] > p) {
      yc <- y[dat.k$id]; xc <- x[dat.k$id]
      
      # compute cumulative sum
      cs <- cum.sum.scalar(yc, xc)
      
      # compute N*sigSq hat
      dat.k$Nsig2.hat<- sum_y2 - (sum_y1-cs$y1)^2/(N-cs$nc) - (cs$y1)^2/(cs$nc) -
        ((sum_xy-cs$xy) - (sum_x1-cs$x1)*(sum_y1-cs$y1)/(N-cs$nc))^2/((sum_x2-cs$x2) - (sum_x1-cs$x1)^2/(N-cs$nc)) -
        ((cs$xy) - (cs$x1)*(cs$y1)/(cs$nc))^2/((cs$x2) - (cs$x1)^2/(cs$nc))
      
      # compute theta1
      theta1 <- ((cs$xy) - (cs$x1)*(cs$y1)/(cs$nc))/((cs$x2) - (cs$x1)^2/(cs$nc)) -
        ((sum_xy-cs$xy) - (sum_x1-cs$x1)*(sum_y1-cs$y1)/(N-cs$nc))/((sum_x2-cs$x2) - (sum_x1-cs$x1)^2/(N-cs$nc))
      
      # compute Fstat.coeff
      dat.k$Fstat.1 <- (theta1^2)/(1/((cs$x2) - (cs$x1)^2/(cs$nc)) + 1/((sum_x2-cs$x2) - (sum_x1-cs$x1)^2/(N-cs$nc)))/
        ((dat.k$Nsig2.hat)/(N-2*p))
      
      # Remove -Inf which is in case of all x's are the same
      dat.k2 <- as.matrix(dat.k[dat.k$dr!=0 & dat.k$Nsig2.hat > 0,])  # Update: 2015/11/11
      
      dim1 <- dim(dat.k2)[1]
      dim2 <- dim(dat.k2)[2]
      
      # Pick the largest Fstat.1. First p-1 is NA or NaN
      kpick <- order(-dat.k2[,dim2])[1]
      mlc <- rbind(mlc, dat.k2[kpick,-c(2,4)])
    }
  }
  
  mlc[,3] <- mlc[,3]/N
  colnames(mlc) <- c("center","radius","sig2.hat","supFstat.1")
  d2 <- dim(mlc)[2]
  
  # Pick the largest Fstat.1.
  Mlc <- mlc[order(-mlc[,d2])[1],]
  
  return(list(Mlc = Mlc, mlc = mlc))
}

###################################################################################################
### Added 2015/07/12
### Updated 2015/11/11
### Fstat.CDL.ID.SL2(y,X,cdataL,ID,overlap):
###      find the Most Likely Cluster in the Simple Linear Regressionfor given potential centroids
###      in the INTercept only: the same slope but different intercept
#### y: The input data(in Vector Format)
#### X: The input data(in Matrix Format), N-by-2 matrix
#### cdataL: Pre-defined cdata list which is from List.C.Data(DMatrix,MR)
#### ID: Indices for potential centroids
#### overlap: Boolean which is TRUE for overlapping clusters / FALSE for non-overlapping clusters

#' @title Fstat.CDL.ID.SL2
#' @description Find the most likely cluster in the simple linear regression for given potential centroids in the intercept only: the same slope but different intercepts.
#' @param y The input data (in vector format).
#' @param X The input data (in matrix format), N-by-2 matrix.
#' @param cdataL Pre-defined cdata list which is from \code{List.C.Data(DMatrix,MR)}.
#' @param ID Indices for potential centroids.
#' @param overlap Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters.
#' @return Most likely cluster (defined by center and radius) and \eqn{\hat{\sigma}^2}.
#' @export
Fstat.CDL.ID.SL2 <- function(y, X, cdataL, ID, overlap) {
  mlc <- NULL
  N <- dim(X)[1]; p <- dim(X)[2]
  
  # compute sum(y_i), sum(y_i^2), sum(x_i), sum(x_i^2), sum(x_i*y_i)
  x <- X[,p]
  sum_y1 <- sum(y);  sum_y2 <- sum(y^2)
  sum_x1 <- sum(x);  sum_x2 <- sum(x^2)
  sum_xy <- sum(x*y)
  
  wholeID <- 1:N                             # Update: 2015/06/22
  nonID <- wholeID[!(wholeID %in% ID)]       # Update: 2015/06/22
  
  for (k in ID) {
    dat.k <- cdataL[[k]]
    dat.k$dr <- round(dat.k$r - c(dat.k$r[-1],0),10) # Update: 2015/08/05
    # For non-overlapping assumption
    # Update: 2015/06/22
    if (!overlap && length(which(dat.k$id %in% nonID)) > 0) {
      # remove all cluster which overlaps the previous clusters  
      # also remove all cluster which is greater than overlapped one
      dat.k <- dat.k[-(min(which(dat.k$id %in% nonID)):length(dat.k$id)),]
    }
    if (!is.null(dim(dat.k)) && dim(dat.k)[1] >= p) {  # the same slope / different intercept
      yc <- y[dat.k$id]; xc <- x[dat.k$id]
      
      # compute cumulative sum
      cs <- cum.sum.scalar(yc, xc)
      
      # compute N*sigSq hat
      dat.k$Nsig2.hat<- sum_y2 -
        (((cs$nc)*sum_x2-cs$x1^2)*sum_y1^2 + (cs$nc)*(N-(cs$nc))*sum_xy^2 + 
           (N*sum_x2-sum_x1^2)*cs$y1^2 - 2*(cs$nc)*(sum_x1-cs$x1)*sum_y1*sum_xy -
           2*((cs$nc)*sum_x2-sum_x1*cs$x1)*sum_y1*cs$y1 + 
           2*((cs$nc)*sum_x1-N*cs$x1)*sum_xy*cs$y1) / 
        ((cs$nc)*((N-(cs$nc))*(sum_x2-cs$x1^2/(cs$nc)) - (sum_x1-cs$x1)^2))
      
      # Remove -Inf which is in case of all x's are the same
      dat.k2 <- as.matrix(dat.k[dat.k$dr!=0 & dat.k$Nsig2.hat > 0,])  # Update: 2015/11/11
      
      dim1 <- dim(dat.k2)[1]
      dim2 <- dim(dat.k2)[2]
      
      # Pick the smallest one. First p-1 is NA(??)
      kpick <- order(dat.k2[,dim2])[1]
      
      kpick1 <- c(k, dat.k2[kpick,dim2-2], dat.k2[kpick,dim2]/N)
      mlc <- rbind(mlc,kpick1)
    }      
  }
  
  colnames(mlc) <- c("center","radius","sig2.hat")     # Update: 2015/07/06
  d2 <- dim(mlc)[2]
  
  # Pick the smallest one
  Mlc <- mlc[order(mlc[,d2])[1],]
  
  return(list(Mlc = Mlc, mlc = mlc))
}


###################################################################################################
# For the ANOVA model with 2 obs (binary covariate) in a site
###################################################################################################
### Added 2015/10/28
### Updated 2015/10/29
### MLC.CDL.ID.1AOVbi(y,X,cdataL,ID,overlap): 
###      find the Most Likely Cluster in the one-way ANOVA with a binary covariate
###      for given potential centroids
#### y: The input data(in Vector Format, 2N-by-1)
#### X: The input data(in Matrix Format, 2N-by-2)
#### cdataL: Pre-defined cdata list which is from List.C.Data(DMatrix,MR)
#### ID: Indices for potential centroids
#### overlap: Boolean which is TRUE for overlapping clusters / FALSE for non-overlapping clusters

#'@title MLC.CDL.ID.1AOVbi
#'@description Find the most likely cluster in the one-way ANOVA with a binary covariate for given potential centroids.
#'@param y The input data (in vector format, 2N-by-1).
#'@param X The input data (in matrix format, 2N-by-2).
#'@param cdataL Pre-defined cdata list which is from \code{List.C.Data(DMatrix,MR)}.
#'@param ID Indices for potential centroids.
#'@param overlap  Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters
#'@return Most likely cluster (defined by center and radius) and \eqn{\hat{\sigma}^2}.
#'@export
#'
MLC.CDL.ID.1AOVbi <- function(y, X, cdataL, ID, overlap) {
  mlc <- NULL
  N <- length(cdataL); p <- dim(X)[2]
  
  # compute sum(y_i), sum(y_i^2), sum(x_i), sum(x_i^2), sum(x_i*y_i)
  # Updated 2015/07/11
  x <- X[,p]
  sum_y1 <- sum(y);  sum_y2 <- sum(y^2)
  sum_xy <- sum(x*y)
  
  wholeID <- 1:N                             # Update: 2015/06/22
  nonID <- wholeID[!(wholeID %in% ID)]       # Update: 2015/06/22
  
  for (k in ID) {
    dat.k <- cdataL[[k]]
    dat.k$dr <- round(dat.k$r - c(dat.k$r[-1],0),10) # Update: 2015/08/05
    # For non-overlapping assumpsion
    # Update: 2015/06/22
    if (!overlap && length(which(dat.k$id %in% nonID)) > 0) {
      # remove all cluster which overlaps the previous clusters  
      # also remove all cluster which is greater than overlapped one
      dat.k <- dat.k[-(min(which(dat.k$id %in% nonID)):length(dat.k$id)),]
    }
    
    if (!is.null(dim(dat.k)) && dim(dat.k)[1] > p) {  # Updated 2015/07/11
      id.k <- c(rbind(2*dat.k$id-1,2*dat.k$id))  # Updated 2015/10/28
      yc <- y[id.k]; xc <- x[id.k]
      
      # compute cumulative sum
      cs <- cum.sum.scalar2(yc, xc)
      cs.id.k <- (1:(length(yc)/2))*2
      cs.y1 <- cs$y1[cs.id.k]
      cs.xy <- cs$xy[cs.id.k]
      cs.nc <- 1:length(cs.id.k)
      
      # compute SSE = 2N*sigSq hat
      dat.k$SSE <- sum_y2 - (sum_y1-cs.y1)^2/(2*(N-cs.nc)) - (cs.y1)^2/(2*cs.nc) -
        ((sum_xy-cs.xy) - (sum_y1-cs.y1)/2)^2/((N-cs.nc)/2) - (cs.xy - cs.y1/2)^2/(cs.nc/2)
      
      # Remove -Inf which is in case of all x's are the same
      dat.k2 <- as.matrix(dat.k[dat.k$dr!=0 & dat.k$SSE > 0,])  # Update: 2015/11/11
      
      dim1 <- dim(dat.k2)[1]
      dim2 <- dim(dat.k2)[2]
      
      # Pick the smallest one.
      kpick <- order(dat.k2[,dim2])[1]
      kpick1 <- c(k, dat.k2[kpick,dim2-2], dat.k2[kpick,dim2]/(2*N))
      mlc <- rbind(mlc,kpick1)
    }      
  }
  
  colnames(mlc) <- c("center","radius","sig2.hat")
  d2 <- dim(mlc)[2]
  
  # Pick the smallest one
  Mlc <- mlc[order(mlc[,d2])[1],]
  
  return(list(Mlc = Mlc, mlc = mlc))
}


###################################################################################################
### Added 2016/02/05
### Updated 2016/02/29
### MLC.CDL.ID.SL.ST(yList, XList, cdataL, ID, overlap): 
###      find the Spatio-Temporal Cluster estimate in the Simple Linear Regression for given potential centroids
#### yList: The input data(in the List of Vectors Format)
#### XList: The input data(in the List of Matrices Format)
#### cdataL: Pre-defined cdata list which is from List.C.Data(DMatrix,MR)
#### ID: Indices for potential centroids
#### overlap: Boolean which is TRUE for overlapping clusters / FALSE for non-overlapping clusters

#'@title MLC.CDL.ID.SL.ST
#'@description Find the spatio-temporal cluster estimate in the simple linear regression for given potential centroids.
#'@param yList The input data (as a list of vectors).
#'@param XList The input data (as a list of matrices).
#'@param cdataL Pre-defined cdata list which is from \code{List.C.Data(DMatrix,MR)}.
#'@param ID Indices for potential centroids.
#'@param overlap  Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters
#'@return Most likely cluster (defined by center and radius) and SSE (error sum of squares).
#'@export
MLC.CDL.ID.SL.ST <- function(yList, XList, cdataL, ID, overlap) {
  mlc <- NULL
  T <- length(XList)
  N <- dim(XList[[1]])[1]; p <- dim(XList[[1]])[2]
  
  wholeID <- 1:N                             # Update: 2015/06/22
  nonID <- wholeID[!(wholeID %in% ID)]       # Update: 2015/06/22
  
  # compute sum(y_it), sum(y_it^2), sum(x_it), sum(x_it^2), sum(x_it*y_it) for each t=1,...T fot all i=1,...,N
  # Updated 2016/02/29
  sum_y1L <- list();  sum_y2L <- list();  sum_x1L <- list();  sum_x2L <- list();  sum_xyL <- list();
  for (t in 1:T) {
    y <- yList[[t]]
    x <- XList[[t]][,p]; x[is.na(y)] <- NA
    sum_y1L[[t]] <- sum(y, na.rm=TRUE);  sum_y2L[[t]] <- sum(y^2, na.rm=TRUE)
    sum_x1L[[t]] <- sum(x, na.rm=TRUE);  sum_x2L[[t]] <- sum(x^2, na.rm=TRUE)
    sum_xyL[[t]] <- sum(x*y, na.rm=TRUE)
  }
  
  for (k in ID) {
    dat.k <- cdataL[[k]]
    dat.k$dr <- round(dat.k$r - c(dat.k$r[-1],0),10) # Update: 2015/08/05
    # For non-overlapping assumption
    # Update: 2015/06/22
    if (!overlap && length(which(dat.k$id %in% nonID)) > 0) {
      # remove all cluster which overlaps the previous clusters  
      # also remove all cluster which is greater than overlapped one
      dat.k <- dat.k[-(min(which(dat.k$id %in% nonID)):length(dat.k$id)),]
    }
    
    if (!is.null(dim(dat.k)) && dim(dat.k)[1] > p) {  # Updated 2015/07/11
      SSE.kList <- list()
      for (t in 1:T) {
        sum_y1 <- sum_y1L[[t]];  sum_y2 <- sum_y2L[[t]]
        sum_x1 <- sum_x1L[[t]];  sum_x2 <- sum_x2L[[t]]
        sum_xy <- sum_xyL[[t]]
        yc <- yList[[t]][dat.k$id]; xc <- XList[[t]][dat.k$id,p]
        
        length_yc <- length(yc)      # Updated 2016/02/29
        nc <- 1:length_yc            # Updated 2016/02/29
        if (sum(is.na(yc)) > 0) {    # Updated 2016/02/29
          na_yc <- (1:length_yc)[is.na(yc)]
          yc[na_yc] <- xc[na_yc] <- 0
          for (j in na_yc) {
            nc[j:length_yc] <- nc[j:length_yc] - 1
          }
        }
        
        # compute cumulative sum
        cs <- cum.sum.scalar(yc, xc)
        cs$nc <- nc                       # Updated 2016/02/29
        num_na<- sum(is.na(yList[[t]]))   # Updated 2016/02/29
        SSE.kList[[t]] <- sum_y2 - (sum_y1-cs$y1)^2/(N-cs$nc-num_na) - (cs$y1)^2/(cs$nc) -
          ((sum_xy-cs$xy) - (sum_x1-cs$x1)*(sum_y1-cs$y1)/(N-cs$nc-num_na))^2/((sum_x2-cs$x2) - (sum_x1-cs$x1)^2/(N-cs$nc-num_na)) -
          ((cs$xy) - (cs$x1)*(cs$y1)/(cs$nc))^2/((cs$x2) - (cs$x1)^2/(cs$nc))   # Updated 2016/02/29
      }
      
      # compute SSE
      dat.k$SSE<- Reduce("+",SSE.kList)
      
      # Remove -Inf which is in case of all x's are the same
      dat.k2 <- as.matrix(dat.k[dat.k$dr!=0 & dat.k$SSE > 0,])
      
      dim1 <- dim(dat.k2)[1]
      dim2 <- dim(dat.k2)[2]
      
      # Pick the smallest one. First p-1 is NA or NaN
      kpick <- order(dat.k2[,dim2])[1]
      kpick1 <- c(k, dat.k2[kpick,dim2-2], dat.k2[kpick,dim2])
      mlc <- rbind(mlc,kpick1)
    }      
  }
  
  colnames(mlc) <- c("center","radius","SSE")
  d2 <- dim(mlc)[2]
  
  # Pick the smallest one
  Mlc <- mlc[order(mlc[,d2])[1],]
  
  return(list(Mlc = Mlc, mlc = mlc))
}


###################################################################################################
### Added 2016/02/05
### Updated 2016/03/01
### MLC.CDL2.ID.SL.ST(yList, xList, cdataL, sum_x1L, sum_x2L, cs.xxL.L, ID, overlap): 
###      find the Spatio-Temporal Cluster estimate in the Simple Linear Regression for given potential centroids
#### yList: The input data(in the List of Vectors Format)
#### xList: The input data(in the List of Vectors Format)
#### cdataL: Pre-defined cdata list which is from List.C.Data(DMatrix,MR)
#### sum_x1L: List of Pre-computed sum(x_i)
#### sum_x2L: List of Pre-computed sum(x_i^2)
#### cs.xxL.L: List of Pre-computed cs.xx list which is from cum.sum.scalar1(xc)
#### ID: Indices for potential centroids
#### overlap: Boolean which is TRUE for overlapping clusters / FALSE for non-overlapping clusters


#'@title MLC.CDL2.ID.SL.ST
#'@description Find the spatio-temporal cluster estimate in the simple linear regression for given potential centroids.
#'@param yList The input data (as a list of vectors).
#'@param xList The input data (as a list of vectors).
#'@param cdataL Pre-defined cdata list which is from \code{List.C.Data(DMatrix,MR)}.
#'@param sum_x1L List of pre-computed \eqn{\sum x_i}.
#'@param sum_x2L List of pre-computed \eqn{\sum x_i^2}.
#'@param cs.xxL.L List of pre-computed \code{cs.xx} list which is from \code{cum.sum.scalar1(xc)}.
#'@param ID Indices for potential centroids.
#'@param overlap  Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters
#'@return Most likely cluster (defined by center and radius) and SSE (error sum of squares).
#'@export
MLC.CDL2.ID.SL.ST <- function(yList, xList, cdataL, sum_x1L, sum_x2L, cs.xxL.L, ID, overlap) {
  mlc <- NULL
  T <- length(XList)
  N <- dim(XList[[1]])[1]; p <- dim(XList[[1]])[2]
  
  wholeID <- 1:N                             # Update: 2015/06/22
  nonID <- wholeID[!(wholeID %in% ID)]       # Update: 2015/06/22
  
  # compute sum(y_it), sum(y_it^2), sum(x_it*y_it) for each t=1,...T fot all i=1,...,N
  # Updated 2016/03/01
  sum_y1L <- list();  sum_y2L <- list();  sum_xyL <- list();
  for (t in 1:T) {
    x <- xList[[t]]
    y <- yList[[t]]
    sum_y1L[[t]] <- sum(y, na.rm=TRUE);  sum_y2L[[t]] <- sum(y^2, na.rm=TRUE)
    sum_xyL[[t]] <- sum(x*y, na.rm=TRUE)
  }
  
  for (k in ID) {
    dat.k <- cdataL[[k]]
    dat.k$dr <- round(dat.k$r - c(dat.k$r[-1],0),10) # Update: 2015/08/05
    # For non-overlapping assumpsion
    # Update: 2015/06/22
    if (!overlap && length(which(dat.k$id %in% nonID)) > 0) {
      # remove all cluster which overlaps the previous clusters  
      # also remove all cluster which is greater than overlapped one
      dat.k <- dat.k[-(min(which(dat.k$id %in% nonID)):length(dat.k$id)),]
    }
    
    if (!is.null(dim(dat.k)) && dim(dat.k)[1] > p) {  # Update: 2015/06/22
      SSE.kList <- list()
      for (t in 1:T) {
        sum_y1 <- sum_y1L[[t]];  sum_y2 <- sum_y2L[[t]]
        sum_x1 <- sum_x1L[[t]];  sum_x2 <- sum_x2L[[t]]
        sum_xy <- sum_xyL[[t]]
        cs.xxL <- cs.xxL.L[[t]]
        yc <- yList[[t]][dat.k$id]; xc <- xList[[t]][dat.k$id]
        
        if (sum(is.na(yc)) > 0) {    # Updated 2016/03/01
          yc[is.na(yc)] <- xc[is.na(yc)] <- 0
        }
        # compute cumulative sum
        tmp<- cs.xxL[[k]]
        cs <- cum.sum.scalar2(yc, xc)
        x1 <- tmp$x1[1:length(yc)]
        x2 <- tmp$x2[1:length(yc)]
        nc <- tmp$nc[1:length(yc)]
        num_na<- sum(is.na(yList[[t]]))   # Updated 2016/03/01
        SSE.kList[[t]] <- sum_y2 - (sum_y1-cs$y1)^2/(N-nc-num_na) - (cs$y1)^2/(nc) -
          ((sum_xy-cs$xy) - (sum_x1-x1)*(sum_y1-cs$y1)/(N-nc-num_na))^2/((sum_x2-x2) - (sum_x1-x1)^2/(N-nc-num_na)) -
          ((cs$xy) - (x1)*(cs$y1)/(nc))^2/((x2) - (x1)^2/(nc))
      }
      
      # compute SSE
      dat.k$SSE<- Reduce("+",SSE.kList)
      
      # Remove -Inf which is in case of all x's are the same
      dat.k2 <- as.matrix(dat.k[dat.k$dr!=0 & dat.k$SSE > 0,])
      
      dim1 <- dim(dat.k2)[1]
      dim2 <- dim(dat.k2)[2]
      
      # Pick the smallest one. First p-1 is NA or NaN
      kpick <- order(dat.k2[,dim2])[1]
      #kpick1 <- c(k, dat.k2[kpick,dim2-2], dat.k2[kpick,dim2]/N)
      kpick1 <- c(k, dat.k2[kpick,dim2-2], dat.k2[kpick,dim2])
      mlc <- rbind(mlc,kpick1)
    }      
  }
  
  colnames(mlc) <- c("center","radius","SSE")
  d2 <- dim(mlc)[2]
  
  # Pick the smallest one  
  Mlc <- mlc[order(mlc[,d2])[1],]
  
  return(list(Mlc = Mlc, mlc = mlc))
}


###################################################################################################
### Added 2016/02/07
### Updated 2016/03/06
### Fstat.CDL.ID.SL1.ST(yList, XList, cdataL, ID, overlap):
###      find the Spatio-Temporal Cluster estimate in the Simple Linear Regressionfor given potential centroids
###      in the Slope: different slope and different intercept
#### yList: The input data(in the List of Vectors Format)
#### XList: The input data(in the List of Matrices Format)
#### cdataL: Pre-defined cdata list which is from List.C.Data(DMatrix,MR)
#### ID: Indices for potential centroids
#### overlap: Boolean which is TRUE for overlapping clusters / FALSE for non-overlapping clusters
Fstat.CDL.ID.SL1.ST <- function(yList, XList, cdataL, ID, overlap) {
  mlc <- NULL
  T <- length(XList)
  N <- dim(XList[[1]])[1]; p <- dim(XList[[1]])[2]
  
  wholeID <- 1:N                             # Update: 2015/06/22
  nonID <- wholeID[!(wholeID %in% ID)]       # Update: 2015/06/22
  
  # compute sum(y_it), sum(y_it^2), sum(x_it), sum(x_it^2), sum(x_it*y_it) for each t=1,...T fot all i=1,...,N
  # Updated 2016/02/29
  sum_y1L <- list();  sum_y2L <- list();  sum_x1L <- list();  sum_x2L <- list();  sum_xyL <- list();
  n_na_y <- rep(0,T)      # number of is.na(y_t)  # Update: 2016/03/05
  for (t in 1:T) {
    y <- yList[[t]]
    x <- XList[[t]][,p]; x[is.na(y)] <- NA
    sum_y1L[[t]] <- sum(y, na.rm=TRUE);  sum_y2L[[t]] <- sum(y^2, na.rm=TRUE)
    sum_x1L[[t]] <- sum(x, na.rm=TRUE);  sum_x2L[[t]] <- sum(x^2, na.rm=TRUE)
    sum_xyL[[t]] <- sum(x*y, na.rm=TRUE)
    n_na_y[t] <- sum(is.na(y))  # Update: 2016/03/05
  }
  
  for (k in ID) {
    dat.k <- cdataL[[k]]
    dat.k$dr <- round(dat.k$r - c(dat.k$r[-1],0),10) # Update: 2015/08/05
    # For non-overlapping assumpsion
    # Update: 2015/06/22
    if (!overlap && length(which(dat.k$id %in% nonID)) > 0) {
      # remove all cluster which overlaps the previous clusters  
      # also remove all cluster which is greater than overlapped one
      dat.k <- dat.k[-(min(which(dat.k$id %in% nonID)):length(dat.k$id)),]
    }
    
    if (!is.null(dim(dat.k)) && dim(dat.k)[1] > p) {
      SSE1.kList <- SSE2.kList <- list()
      for (t in 1:T) {
        sum_y1 <- sum_y1L[[t]];  sum_y2 <- sum_y2L[[t]]
        sum_x1 <- sum_x1L[[t]];  sum_x2 <- sum_x2L[[t]]
        sum_xy <- sum_xyL[[t]]
        yc <- yList[[t]][dat.k$id]; xc <- XList[[t]][dat.k$id,p]
        
        length_yc <- length(yc)      # Updated 2016/02/29
        nc <- 1:length_yc            # Updated 2016/02/29
        if (sum(is.na(yc)) > 0) {    # Updated 2016/02/29
          na_yc <- (1:length_yc)[is.na(yc)]
          yc[na_yc] <- xc[na_yc] <- 0
          for (j in na_yc) {
            nc[j:length_yc] <- nc[j:length_yc] - 1
          }
        }
        
        # compute cumulative sum
        cs <- cum.sum.scalar(yc, xc)
        cs$nc <- nc                 # Updated 2016/02/29
        num_na<- n_na_y[t]          # Updated 2016/03/05
        
        # compute SSEs(t)
        detM <- (N-num_na)*(cs$nc*sum_x2-(cs$x1)^2) - cs$nc*sum_x1*(sum_x1-cs$x1) + cs$nc*(sum_x1*cs$x1-cs$nc*sum_x2)
        m11 <- cs$nc*sum_x2-(cs$x1)^2
        m12 <- -(cs$nc*(sum_x1-cs$x1))
        m13 <- sum_x1*cs$x1-cs$nc*sum_x2
        m22 <- (N-num_na)*cs$nc - cs$nc^2
        m23 <- -((N-num_na)*cs$x1 - cs$nc*sum_x1)
        m33 <- (N-num_na)*sum_x2 - (sum_x1)^2
        
        SSE1.kList[[t]] <- sum_y2 - (m11*(sum_y1)^2 + m22*(sum_xy)^2 + m33*(cs$y1)^2 + 
                                       2*m12*(sum_y1)*(sum_xy) + 2*m13*(sum_y1)*(cs$y1) + 2*m23*(sum_xy)*(cs$y1))/detM
        SSE2.kList[[t]] <- sum_y2 - (sum_y1-cs$y1)^2/(N-cs$nc-num_na) - (cs$y1)^2/(cs$nc) -
          ((sum_xy-cs$xy) - (sum_x1-cs$x1)*(sum_y1-cs$y1)/(N-cs$nc-num_na))^2/((sum_x2-cs$x2) - (sum_x1-cs$x1)^2/(N-cs$nc-num_na)) -
          ((cs$xy) - (cs$x1)*(cs$y1)/(cs$nc))^2/((cs$x2) - (cs$x1)^2/(cs$nc))   # Updated 2016/02/29
      }
      
      # compute SSEs and Fslope
      dat.k$SSE1 <- Reduce("+",SSE1.kList)
      dat.k$SSE2 <- Reduce("+",SSE2.kList)
      dat.k$Fslope <- ((dat.k$SSE1 - dat.k$SSE2)/(T))/(dat.k$SSE2/((N-4)*T-sum(n_na_y)))  # Updated 2016/03/05
      
      # Remove -Inf which is in case of all x's are the same
      dat.k2 <- as.matrix(dat.k[dat.k$dr!=0 & dat.k$Fslope > 0,])
      
      dim1 <- dim(dat.k2)[1]
      dim2 <- dim(dat.k2)[2]
      
      # Pick the largest Fslope. First p-1 is NA or NaN
      kpick <- order(-dat.k2[,dim2])[1]
      kpick1 <- c(k, dat.k2[kpick,3], dat.k2[kpick,(dim2-2):dim2])        # Updated 2016/03/05
      mlc <- rbind(mlc,kpick1)
    }
  }
  
  colnames(mlc) <- c("center","radius","SSE1","SSE2","maxFslope")        # Updated 2016/03/06
  d2 <- dim(mlc)[2]
  
  # Pick the largest maxFslope
  Mlc <- mlc[order(-mlc[,d2])[1],]
  
  # unbiased est of sigSq under Null Hypothesis (H1)
  sigSq1 <- Mlc[3]/((N-3)*T-sum(n_na_y))                             # Updated 2016/03/05
  names(sigSq1) <- NULL
  
  return(list(Mlc = Mlc[-(3:4)], mlc = mlc, sigSq1 = sigSq1))        # Updated 2016/03/05
}


###################################################################################################
### Added 2016/02/09
### Updated 2016/03/06
### Fstat.CDL.ID.SL2.ST(yList, XList, cdataL, ID, overlap):
###      find the Spatio-Temporal Cluster estimate in the Simple Linear Regressionfor given potential centroids
###      in the INTercept only: the same slope but different intercept
#### yList: The input data(in the List of Vectors Format)
#### XList: The input data(in the List of Matrices Format)
#### cdataL: Pre-defined cdata list which is from List.C.Data(DMatrix,MR)
#### ID: Indices for potential centroids
#### overlap: Boolean which is TRUE for overlapping clusters / FALSE for non-overlapping clusters
Fstat.CDL.ID.SL2.ST <- function(yList, XList, cdataL, ID, overlap) {
  mlc <- NULL
  T <- length(XList)
  N <- dim(XList[[1]])[1]; p <- dim(XList[[1]])[2]
  
  wholeID <- 1:N                             # Update: 2015/06/22
  nonID <- wholeID[!(wholeID %in% ID)]       # Update: 2015/06/22
  
  # compute sum(y_it), sum(y_it^2), sum(x_it), sum(x_it^2), sum(x_it*y_it) for each t=1,...T fot all i=1,...,N
  # Updated 2016/02/29
  sum_y1L <- list();  sum_y2L <- list();  sum_x1L <- list();  sum_x2L <- list();  sum_xyL <- list();
  n_na_y <- rep(0,T)      # number of is.na(y_t)  # Update: 2016/03/05
  for (t in 1:T) {
    y <- yList[[t]]
    x <- XList[[t]][,p]; x[is.na(y)] <- NA
    sum_y1L[[t]] <- sum(y, na.rm=TRUE);  sum_y2L[[t]] <- sum(y^2, na.rm=TRUE)
    sum_x1L[[t]] <- sum(x, na.rm=TRUE);  sum_x2L[[t]] <- sum(x^2, na.rm=TRUE)
    sum_xyL[[t]] <- sum(x*y, na.rm=TRUE)
    n_na_y[t] <- sum(is.na(y))  # Update: 2016/03/05
  }
  
  for (k in ID) {
    dat.k <- cdataL[[k]]
    dat.k$dr <- round(dat.k$r - c(dat.k$r[-1],0),10) # Update: 2015/08/05
    # For non-overlapping assumpsion
    # Update: 2015/06/22
    if (!overlap && length(which(dat.k$id %in% nonID)) > 0) {
      # remove all cluster which overlaps the previous clusters  
      # also remove all cluster which is greater than overlapped one
      dat.k <- dat.k[-(min(which(dat.k$id %in% nonID)):length(dat.k$id)),]
    }
    
    if (!is.null(dim(dat.k)) && dim(dat.k)[1] >= p) {  # the same slope / different intercept
      SSE1.kList <- list()
      for (t in 1:T) {
        sum_y1 <- sum_y1L[[t]];  sum_y2 <- sum_y2L[[t]]
        sum_x1 <- sum_x1L[[t]];  sum_x2 <- sum_x2L[[t]]
        sum_xy <- sum_xyL[[t]]
        yc <- yList[[t]][dat.k$id]; xc <- XList[[t]][dat.k$id,p]
        
        length_yc <- length(yc)      # Updated 2016/02/29
        nc <- 1:length_yc            # Updated 2016/02/29
        if (sum(is.na(yc)) > 0) {    # Updated 2016/02/29
          na_yc <- (1:length_yc)[is.na(yc)]
          yc[na_yc] <- xc[na_yc] <- 0
          for (j in na_yc) {
            nc[j:length_yc] <- nc[j:length_yc] - 1
          }
        }
        
        # compute cumulative sum
        cs <- cum.sum.scalar(yc, xc)
        cs$nc <- nc                 # Updated 2016/02/29
        num_na<- n_na_y[t]          # Updated 2016/03/05
        
        # compute SSE1(t)
        detM <- (N-num_na)*(cs$nc*sum_x2-(cs$x1)^2) - cs$nc*sum_x1*(sum_x1-cs$x1) + cs$nc*(sum_x1*cs$x1-cs$nc*sum_x2)
        m11 <- cs$nc*sum_x2-(cs$x1)^2
        m12 <- -(cs$nc*(sum_x1-cs$x1))
        m13 <- sum_x1*cs$x1-cs$nc*sum_x2
        m22 <- (N-num_na)*cs$nc - cs$nc^2
        m23 <- -((N-num_na)*cs$x1 - cs$nc*sum_x1)
        m33 <- (N-num_na)*sum_x2 - (sum_x1)^2
        
        SSE1.kList[[t]] <- sum_y2 - (m11*(sum_y1)^2 + m22*(sum_xy)^2 + m33*(cs$y1)^2 + 
                                       2*m12*(sum_y1)*(sum_xy) + 2*m13*(sum_y1)*(cs$y1) + 2*m23*(sum_xy)*(cs$y1))/detM
      }
      
      # compute SSE1
      dat.k$SSE1 <- Reduce("+",SSE1.kList)
      
      # Remove -Inf which is in case of all x's are the same
      dat.k2 <- as.matrix(dat.k[dat.k$dr!=0 & dat.k$SSE1 > 0,])
      
      dim1 <- dim(dat.k2)[1]
      dim2 <- dim(dat.k2)[2]
      
      # Pick the smallest one. First p-1 is NA or NaN
      kpick <- order(dat.k2[,dim2])[1]
      kpick1 <- c(k, dat.k2[kpick,3], dat.k2[kpick,dim2])
      mlc <- rbind(mlc,kpick1)
    }      
  }
  
  colnames(mlc) <- c("center","radius","SSE1")
  d2 <- dim(mlc)[2]
  
  # Pick the smallest one
  Mlc <- mlc[order(mlc[,d2])[1],]
  
  return(list(Mlc = Mlc, mlc = mlc))
}


###################################################################################################
### Added 2016/04/29
### Updated 2016/04/29
### TStg.CDL.ID.SL.ST(yList, XList, cdataL, ID, overlap):
###      find the Spatio-Temporal Cluster estimate in the Simple Linear Regressionfor given potential centroids
###      in the Slope (1st Stg) and the Intercept (2nd Stg). 
#### yList: The input data(in the List of Vectors Format)
#### XList: The input data(in the List of Matrices Format)
#### cdataL: Pre-defined cdata list which is from List.C.Data(DMatrix,MR)
#### ID: Indices for potential centroids
#### overlap: Boolean which is TRUE for overlapping clusters / FALSE for non-overlapping clusters
TStg.CDL.ID.SL.ST <- function(yList, XList, cdataL, ID, overlap) {
  mlc_slp <- mlc_int <- NULL
  T <- length(XList)
  N <- dim(XList[[1]])[1]; p <- dim(XList[[1]])[2]
  
  wholeID <- 1:N                             # Update: 2015/06/22
  nonID <- wholeID[!(wholeID %in% ID)]       # Update: 2015/06/22
  
  # compute sum(y_it), sum(y_it^2), sum(x_it), sum(x_it^2), sum(x_it*y_it) for each t=1,...T fot all i=1,...,N
  # Updated 2016/02/29
  sum_y1L <- list();  sum_y2L <- list();  sum_x1L <- list();  sum_x2L <- list();  sum_xyL <- list();
  n_na_y <- rep(0,T)      # number of is.na(y_t)  # Update: 2016/03/05
  for (t in 1:T) {
    y <- yList[[t]]
    x <- XList[[t]][,p]; x[is.na(y)] <- NA
    sum_y1L[[t]] <- sum(y, na.rm=TRUE);  sum_y2L[[t]] <- sum(y^2, na.rm=TRUE)
    sum_x1L[[t]] <- sum(x, na.rm=TRUE);  sum_x2L[[t]] <- sum(x^2, na.rm=TRUE)
    sum_xyL[[t]] <- sum(x*y, na.rm=TRUE)
    n_na_y[t] <- sum(is.na(y))  # Update: 2016/03/05
  }
  
  for (k in ID) {
    dat.k <- cdataL[[k]]
    dat.k$dr <- round(dat.k$r - c(dat.k$r[-1],0),10) # Update: 2015/08/05
    # For non-overlapping assumpsion
    # Update: 2015/06/22
    if (!overlap && length(which(dat.k$id %in% nonID)) > 0) {
      # remove all cluster which overlaps the previous clusters  
      # also remove all cluster which is greater than overlapped one
      dat.k <- dat.k[-(min(which(dat.k$id %in% nonID)):length(dat.k$id)),]
    }
    
    if (!is.null(dim(dat.k)) && dim(dat.k)[1] > p) {
      SSE1_slp.kList <- SSE1_int.kList <- SSE2.kList <- list()
      for (t in 1:T) {
        sum_y1 <- sum_y1L[[t]];  sum_y2 <- sum_y2L[[t]]
        sum_x1 <- sum_x1L[[t]];  sum_x2 <- sum_x2L[[t]]
        sum_xy <- sum_xyL[[t]]
        yc <- yList[[t]][dat.k$id]; xc <- XList[[t]][dat.k$id,p]
        
        length_yc <- length(yc)      # Updated 2016/02/29
        nc <- 1:length_yc            # Updated 2016/02/29
        if (sum(is.na(yc)) > 0) {    # Updated 2016/02/29
          na_yc <- (1:length_yc)[is.na(yc)]
          yc[na_yc] <- xc[na_yc] <- 0
          for (j in na_yc) {
            nc[j:length_yc] <- nc[j:length_yc] - 1
          }
        }
        
        # compute cumulative sum
        cs <- cum.sum.scalar(yc, xc)
        cs$nc <- nc                 # Updated 2016/02/29
        num_na<- n_na_y[t]          # Updated 2016/03/05
        
        # compute SSEs(t)
        detM_slp <- (N-num_na)*(cs$nc*sum_x2-(cs$x1)^2) - cs$nc*sum_x1*(sum_x1-cs$x1) + cs$nc*(sum_x1*cs$x1-cs$nc*sum_x2)
        m11_slp <- cs$nc*sum_x2-(cs$x1)^2
        m12_slp <- -(cs$nc*(sum_x1-cs$x1))
        m13_slp <- sum_x1*cs$x1-cs$nc*sum_x2
        m22_slp <- (N-num_na)*cs$nc - cs$nc^2
        m23_slp <- -((N-num_na)*cs$x1 - cs$nc*sum_x1)
        m33_slp <- (N-num_na)*sum_x2 - (sum_x1)^2
        
        detM_int <- (N-num_na)*(sum_x2-cs$x2)*(cs$x2) - (sum_x1)*(sum_x1-cs$x1)*(cs$x2) + (cs$x1)*(sum_x1*cs$x2-cs$x1*sum_x2)
        m11_int <- (sum_x2-cs$x2)*(cs$x2)
        m12_int <- -(sum_x1-cs$x1)*(cs$x2)
        m13_int <- (sum_x1*cs$x2-cs$x1*sum_x2)
        m22_int <- (N-num_na)*cs$x2 - (cs$x1)^2
        m23_int <- -((N-num_na)*cs$x2 - sum_x1*cs$x1)
        m33_int <- (N-num_na)*sum_x2 - (sum_x1)^2
        
        SSE1_slp.kList[[t]] <- sum_y2 - (m11_slp*(sum_y1)^2 + m22_slp*(sum_xy)^2 + m33_slp*(cs$y1)^2 + 
                                           2*m12_slp*(sum_y1)*(sum_xy) + 2*m13_slp*(sum_y1)*(cs$y1) + 2*m23_slp*(sum_xy)*(cs$y1))/detM_slp
        
        SSE1_int.kList[[t]] <- sum_y2 - (m11_int*(sum_y1)^2 + m22_int*(sum_xy)^2 + m33_int*(cs$xy)^2 + 
                                           2*m12_int*(sum_y1)*(sum_xy) + 2*m13_int*(sum_y1)*(cs$xy) + 2*m23_int*(sum_xy)*(cs$xy))/detM_int        
        
        SSE2.kList[[t]] <- sum_y2 - (sum_y1-cs$y1)^2/(N-cs$nc-num_na) - (cs$y1)^2/(cs$nc) -
          ((sum_xy-cs$xy) - (sum_x1-cs$x1)*(sum_y1-cs$y1)/(N-cs$nc-num_na))^2/((sum_x2-cs$x2) - (sum_x1-cs$x1)^2/(N-cs$nc-num_na)) -
          ((cs$xy) - (cs$x1)*(cs$y1)/(cs$nc))^2/((cs$x2) - (cs$x1)^2/(cs$nc))   # Updated 2016/02/29
      }
      
      # compute SSEs and Fslope
      dat.k$SSE1_slp <- Reduce("+",SSE1_slp.kList)
      dat.k$SSE1_int <- Reduce("+",SSE1_int.kList)
      dat.k$SSE2 <- Reduce("+",SSE2.kList)	  
      dat.k$F_slp <- ((dat.k$SSE1_slp - dat.k$SSE2)/(T))/(dat.k$SSE2/((N-4)*T-sum(n_na_y)))  # Updated 2016/04/29
      dat.k$F_int <- ((dat.k$SSE1_int - dat.k$SSE2)/(T))/(dat.k$SSE2/((N-4)*T-sum(n_na_y)))  # Updated 2016/04/29
      
      # Remove -Inf which is in case of all x's are the same
      dat.k2 <- as.matrix(dat.k[dat.k$dr!=0 & dat.k$SSE1_slp > 0,])
      
      dim1 <- dim(dat.k2)[1]
      dim2 <- dim(dat.k2)[2]
      
      # Pick the largest Fslp and Fint. First p-1 is NA or NaN
      kpick_slp <- order(-dat.k2[,(dim2-1)])[1]
      kpick1_slp <- c(k, dat.k2[kpick_slp,3], dat.k2[kpick_slp,c(5,7,(dim2-1))])
      mlc_slp <- rbind(mlc_slp,kpick1_slp)
      
      kpick_int <- order(-dat.k2[,dim2])[1]
      kpick1_int <- c(k, dat.k2[kpick_int,3], dat.k2[kpick_int,c(6,7,dim2)])
      mlc_int <- rbind(mlc_int,kpick1_int)
    }
  }
  
  colnames(mlc_slp) <- c("center","radius","SSE1_slp","SSE2","maxF_slp")        # Updated 2016/04/29
  colnames(mlc_int) <- c("center","radius","SSE1_int","SSE2","maxF_int")        # Updated 2016/04/29
  d2_slp <- dim(mlc_slp)[2]
  d2_int <- dim(mlc_int)[2]
  
  # Pick the largest maxF_slp and maxF_int
  Mlc_slp <- mlc_slp[order(-mlc_slp[,d2_slp])[1],]
  Mlc_int <- mlc_int[order(-mlc_int[,d2_int])[1],]
  
  # unbiased est of sigSq under Null Hypothesis (H1)
  sigSq1_slp <- Mlc_slp[3]/((N-3)*T-sum(n_na_y))                             # Updated 2016/04/29
  sigSq1_int <- Mlc_int[3]/((N-3)*T-sum(n_na_y))                             # Updated 2016/04/29
  
  names(sigSq1_slp) <- NULL
  names(sigSq1_int) <- NULL
  return(list(Mlc_slp = Mlc_slp[-(3:4)], mlc_slp = mlc_slp, sigSq1_slp = sigSq1_slp,
              Mlc_int = Mlc_int[-(3:4)], mlc_int = mlc_int, sigSq1_int = sigSq1_int))        # Updated 2016/04/29
}
