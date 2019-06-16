
#'@title List.C.Data
#'@description Construct  a list of all possible clusters where r <= MR for each center.
#'@param DMatrix Distance matrix
#'@param MR Maximum radius.
#'@return Matrix
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


#'@title cum.sum.scalar
#'@description Compute cumulative sums.
#'@param yc y vector
#'@param xc x vector
#'@return a list
#'@export
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
#'@description TODO
#'@param xc x vector
#'@return a list
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
#'@description TODO
#'@param yc y vector
#'@param xc x vector
#'@return a list

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


#'@title DC.Simul.SL
#'@description Detect the Cluster in the Simple Linear Regression for given potential centroids.
#'@param y The input data(as a vector)
#'@param X The input data(as a matrix)
#'@param cdataL Pre-defined cdata list which is from List.C.Data(DMatrix,MR)
#'@param ID Indices for potential centroids
#'@param overlap Boolean which is TRUE for overlapping clusters / FALSE for non-overlapping clusters
#'@return List of most likely clusters
#'@export
#'

DC.Simul.SL <- function(y, X, cdataL, ID, overlap) {
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
    dat.k$dr <- round(dat.k$r - c(dat.k$r[-1],0),10)
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

      # Remove -Inf which is in case of all x's are the same
      dat.k2 <- as.matrix(dat.k[dat.k$dr!=0 & dat.k$Nsig2.hat > 0,])
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


#'@title DC.Simul.SL2
#'@description Detect the cluster in the simple linear regression for given potential centroids via simultaneous detection
#'@param y The input data(as a vector)
#'@param X The input data(as a matrix)
#'@param cdataL Pre-defined cdata list which is from List.C.Data(DMatrix,MR)
#'@param sum_x1 Pre-computed \eqn{\sum(x_i)}
#'@param sum_x2 Pre-computed \eqn{\sum(x_i^2)}
#'@param cs.xxL Pre-computed cs.xx list which is from \code{cum.sum.scalar1(xc)}
#'@param ID Indices for potential centroids
#'@param overlap  Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters
#'@return List of most likely clusters
#'@export
DC.Simul.SL2 <- function(y, X, cdataL, sum_x1, sum_x2, cs.xxL, ID, overlap) {
  mlc <- NULL
  N <- length(y); p <- 2

  # compute sum(y_i), sum(y_i^2), sum(x_i*y_i)
  sum_y1 <- sum(y);  sum_y2 <- sum(y^2)
  sum_xy <- sum(X*y)

  wholeID <- 1:N
  nonID <- wholeID[!(wholeID %in% ID)]

  for (k in ID) {
    dat.k <- cdataL[[k]]
    dat.k$dr <- round(dat.k$r - c(dat.k$r[-1],0),10)
    # For non-overlapping assumpsion
    if (!overlap && length(which(dat.k$id %in% nonID)) > 0) {
      # remove all cluster which overlaps the previous clusters
      # also remove all cluster which is greater than overlapped one
      dat.k <- dat.k[-(min(which(dat.k$id %in% nonID)):length(dat.k$id)),]
    }

    if (!is.null(dim(dat.k)) && dim(dat.k)[1] > p) {
      yc <- y[dat.k$id]; xc <- X[dat.k$id]

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
      dat.k2 <- as.matrix(dat.k[dat.k$dr!=0 & dat.k$Nsig2.hat > 0,])
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


#'@title DC.TStg1.SL
#'@description Detect the cluster in the simple linear regression for given potential centroids via the 1st stage in Two-stage detection: different slope and different intercept.
#'@param y The input data(as a vector)
#'@param X The input data(as a matrix)
#'@param cdataL Pre-defined cdata list which is from \code{List.C.Data(DMatrix,MR)}
#'@param ID Indices for potential centroids
#'@param overlap  Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters
#'@return List of most likely clusters
#'@export
DC.TStg1.SL <- function(y, X, cdataL, ID, overlap) {
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
    dat.k$dr <- round(dat.k$r - c(dat.k$r[-1],0),10)
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
  colnames(mlc) <- c("center","radius","sig2.hat","maxFstat.1")
  d2 <- dim(mlc)[2]

  # Pick the largest Fstat.1.
  Mlc <- mlc[order(-mlc[,d2])[1],]

  return(list(Mlc = Mlc, mlc = mlc))
}



#'@title DC.TStg2.SL
#'@description  Detect the cluster in the simple linear regression for given potential centroids via the 2nd stage in Two-stage detection: the same slope but different intercept.
#'@param y The input data(as a vector)
#'@param X The input data(as a matrix), N-by-2 matrix
#'@param cdataL Pre-defined cdata list which is from \code{List.C.Data(DMatrix,MR)}
#'@param ID Indices for potential centroids
#'@param overlap Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters
#'@return List of most likely clusters
#'@export

DC.TStg2.SL <- function(y, X, cdataL, ID, overlap) {
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
    dat.k$dr <- round(dat.k$r - c(dat.k$r[-1],0),10)
    # For non-overlapping assumpsion
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
      dat.k2 <- as.matrix(dat.k[dat.k$dr!=0 & dat.k$Nsig2.hat > 0,])
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

