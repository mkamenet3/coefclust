
#'@title Test.Cluster.TStg1.SL
#'@description Find and test the cluster in the simple linear regression for given potential centroids via the 1st stage in Two-stage detection: different slope and different intercept.
#'@param y The input data(as a vector)
#'@param X The input data(as a matrix)
#'@param cdataL Pre-defined cdata list which is from \code{List.C.Data(DMatrix,MR)}
#'@param M number of simulations
#'@param ID Indices for potential centroids
#'@param overlap Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters
#'@return Most likely cluster, maximum F-statistic (of all simulations), and associated p-value.
#'@export
Test.Cluster.TStg1.SL <- function(y, X, cdataL, M, ID, overlap) {
  T <- rep(NA,M+1)
  N <- dim(X)[1]; p <- dim(X)[2]

  Mlc <- DC.TStg1.SL(y, X, cdataL, ID, overlap)$Mlc
  l <- length(Mlc)

  # F statistic of the Slope with (df1,df2) = (1,N-2p)
  T[1]  <- Mlc[l]

  P <- X%*%solve(t(X)%*%X)%*%t(X)
  IP <- diag(N) - P
  Ey <- P%*%y
  e_vec <- y - Ey
  sigSq0 <- t(e_vec)%*%e_vec/N

  for (k in 1:M) {
    e_k   <- rnorm(N, mean = 0, sd = sqrt(sigSq0))
    Fstat_k <- DC.TStg1.SL(e_k, X, cdataL, ID, overlap)$Mlc[l]
    T[k+1]<- Fstat_k
  }

  pval <- (rank(-T)[1])/(M+1)
  c(Mlc, pval=pval)
}


#'@title Find.Clusters.TStg1
#'@description Find multiple clusters sequentially via simulataneous detection.
#'Find and test the cluster in the simple linear regression for given potential centroids via
#'the 1st stage in Two-stage detection: different slope and different intercept.
#'@param y The input data(as a vector)
#'@param X The input data(as a matrix)
#'@param long longitude
#'@param lat latitude
#'@param MR Maximum radius
#'@param M number of simulations
#'@param overlap  Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters
#'@param alpha significance level
#'@return list of cluster, coefficient
#'@export
Find.Clusters.TStg1 <- function(y, X, long, lat, MR, M, overlap, alpha) {
  ID <- 1:length(y)
  N <- dim(X)[1]; p <- dim(X)[2]
  LL <- cbind(long, lat)
  DMatrix <- distm(LL)/1000
  cdataL <- List.C.Data(DMatrix,MR)

  b_tmp <- solve(t(X)%*%X)%*%t(X)%*%y
  coef_tmp <- c(b_tmp,rep(NA,length(b_tmp)))

  print("Finding 1st Cluster")
  time_tmp <- system.time(C_tmp <- Test.Cluster.TStg1.SL(y, X, cdataL, M, ID, overlap))
  # Clusters <- c(C_tmp,time_tmp[3]/60)
  Clusters <- rbind(c(C_tmp,time_tmp[3]/60), NULL)    # Update: 2019/06/16
  pval_tmp <- C_tmp[5]
  cent_tmp <- C_tmp[1]
  r_tmp    <- C_tmp[2]

  clsL <- list()    # a list of indicator for each cluster
  n_cls <- 1        # number of clusters

  d_tmp <- distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
  clsL[[n_cls]]  <- as.vector(d_tmp <= r_tmp)
  if (overlap) {
    ID_tmp <- ID[ID != cent_tmp]
  } else {
    ID_tmp <- ID[!(d_tmp <= r_tmp)]
  }

  X_cls <- X*(clsL[[n_cls]])
  W <- cbind(X,X_cls)
  b_tmp <- solve(t(W)%*%W)%*%t(W)%*%y
  coef_tmp <- cbind(coef_tmp,b_tmp)
  y_tmp <- y - X_cls%*%b_tmp[(p+1):(2*p)]

  while (pval_tmp < alpha) {
    print(paste("Finding ", n_cls + 1, "th Cluster", sep=""))
    time_tmp <- system.time(C_tmp <- Test.Cluster.TStg1.SL(y_tmp, X, cdataL, M, ID_tmp, overlap))
    Clusters <- rbind(Clusters,c(C_tmp,time_tmp[3]/60))
    pval_tmp <- C_tmp[5]
    cent_tmp <- C_tmp[1]
    r_tmp    <- C_tmp[2]

    n_cls <- n_cls + 1

    d_tmp1 <- distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
    d_tmp2 <- distm(cbind(long[ID_tmp],lat[ID_tmp]), c(long[cent_tmp],lat[cent_tmp]))/1000
    clsL[[n_cls]]  <- as.vector(d_tmp1 <= r_tmp)
    if (overlap) {
      ID_tmp <- ID_tmp[ID_tmp != cent_tmp]
    } else {
      ID_tmp <- ID_tmp[!(d_tmp2 <= r_tmp)]
    }

    X_cls <- X*(clsL[[n_cls]])
    W <- cbind(X,X_cls)
    b_tmp <- solve(t(W)%*%W)%*%t(W)%*%y_tmp
    coef_tmp <- cbind(coef_tmp,b_tmp)
    y_tmp <- y_tmp - X_cls%*%b_tmp[(p+1):(2*p)]
  }

  n_obs <- sapply(clsL, sum)
  # if (length(n_obs) > 1) {
  #   rownames(Clusters) <- NULL
  #   Clusters <- cbind(Clusters, n_obs=n_obs)
  # } else {
  #   Clusters <- c(Clusters, n_obs=n_obs)
  # }
  rownames(Clusters) <- NULL                   # Update: 2019/06/17
  Clusters <- cbind(Clusters, n_obs=n_obs)     # Update: 2019/06/17

  # Coef <- coef_tmp[1:p,1:n_cls]
  Coef <- cbind(coef_tmp[1:p,1:n_cls],NULL)    # Update: 2019/06/16
  coef_names <- c("beta_0","beta_1")
  if (n_cls > 1) {                             # Update: 2019/06/16
    for (j in 1:(n_cls-1)) {
      Coef <- rbind(Coef,
                    cbind(matrix(rep(NA,p*j),p,j),
                          matrix(rep(coef_tmp[(p+1):(2*p),(j+1)],(n_cls-j)),p,(n_cls-j))))
      coef_names <- c(coef_names, paste("theta_", j, ",", 0:(p-1), sep=""))
    }
  }
  colnames(Coef) <- 0:(n_cls-1)
  rownames(Coef) <- coef_names

  return(list(Clusters = Clusters, Coef = Coef, clsL = clsL))
}


#'@title Test.Cluster.TStg2.SL
#'@description Find and test the cluster in the simple linear regression for given potential centroids via the 2nd stage in Two-stage detection: the same slope but different intercept.
#'@param y The input data(as a vector)
#'@param X The input data(as a matrix)
#'@param cdataL Pre-defined cdata list which is from \code{List.C.Data(DMatrix,MR)}
#'@param M number of simulations
#'@param ID Indices for potential centroids
#'@param overlap Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters
#'@return Most likely cluster, maximum F-statistic (of all simulations), and associated p-value.
#'@export

Test.Cluster.TStg2.SL <- function(y, X, cdataL, M, ID, overlap) {
  T <- rep(NA,M+1)
  N <- dim(X)[1]; p <- dim(X)[2]
  Mlc <- DC.TStg2.SL(y, X, cdataL, ID, overlap)$Mlc
  l <- length(Mlc)

  P <- X%*%solve(t(X)%*%X)%*%t(X)
  IP <- diag(N) - P
  Ey <- P%*%y
  e_vec <- y - Ey
  sigSq0 <- t(e_vec)%*%e_vec/N

  # F statistic with (df1,df2) = (1,N-3)
  T[1]  <- ((sigSq0 - Mlc[l])/1)/(Mlc[l]/(N-3))

  for (k in 1:M) {
    e_k   <- rnorm(N, mean = 0, sd = sqrt(sigSq0))
    s2_k <- DC.TStg2.SL(e_k, X, cdataL, ID, overlap)$Mlc[l]
    T[k+1]<- ((t(e_k)%*%IP%*%e_k/N - s2_k)/1)/(s2_k/(N-3))
  }

  pval <- (rank(-T)[1])/(M+1)
  c(Mlc, maxFstat.0=T[1], pval=pval)
}


#'@title Find.Clusters.TStg2
#'@description Find multiple clusters sequentially via SImulataneous detection.
#' Find and test the cluster in the simple linear regression for given potential centroids via the 2nd stage in Two-stage detection: the same slope but different intercept.
#'@param y The input data(as a vector)
#'@param X The input data(as a matrix)
#'@param long longitude
#'@param lat latitude
#'@param MR Maximum radius
#'@param M number of simulations
#'@param Cls1st the output from \code{Find.Clusters.TStg1} in the 1st-Stage
#'@param overlap  Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters
#'@param alpha significance level
#'@return list of cluster, coefficient

Find.Clusters.TStg2 <- function(y, X, long, lat, MR, M, Cls1st, overlap, alpha) {
  ID <- 1:length(y)
  N <- dim(X)[1]; p <- dim(X)[2]
  LL <- cbind(long, lat)
  DMatrix <- distm(LL)/1000
  cdataL <- List.C.Data(DMatrix,MR)

  clsL <- list()
  n_cls <- dim(Cls1st$Clusters)[1] - 1   # number of clusters
  n_cls1<- n_cls                         # number of significant clusters in the 1st-Stage

  ID_tmp <- ID
  y_tmp <- y
  W <- X
  if (n_cls1 > 0) {                            # Update: 2019/06/17
    for (j in 1:n_cls1) {
      cent_tmp <- Cls1st$Clusters[j,1]
      r_tmp <- Cls1st$Clusters[j,2]
      d_tmp1 <- distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
      d_tmp2 <- distm(cbind(long[ID_tmp],lat[ID_tmp]), c(long[cent_tmp],lat[cent_tmp]))/1000

      clsL[[j]]  <- as.vector(d_tmp1 <= r_tmp)
      if (overlap) {
        ID_tmp <- ID_tmp[ID_tmp != cent_tmp]
      } else {
        ID_tmp <- ID_tmp[!(d_tmp2 <= r_tmp)]
      }
      W <- cbind(W,X*(clsL[[j]]))
    }
    y_tmp <- y_tmp - W[,-(1:p)]%*%Cls1st$Coef[-(1:p),(n_cls1+1)]
  }

  pval_tmp <- 0
  Clusters <- NULL
  coef_tmp <- NULL
  while (pval_tmp < alpha) {
    print(paste("Finding ", n_cls + 1, "th Cluster", sep=""))
    time_tmp <- system.time(C_tmp <- Test.Cluster.TStg2.SL(y_tmp, X, cdataL, M, ID_tmp, overlap))
    Clusters <- rbind(Clusters,c(C_tmp,time_tmp[3]/60))
    pval_tmp <- C_tmp[5]
    cent_tmp <- C_tmp[1]
    r_tmp    <- C_tmp[2]

    n_cls <- n_cls + 1

    d_tmp1 <- distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
    d_tmp2 <- distm(cbind(long[ID_tmp],lat[ID_tmp]), c(long[cent_tmp],lat[cent_tmp]))/1000
    clsL[[n_cls]]  <- as.vector(d_tmp1 <= r_tmp)
    if (overlap) {
      ID_tmp <- ID_tmp[ID_tmp != cent_tmp]
    } else {
      ID_tmp <- ID_tmp[!(d_tmp2 <= r_tmp)]
    }

    X_cls <- 1*(clsL[[n_cls]])
    W <- cbind(X,X_cls)
    b_tmp <- solve(t(W)%*%W)%*%t(W)%*%y_tmp
    coef_tmp <- cbind(coef_tmp,b_tmp)
    y_tmp <- y_tmp - X_cls*b_tmp[(p+1)]
  }

  # n_obs <- sapply(clsL, sum)[-(1:n_cls1)]
  n_obs <- sapply(clsL, sum)                    # Update: 2019/06/17
  if (n_cls1 > 0) {                             # Update: 2019/06/17
    n_obs <- n_obs[-(1:n_cls1)]
  }

  # if (length(n_obs) > 1) {
  #   rownames(Clusters) <- NULL
  #   Clusters <- cbind(Clusters, n_obs=n_obs)
  # } else {
  #   Clusters <- c(Clusters, n_obs=n_obs)
  # }
  rownames(Clusters) <- NULL                   # Update: 2019/06/17
  Clusters <- cbind(Clusters, n_obs=n_obs)     # Update: 2019/06/17

  n_cls2 <- n_cls - n_cls1 - 1
  if (n_cls2 > 0) {                            # Update: 2019/06/17
    Coef1 <- cbind(Cls1st$Coef,
                   rbind(as.matrix(coef_tmp[1:p,1:n_cls2]),
                         matrix(rep(Cls1st$Coef[-(1:p),(n_cls1+1)],n_cls2),(n_cls1*p),n_cls2)))
    Coef2 <- NULL
    coef_names <- NULL
    for (j in 1:n_cls2) {
      Coef2 <- rbind(Coef2, c(rep(NA,(n_cls1+j)), rep(coef_tmp[(p+1),j],(n_cls2-j+1))))
      coef_names <- c(coef_names, paste("theta_", (n_cls1+j), ",", 0, sep=""))
    }
    rownames(Coef2) <- coef_names
    Coef <- rbind(Coef1,Coef2)
  } else {
    Coef <- Cls1st$Coef
  }
  colnames(Coef) <- 0:(n_cls-1)

  if (n_cls2 > 0) {                            # Update: 2019/06/17
    if (n_cls1 > 0) {
      Clusters <- rbind(Cls1st$Clusters[1:n_cls1,], Clusters)
    }
    Clusters <- cbind(Clusters, stage=c(rep(1,n_cls1),rep(2,(n_cls-n_cls1))))
  } else {
    Clusters <- cbind(Cls1st$Clusters, stage=1)
  }

  # Clusters <- Clusters[,-c(3,4)]
  Clusters <- rbind(Clusters[,-c(3,4)], NULL)  # Update: 2019/06/17

  return(list(Clusters = Clusters, Coef = Coef, clsL = clsL))
}

