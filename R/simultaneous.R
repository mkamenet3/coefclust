
#'@title Test.Cluster.Simul.SL
#'@description Find and test the cluster in the simple linear regression for given potential centroids via the simultaneous detection.
#'@param y The input data(as a vector)
#'@param X The input data(as a matrix)
#'@param cdataL Pre-defined cdata list which is from \code{List.C.Data(DMatrix,MR)}
#'@param M number of simulations
#'@param ID Indices for potential centroids
#'@param overlap Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters
#'@return Most likely cluster, maximum F-statistic (of all simulations), and associated p-value.
#'@export
Test.Cluster.Simul.SL <- function(y, X, cdataL, M, ID, overlap) {
  T <- rep(NA,M+1)
  N <- dim(X)[1]; p <- dim(X)[2]

  Mlc <- DC.Simul.SL(y, X, cdataL, ID, overlap)$Mlc
  l <- length(Mlc)

  P <- X%*%solve(t(X)%*%X)%*%t(X)
  IP <- diag(N) - P
  Ey <- P%*%y
  e_vec <- y - Ey
  sigSq0 <- t(e_vec)%*%e_vec/N

  # F statistic with (df1,df2) = (p,N-2p)
  T[1]  <- ((sigSq0 - Mlc[l])/p)/(Mlc[l]/(N-2*p))

  x <- X[,p]
  sum_x1 <- sum(x);  sum_x2 <- sum(x^2)
  cs.xxL <- list()
  for (i in 1:N) {
    dat.i <- cdataL[[i]]
    xc <- x[dat.i$id]
    cs.xxL[[i]] <- cum.sum.scalar1(xc)
  }

  for (k in 1:M) {
    e_k   <- rnorm(N, mean = 0, sd = sqrt(sigSq0))
    s2_k  <- DC.Simul.SL2(e_k, x, cdataL, sum_x1, sum_x2, cs.xxL, ID, overlap)$Mlc[l]
    T[k+1]<- ((t(e_k)%*%IP%*%e_k/N - s2_k)/p)/(s2_k/(N-2*p))
  }

  pval <- (rank(-T)[1])/(M+1)
  c(Mlc, maxF=T[1], pval=pval)
}


#'@title Find.Clusters.Simul
#'@description Find multiple clusters sequentially via simulataneous detection.
#'Find and test the cluster in the simple linear regression for given potential centroids via the simultaneous detection.
#'@param y The input data(as a vector)
#'@param X The input data(as a matrix)
#'@param long longitude
#'@param lat latitude
#'@param MR Maximum radius
#'@param M number of simulations
#'@param overlap  Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters
#'@param alpha significance level
#'@return list of cluster, coefficient, and indicator of cluster membership.
#'@export
Find.Clusters.Simul <- function(y, X, long, lat, MR, M, overlap, alpha) {
  ID <- 1:length(y)
  N <- dim(X)[1]; p <- dim(X)[2]
  LL <- cbind(long, lat)
  DMatrix <- geosphere::distm(LL)/1000
  cdataL <- List.C.Data(DMatrix,MR)

  b_tmp <- solve(t(X)%*%X)%*%t(X)%*%y
  coef_tmp <- c(b_tmp,rep(NA,length(b_tmp)))

  message("Finding 1st Cluster")
  time_tmp <- system.time(C_tmp <- Test.Cluster.Simul.SL(y, X, cdataL, M, ID, overlap))
  # Clusters <- c(C_tmp,time_tmp[3]/60)
  Clusters <- rbind(c(C_tmp,time_tmp[3]/60), NULL)    # Update: 2019/06/17
  pval_tmp <- C_tmp[5]
  cent_tmp <- C_tmp[1]

  clsL <- list()    # a list of indicator for each cluster
  n_cls <- 1        # number of clusters

  d_tmp <- geosphere::distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
  clsL[[n_cls]]  <- as.vector(d_tmp <= C_tmp[2])
  if (overlap) {
    ID_tmp <- ID[ID != cent_tmp]
  } else {
    ID_tmp <- ID[!(d_tmp <= C_tmp[2])]
  }

  X_cls <- X*(clsL[[n_cls]])
  W <- cbind(X,X_cls)
  b_tmp <- solve(t(W)%*%W)%*%t(W)%*%y
  coef_tmp <- cbind(coef_tmp,b_tmp)
  y_tmp <- y - X_cls%*%b_tmp[(p+1):(2*p)]

  while (pval_tmp < alpha) {
    message(paste("Finding ", n_cls + 1, "th Cluster", sep=""))
    time_tmp <- system.time(C_tmp <- Test.Cluster.Simul.SL(y_tmp, X, cdataL, M, ID_tmp, overlap))
    Clusters <- rbind(Clusters,c(C_tmp,time_tmp[3]/60))
    pval_tmp <- C_tmp[5]
    cent_tmp <- C_tmp[1]

    n_cls <- n_cls + 1

    d_tmp1 <- geosphere::distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
    d_tmp2 <- geosphere::distm(cbind(long[ID_tmp],lat[ID_tmp]), c(long[cent_tmp],lat[cent_tmp]))/1000
    clsL[[n_cls]]  <- as.vector((d_tmp1 <= C_tmp[2]))
    if (overlap) {
      ID_tmp <- ID_tmp[ID_tmp != cent_tmp]
    } else {
      ID_tmp <- ID_tmp[!(d_tmp2 <= C_tmp[2])]
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
  Coef <- cbind(coef_tmp[1:p,1:n_cls],NULL)    # Update: 2019/06/17
  coef_names <- c("beta_0","beta_1")
  if (n_cls > 1) {                             # Update: 2019/06/17
    for (j in 1:(n_cls-1)) {
      Coef <- rbind(Coef,
                    cbind(matrix(rep(NA,p*j),p,j),
                          matrix(rep(coef_tmp[(p+1):(2*p),(j+1)],(n_cls-j)),p,(n_cls-j))))
      coef_names <- c(coef_names, paste("theta_", j, ",", 0:(p-1), sep=""))
    }
  }
  colnames(Coef) <- 0:(n_cls-1)
  rownames(Coef) <- coef_names
  # Clusters <- Clusters[,-c(3,4)]
  Clusters <- rbind(Clusters[,-c(3,4)], NULL)  # Update: 2019/06/17
  return(list(Clusters = Clusters, Coef = Coef, clsL = clsL))
}
