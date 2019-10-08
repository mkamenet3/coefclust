#'@title Test.MLC.SI.ID.ST
#'@description Find and test the cylindrical spatio-temporal cluster for given ID via simulataneous detection.
#'@param yList The input data (as a list of vectors).
#'@param XList The input data (as a list of matrices).
#'@param cdataL Pre-defined cdata list which is from \code{List.C.Data(DMatrix,MR)}.
#'@param M Number of simulations
#'@param ID Indices for potential centroids.
#'@param overlap  Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters
#'@return Most likely cluster, maximum F-statistic (of all simulations), and associated p-value.
Test.MLC.SI.ID.ST <- function(yList, XList, cdataL, M, ID, overlap) {
  TestStat <- rep(NA,M+1)
  N <- dim(XList[[1]])[1]; p <- dim(XList[[1]])[2]
  T <- length(XList)

  IPList <- list()
  SSE_0  <- 0
  n_na_y <- rep(0,T)      # number of is.na(y_t)
  for (t in 1:T) {
    y <- yList[[t]]
    Xt <- XList[[t]]
    y_valid <- y[!is.na(y)]
    Xt_valid <- Xt[!is.na(y),]
    Pt <- Xt_valid%*%solve(t(Xt_valid)%*%Xt_valid)%*%t(Xt_valid)
    IPList[[t]] <- diag(dim(Pt)[1]) - Pt
    SSE_0 <- SSE_0 + t(y_valid)%*%IPList[[t]]%*%(y_valid)
    n_na_y[t] <- sum(is.na(y))
  }
  df_den <- (N-2*p)*T - sum(n_na_y)    # df of the denominator in the F statistic

  Mlc <- MLC.CDL.ID.SL.ST(yList, XList, cdataL, ID, overlap)$Mlc
  l <- length(Mlc)
  SSE_2 <- Mlc[l]
  maxF <- ((SSE_0 - SSE_2)/(p*T))/(SSE_2/(df_den))
  sigSq0 <- SSE_0/((N-p)*T - sum(n_na_y))

  # F statistic with (df1,df2) = (p*T,df_den)
  TestStat[1]  <- maxF

  xList <- list();  sum_x1L <- list();  sum_x2L <- list();  cs.xxL.L <- list();
  for (t in 1:T) {
    y <- yList[[t]]
    x <- XList[[t]][,p]; x[is.na(y)] <- NA
    xList[[t]] <- x
    sum_x1L[[t]] <- sum(x, na.rm=TRUE);  sum_x2L[[t]] <- sum(x^2, na.rm=TRUE)
    cs.xxL <- list()
    for (i in 1:N) {
      dat.i <- cdataL[[i]]
      xc <- x[dat.i$id]
      length_xc <- length(xc)
      nc <- 1:length_xc
      if (sum(is.na(xc)) > 0) {
        na_xc <- (1:length_xc)[is.na(xc)]
        xc[na_xc] <- 0
        for (j in na_xc) {
          nc[j:length_xc] <- nc[j:length_xc] - 1
        }
      }
      cs <- cum.sum.scalar1(xc)
      cs$nc <- nc
      cs.xxL[[i]] <- cs
    }
    cs.xxL.L[[t]] <- cs.xxL
  }

  for (k in 1:M) {
    eList_k <- list()
    SSE_0 <- 0
    for (t in 1:T) {
      e_k <- rnorm(N, mean = 0, sd = sqrt(sigSq0))
      e_k[is.na(yList[[t]])] <- NA
      eList_k[[t]] <- e_k
      e_valid <- e_k[!is.na(e_k)]
      SSE_0 <- SSE_0 + t(e_valid)%*%IPList[[t]]%*%(e_valid)
    }
    SSE_2  <- MLC.CDL2.ID.SL.ST(eList_k, xList, cdataL, sum_x1L, sum_x2L, cs.xxL.L, ID, overlap)$Mlc[l]
    TestStat[k+1]<- ((SSE_0 - SSE_2)/(p*T))/(SSE_2/(df_den))
  }

  pval <- (rank(-TestStat)[1])/(M+1)
  c(Mlc, maxF=TestStat[1], pval=pval)
}


#'@title Find.Clusters.SI.ST
#'@description Find multiple cylindrical spatio-temporal clusters sequentially via simultaneous detection. Find and test the cluster in the simple linear regression for given potential centroids.
#'@param yList The input data (as a list of vectors).
#'@param XList The input data (as a list of matrices).
#'@param long longitude.
#'@param lat latitude.
#'@param MR Maximum radius.
#'@param M Number of simulations.
#'@param ID Indices for potential centroids.
#'@param overlap  Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters.
#'@param alpha Significance level
#'@return List of clusters, coefficients, and indicator of cluster membership.
#'@export
Find.Clusters.SI.ST <- function(yList, XList, long, lat, MR, M, overlap, alpha) {
  N <- dim(XList[[1]])[1]; p <- dim(XList[[1]])[2]
  T <- length(XList)
  ID <- 1:N
  LL <- cbind(long, lat)
  DMatrix <- geosphere::distm(LL)/1000
  cdataL <- List.C.Data(DMatrix,MR)

  print("Finding 1st Cluster")
  time_tmp <- system.time(C_tmp <- Test.MLC.SI.ID.ST(yList, XList, cdataL, M, ID, overlap))
  Clusters <- rbind(c(C_tmp,time_tmp[3]/60), NULL)    # Update: 2019/09/07
  pval_tmp <- C_tmp[5]
  cent_tmp <- C_tmp[1]

  clsL <- list()    # a list of indicator for each cluster
  n_cls <- 1        # number of clusters

  d_tmp <- geosphere::distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
  clsL[[n_cls]]  <- as.vector(d_tmp <= C_tmp[2])
  ID_tmp <- ID[ID != cent_tmp]

  yList_tmp <- list()
  for (t in 1:T) {
    yt <- yList[[t]]
    Xt <- XList[[t]]
    Xt_cls <- Xt*(clsL[[n_cls]])
    Wt <- cbind(Xt,Xt_cls)
    yt_valid <- yt[!is.na(yt)]
    Wt_valid <- Wt[!is.na(yt),]
    bt_tmp <- solve(t(Wt_valid)%*%Wt_valid)%*%t(Wt_valid)%*%yt_valid
    yList_tmp[[t]] <- yt - Xt_cls%*%bt_tmp[(p+1):(2*p)]
  }

  while (pval_tmp < alpha) {
    print(paste("Finding ", n_cls + 1, "th Cluster", sep=""))
    time_tmp <- system.time(C_tmp <- Test.MLC.SI.ID.ST(yList_tmp, XList, cdataL, M, ID_tmp, overlap))
    Clusters <- rbind(Clusters,c(C_tmp,time_tmp[3]/60))
    pval_tmp <- C_tmp[5]
    cent_tmp <- C_tmp[1]

    n_cls <- n_cls + 1

    d_tmp1 <- geosphere::distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
    d_tmp2 <- geosphere::distm(cbind(long[ID_tmp],lat[ID_tmp]), c(long[cent_tmp],lat[cent_tmp]))/1000
    clsL[[n_cls]]  <- as.vector((d_tmp1 <= C_tmp[2]))
    ID_tmp <- ID_tmp[ID_tmp != cent_tmp]

    for (t in 1:T) {
      yt <- yList_tmp[[t]]
      Xt <- XList[[t]]
      Xt_cls <- Xt*(clsL[[n_cls]])
      Wt <- cbind(Xt,Xt_cls)
      yt_valid <- yt[!is.na(yt)]
      Wt_valid <- Wt[!is.na(yt),]
      bt_tmp <- solve(t(Wt_valid)%*%Wt_valid)%*%t(Wt_valid)%*%yt_valid
      yList_tmp[[t]] <- yt - Xt_cls%*%bt_tmp[(p+1):(2*p)]
    }
  }

  n_obs <- sapply(clsL, sum)
  rownames(Clusters) <- NULL                   # Update: 2019/09/07
  Clusters <- cbind(Clusters, n_obs=n_obs)     # Update: 2019/09/07
  Clusters <- rbind(Clusters[,-c(3,4)], NULL)  # Update: 2019/09/07

  Coef <- Est.Coeff.SI.ST(yList, XList, long, lat, Clusters)  # Update: 2019/09/07

  return(list(Clusters = Clusters, Coef = Coef, Indicator = clsL))
}



#'@title Fit.Model.Clusters.ST
#'@description Fit a simple linear regression model with detected clusters.
#'@param yList The input data (as a list of vectors).
#'@param XList The input data (as a list of matrices).
#'@param long longitude.
#'@param lat latitude.
#'@param Clusters Inherited output from \code{Find.Clusters.SI.ST()}.
#'@export
Fit.Model.Clusters.ST <- function(yList, XList, long, lat, Clusters) {
  WList <- XList
  clsL <- list()
  n_cls <- dim(Clusters)[1] - 1       # number of clusters
  p <- dim(XList[[1]])[2]; T <- length(XList)

  coeff_namesList <- list()
  modelList <- list()
  for (t in 1:T) {
    coeff_namesList[[t]] <- paste("b0_", 0:(p-1), "_t", t, sep="")   # names for the parameter estimates
    modelList[[t]] <- paste(" + b0_", 0:(p-1), "_t", t, sep="", collapse = "")  # model to fit
  }

  if (n_cls > 0) {
    for (j in 1:n_cls) {
      for (t in 1:T) {
        coeff_namesList[[t]][(p*j+1):(p*j+p)] <- paste("c", j, "_", 0:(p-1), "_t", t, sep="")
        modelList[[t]] <- paste(modelList[[t]], paste(" + c", j, "_", 0:(p-1), "_t", t, sep="", collapse = ""), sep="")
        cent_tmp <- Clusters[j,1]
        r_tmp <- Clusters[j,2]
        d_tmp <- geosphere::distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
        clsL[[j]]  <- as.vector(d_tmp <= r_tmp)
        WList[[t]] <- cbind(WList[[t]],XList[[t]]*(clsL[[j]]))
      }
    }
  }

  coeff_names <- NULL
  model <- "y ~ -1"
  y <- NULL
  W <- 181818
  for (t in 1:T) {
    coeff_names <- c(coeff_names, coeff_namesList[[t]])
    model <- paste(model, modelList[[t]], sep="")
    y <- c(y,yList[[t]])
    W <- as.matrix(bdiag(W,WList[[t]]))
  }
  W <- W[-1,-1]

  data_cls <- data.frame(y,W)
  colnames(data_cls) <- c("y", coeff_names)
  lm(model, data = data_cls)
}


#'@title Est.Coeff.SI.ST
#'@description Estimate coefficients via simultaneous detection.
#'@param yList The input data (as a list of vectors).
#'@param XList The input data (as a list of matrices).
#'@param long longitude.
#'@param lat latitude.
#'@param Clusters Inherited output from \code{Find.Clusters.SI.ST()}.
#'@return List of coefficients
Est.Coeff.SI.ST <- function(yList, XList, long, lat, Clusters) {
  WList <- list()
  clsL <- list()
 # number of clusters
  n_cls <- dim(Clusters)[1]             # Update: 2019/09/07
  p <- dim(XList[[1]])[2]; T <- length(XList)

  b_tmp <- list()
  coef_tmp <- list()
  for (t in 1:T) {
    b_tmp[[t]] <- solve(t(XList[[t]])%*%XList[[t]])%*%t(XList[[t]])%*%yList[[t]]
    coef_tmp[[t]] <- c(b_tmp[[t]],rep(NA,length(b_tmp[[t]])))
  }

  XList_cls <- list()
  yList_tmp <- yList
  if (n_cls > 0) {
    for (j in 1:n_cls) {
      cent_tmp <- Clusters[j,1]
      r_tmp <- Clusters[j,2]
      d_tmp <- geosphere::distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
      clsL[[j]]  <- as.vector(d_tmp <= r_tmp)
      for (t in 1:T) {
        XList_cls[[t]] <- XList[[t]]*(clsL[[j]])
        WList[[t]] <- cbind(XList[[t]],XList_cls[[t]])
        b_tmp[[t]] <- solve(t(WList[[t]])%*%WList[[t]])%*%t(WList[[t]])%*%yList_tmp[[t]]
        coef_tmp[[t]] <- cbind(coef_tmp[[t]], b_tmp[[t]])
        yList_tmp[[t]] <- yList_tmp[[t]] - XList_cls[[t]]%*%b_tmp[[t]][(p+1):(2*p)]
      }
    }
  }

  Coef <- list()
  for (t in 1:T) {
    Coef[[t]] <- cbind(coef_tmp[[t]][1:p,1:n_cls],NULL)
    coef_names <- c("beta_0","beta_1")
    if (n_cls > 1) {                             # Update: 2019/09/07
      for (j in 1:(n_cls-1)) {
        Coef[[t]] <- rbind(Coef[[t]],
                           cbind(matrix(rep(NA,p*j),p,j),
                                 matrix(rep(coef_tmp[[t]][(p+1):(2*p),(j+1)],(n_cls-j)),p,(n_cls-j))))
        coef_names <- c(coef_names, paste("theta_", j, ",", 0:(p-1), sep=""))
      }
    }
    colnames(Coef[[t]]) <- 0:(n_cls-1)      # Update: 2019/09/07
    rownames(Coef[[t]]) <- coef_names
  }

  Coeff_Table <- NULL
  for (t in 1:length(Coef)) {
    Coeff_Table <- cbind(Coeff_Table, Coef[[t]][,(dim(Coef[[t]])[2])])
  }

  return(list(Coeff_History = Coef, Coeff_Table = Coeff_Table))
}

