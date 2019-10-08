#'@title Test.MLC.TS.ID.ST1
#'@description Find and test the cylindrical spatio-temporal cluster in the slopes for given ID via two-stage detection (in the 1st Stage).
#'@param yList The input data (as a list of vectors).
#'@param XList The input data (as a list of matrices).
#'@param cdataL Pre-defined cdata list which is from \code{List.C.Data(DMatrix,MR)}.
#'@param M Number of simulations
#'@param ID Indices for potential centroids.
#'@param overlap  Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters
#'@return Most likely cluster, maximum F-statistic (of all simulations), and associated p-value.
#'
Test.MLC.TS.ID.ST1 <- function(yList, XList, cdataL, M, ID, overlap) {
  TestStat <- rep(NA,M+1)
  N <- dim(XList[[1]])[1]; p <- dim(XList[[1]])[2]
  T <- length(XList)

  MlcSlope <- TStg.CDL.ID.SL.ST(yList, XList, cdataL, ID, overlap)
  Mlc <- MlcSlope$Mlc_slp
  l <- length(Mlc)
  sigSq1 <- MlcSlope$sigSq1_slp
  # F statistic for the slope
  TestStat[1]  <- MlcSlope$Mlc_slp[l]

  cent_tmp <- Mlc[1]
  r_tmp <- Mlc[2]
  d_tmp1 <- geosphere::distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
  cluster  <- as.vector((d_tmp1 <= r_tmp))

  EyList <- list()
  for (t in 1:T) {
    yt <- yList[[t]]
    Wt <- cbind(XList[[t]],cluster*1)
    yt_valid <- yt[!is.na(yt)]
    Wt_valid <- Wt[!is.na(yt),]
    bt_tmp <- solve(t(Wt_valid)%*%Wt_valid)%*%t(Wt_valid)%*%yt_valid
    EyList[[t]] <- Wt%*%bt_tmp
  }

  for (k in 1:M) {
    yList_k <- list()
    for (t in 1:T) {
      e_k <- stats::rnorm(N, mean = 0, sd = sqrt(sigSq1))
      e_k[is.na(yList[[t]])] <- NA
      yList_k[[t]] <- EyList[[t]] + e_k
    }
    TestStat[k+1] <- TStg.CDL.ID.SL.ST(yList_k, XList, cdataL, ID, overlap)$Mlc_slp[l]
  }

  pval <- (rank(-TestStat)[1])/(M+1)
  c(Mlc, pval=pval)
}


#'@title Find.Clusters.TS.ST1
#'@description Find multiple (overlapping) cylindrical spatio-temporal clusters sequentially in the slopes via two-stage detection (in the 1st stage).
#'@param yList The input data (as a list of vectors).
#'@param XList The input data (as a list of matrices).
#'@param long longitude.
#'@param lat latitude.
#'@param MR Maximum radius.
#'@param M Number of simulations.
#'@param overlap  Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters.
#'@param alpha Significance level
#'@return List of clusters, coefficients, and indicator of cluster membership.
#'@export
#'
Find.Clusters.TS.ST1 <- function(yList, XList, long, lat, MR, M, overlap, alpha) {
  N <- dim(XList[[1]])[1]; p <- dim(XList[[1]])[2]
  T <- length(XList)
  ID <- 1:N
  LL <- cbind(long, lat)
  DMatrix <- geosphere::distm(LL)/1000
  cdataL <- List.C.Data(DMatrix,MR)

  print("Finding 1st Cluster")
  time_tmp <- system.time(C_tmp <- Test.MLC.TS.ID.ST1(yList, XList, cdataL, M, ID, overlap))
  Clusters <- rbind(c(C_tmp,time_tmp[3]/60), NULL)    # Update: 2019/09/08

  l <- length(C_tmp)
  pval_tmp <- C_tmp[l]
  cent_tmp <- C_tmp[1]
  r_tmp    <- C_tmp[2]

  clsL <- list()    # a list of indicator for each cluster
  n_cls <- 1        # number of clusters

  d_tmp <- geosphere::distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
  clsL[[n_cls]]  <- as.vector(d_tmp <= r_tmp)
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
    time_tmp <- system.time(C_tmp <- Test.MLC.TS.ID.ST1(yList_tmp, XList, cdataL, M, ID_tmp, overlap))
    Clusters <- rbind(Clusters,c(C_tmp,time_tmp[3]/60))
    pval_tmp <- C_tmp[l]
    cent_tmp <- C_tmp[1]
    r_tmp    <- C_tmp[2]

    n_cls <- n_cls + 1

    d_tmp1 <- geosphere::distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
    d_tmp2 <- geosphere::distm(cbind(long[ID_tmp],lat[ID_tmp]), c(long[cent_tmp],lat[cent_tmp]))/1000
    clsL[[n_cls]]  <- as.vector((d_tmp1 <= r_tmp))
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
  rownames(Clusters) <- NULL                   # Update: 2019/09/08
  Clusters <- cbind(Clusters, n_obs=n_obs)     # Update: 2019/09/08
  colnames(Clusters)[3] <- "TestStat"          # Update: 2019/09/08
  # Update: 2019/09/09
  Coef <- Est.Coeff.TS.ST(yList, XList, Cls1stIndicator=clsL, Cls2ndIndicator=clsL)

  return(list(Clusters = Clusters, Coef = Coef, Indicator = clsL))
}



#'@title Test.MLC.TS.ID.ST2
#'@description Find and test the cylindrical spatio-temporal cluster in the intercepts for given ID via two-stages detection (in the 2nd Stage).
#'@param yList The input data (as a list of vectors).
#'@param XList The input data (as a list of matrices).
#'@param cdataL Pre-defined cdata list which is from \code{List.C.Data(DMatrix,MR)}.
#'@param M Number of simulations
#'@param ID Indices for potential centroids.
#'@param overlap  Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters
#'@return Most likely cluster, maximum F-statistic (of all simulations), and associated p-value.
Test.MLC.TS.ID.ST2 <- function(yList, XList, cdataL, M, ID, overlap) {
  TestStat <- rep(NA,M+1)
  N <- dim(XList[[1]])[1]; p <- dim(XList[[1]])[2]
  T <- length(XList)

  MlcSlope <- TStg.CDL.ID.SL.ST(yList, XList, cdataL, ID, overlap)
  Mlc <- MlcSlope$Mlc_int
  l <- length(Mlc)
  sigSq1 <- MlcSlope$sigSq1_int
  # F statistic for the slope
  TestStat[1]  <- MlcSlope$Mlc_int[l]

  cent_tmp <- Mlc[1]
  r_tmp <- Mlc[2]
  d_tmp1 <- geosphere::distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
  cluster  <- as.vector((d_tmp1 <= r_tmp))

  EyList <- list()
  for (t in 1:T) {
    yt <- yList[[t]]
    Wt <- cbind(XList[[t]],XList[[t]][,2]*cluster)
    yt_valid <- yt[!is.na(yt)]
    Wt_valid <- Wt[!is.na(yt),]
    bt_tmp <- solve(t(Wt_valid)%*%Wt_valid)%*%t(Wt_valid)%*%yt_valid
    EyList[[t]] <- Wt%*%bt_tmp
  }

  for (k in 1:M) {
    yList_k <- list()
    for (t in 1:T) {
      e_k <- stats::rnorm(N, mean = 0, sd = sqrt(sigSq1))
      e_k[is.na(yList[[t]])] <- NA
      yList_k[[t]] <- EyList[[t]] + e_k
    }
    TestStat[k+1] <- TStg.CDL.ID.SL.ST(yList_k, XList, cdataL, ID, overlap)$Mlc_int[l]
  }

  pval <- (rank(-TestStat)[1])/(M+1)
  c(Mlc, pval=pval)
}


#'@title Find.Clusters.TS.ST2
#'@description Find multiple (overlapping) cylindrical spatio-temporal clusters sequentially in the intercepts via two-stages detection (in the 2nd stage)
#'@param yList The input data (as a list of vectors).
#'@param XList The input data (as a list of matrices).
#'@param long longitude.
#'@param lat latitude.
#'@param MR Maximum radius.
#'@param M Number of simulations.
#'@param Cls1st Output from \code{Find.Clusters.TS.ST1()} function in the 1st Stage.
#'@param overlap  Boolean which is \code{TRUE} for overlapping clusters / \code{FALSE} for non-overlapping clusters.
#'@param alpha Significance level
#'@return List of clusters, coefficients, and indicator of cluster membership.
#'@export
#'
Find.Clusters.TS.ST2 <- function(yList, XList, long, lat, MR, M, Cls1st, overlap, alpha) {
  N <- dim(XList[[1]])[1]; p <- dim(XList[[1]])[2]
  T <- length(XList)
  ID <- 1:N
  LL <- cbind(long, lat)
  DMatrix <- geosphere::distm(LL)/1000
  cdataL <- List.C.Data(DMatrix,MR)

  clsL <- Cls1st$Indicator
  n_cls <- length(clsL) - 1       # number of clusters
  n_cls1<- n_cls                  # number of significant clusters in the 1st-Stage / Update: 2015/06/19

  ID_tmp <- ID
  yList_tmp <- yList
  WList <- XList
  if (n_cls1 > 0) {
    for (j in 1:n_cls1) {
      cent_tmp <- Cls1st$Clusters[j,1]
      r_tmp <- Cls1st$Clusters[j,2]

      d_tmp1 <- geosphere::distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
      d_tmp2 <- geosphere::distm(cbind(long[ID_tmp],lat[ID_tmp]), c(long[cent_tmp],lat[cent_tmp]))/1000
      ID_tmp <- ID_tmp[ID_tmp != cent_tmp]
      for (t in 1:T) {
        WList[[t]] <- cbind(WList[[t]],XList[[t]]*(clsL[[j]]))
      }
    }

    for (t in 1:T) {
      yt <- yList_tmp[[t]]
      Wt <- WList[[t]]
      yt_valid <- yt[!is.na(yt)]
      Wt_valid <- Wt[!is.na(yt),]
      bt_tmp <- solve(t(Wt_valid)%*%Wt_valid)%*%t(Wt_valid)%*%yt_valid
      yList_tmp[[t]] <- yt - Wt[,(p+1):length(bt_tmp)]%*%bt_tmp[(p+1):length(bt_tmp)]
    }
  }

  pval_tmp <- 0
  Clusters <- NULL
  while (pval_tmp < alpha) {
    print(paste("Finding ", n_cls + 1, "th Cluster", sep=""))
    time_tmp <- system.time(C_tmp <- Test.MLC.TS.ID.ST2(yList_tmp, XList, cdataL, M, ID_tmp, overlap))
    Clusters <- rbind(Clusters,c(C_tmp,time_tmp[3]/60))
    l <- length(C_tmp)
    pval_tmp <- C_tmp[l]
    cent_tmp <- C_tmp[1]
    r_tmp    <- C_tmp[2]

    n_cls <- n_cls + 1

    d_tmp1 <- geosphere::distm(cbind(long,lat), c(long[cent_tmp],lat[cent_tmp]))/1000
    d_tmp2 <- geosphere::distm(cbind(long[ID_tmp],lat[ID_tmp]), c(long[cent_tmp],lat[cent_tmp]))/1000
    clsL[[n_cls]]  <- as.vector((d_tmp1 <= r_tmp))
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

  if (n_cls1 > 0) {
    n_obs <- sapply(clsL, sum)[-(1:n_cls1)]
  } else {
    n_obs <- sapply(clsL, sum)
  }

  rownames(Clusters) <- NULL                   # Update: 2019/09/08
  Clusters <- cbind(Clusters, n_obs=n_obs)     # Update: 2019/09/08

  n_cls2 <- n_cls - n_cls1 - 1                 # Update: 2019/09/08
  if (n_cls2 > 0) {                            # Update: 2019/09/08
    if (n_cls1 > 0) {
      Clusters <- rbind(Cls1st$Clusters[1:n_cls1,], Clusters)
    }
    Clusters <- cbind(Clusters, stage=c(rep(1,n_cls1),rep(2,(n_cls-n_cls1))))
  } else {
    Clusters <- cbind(Cls1st$Clusters, stage=1)
  }

  Clusters <- rbind(Clusters[,-3], NULL)      # Update: 2019/09/08

  # Update: 2019/09/09
  Coef <- Est.Coeff.TS.ST(yList, XList, Cls1stIndicator=Cls1st$Indicator, Cls2ndIndicator=clsL)

  return(list(Clusters = Clusters, Coef = Coef, Indicator = clsL))
}


#'@title Fit.Model.Clusters.ST
#'@description Fit a simple linear regression model with detected clusters.
#'@param yList The input data (as a list of vectors).
#'@param XList The input data (as a list of matrices).
#'@param Cls1st Output from \code{Find.Clusters.TS.ST1()} function in the 1st Stage.
#'@param Cls2nd Output from \code{Find.Clusters.TS.ST2()} function in the 2nd Stage.
#'@export
#'
Fit.Model.Clusters.TS.ST <- function(yList, XList, Cls1st, Cls2nd) {
  WList <- XList
  clsL <- Cls2nd$Indicator
  n_cls1 <- length(Cls1st$Indicator) - 1       # number of clusters in the 1st Stage
  n_cls2 <- length(clsL) - n_cls1- 1              # number of clusters in the 2nd Stage
  n_cls <- n_cls1 + n_cls2
  p <- dim(XList[[1]])[2]; T <- length(XList)

  coeff_namesList <- list()
  modelList <- list()
  for (t in 1:T) {
    coeff_namesList[[t]] <- paste("b0_", 0:(p-1), "_t", t, sep="")   # names for the prameter estimates
    modelList[[t]] <- paste(" + b0_", 0:(p-1), "_t", t, sep="", collapse = "")  # model to fit
  }

  if (n_cls > 0) {
    for (j in 1:n_cls) {
      for (t in 1:T) {
        coeff_namesList[[t]][(p*j+1):(p*j+p)] <- paste("c", j, "_", 0:(p-1), "_t", t, sep="")
        modelList[[t]] <- paste(modelList[[t]], paste(" + c", j, "_", 0:(p-1), "_t", t, sep="", collapse = ""), sep="")
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
    W <- as.matrix(Matrix::bdiag(W,WList[[t]]))
  }
  W <- W[-1,-1]

  data_cls <- data.frame(y,W)
  colnames(data_cls) <- c("y", coeff_names)
  stats::lm(model, data = data_cls)
}


#'@title Est.Coeff.SI.ST
#'@description Estimate coefficients via two-stage detection.
#'@param yList The input data (as a list of vectors).
#'@param XList The input data (as a list of matrices).
#'@param Cls1stIndicator Indicator of clusters identified in stage 1.
#'@param Cls2ndIndicator Indicator of clusters identified in stage 2.
#'@return List of coefficients
Est.Coeff.TS.ST <- function(yList, XList, Cls1stIndicator, Cls2ndIndicator) {
  WList <- XList
  clsL <- Cls2ndIndicator
  n_cls1 <- length(Cls1stIndicator) - 1       # number of clusters in the 1st Stage
  n_cls2 <- length(clsL) - n_cls1               # number of clusters in the 2nd Stage # Update: 2019/09/09
  p <- dim(XList[[1]])[2]; T <- length(XList)

  b_tmp <- list()
  coef_tmp <- list()
  for (t in 1:T) {
    b_tmp[[t]] <- solve(t(XList[[t]])%*%XList[[t]])%*%t(XList[[t]])%*%yList[[t]]
    coef_tmp[[t]] <- c(b_tmp[[t]],rep(NA,length(b_tmp[[t]])))
  }

  XList_cls <- list()
  yList_tmp <- yList
  if (n_cls1 > 0) {
    for (j in 1:n_cls1) {
      for (t in 1:T) {
        XList_cls[[t]] <- XList[[t]]*(clsL[[j]])
        WList[[t]] <- cbind(XList[[t]],XList_cls[[t]])
        b_tmp[[t]] <- solve(t(WList[[t]])%*%WList[[t]])%*%t(WList[[t]])%*%yList_tmp[[t]]
        coef_tmp[[t]] <- cbind(coef_tmp[[t]], b_tmp[[t]])
        yList_tmp[[t]] <- yList_tmp[[t]] - XList_cls[[t]]%*%b_tmp[[t]][(p+1):(2*p)]
      }
    }
  }

  if (n_cls2 > 0) {
    for (j in 1:n_cls2) {
      for (t in 1:T) {
        WList[[t]] <- cbind(XList[[t]],1*(clsL[[j+n_cls1]]))
        b_tmp[[t]] <- solve(t(WList[[t]])%*%WList[[t]])%*%t(WList[[t]])%*%yList_tmp[[t]]
        coef_tmp[[t]] <- cbind(coef_tmp[[t]], c(b_tmp[[t]],0))
        yList_tmp[[t]] <- yList_tmp[[t]] - (clsL[[j+n_cls1]])*b_tmp[[t]][(p+1)]
      }
    }
  }

  n_cls <- n_cls1 + n_cls2
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


