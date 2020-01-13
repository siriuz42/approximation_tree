library('rpart')
library('randomForest')
source('ptestu.R')

### Generator more covariantes from a given set. NO pertubation.
cov.generator <- function(xData)
{
  # Define the pertubation
  nRowX <- nrow(xData)
  nColX <- ncol(xData)
  
  # Generator samples
  cov.generator.kernel <- function(nSample)
  {
    newData <- xData[sample(nRowX, nSample, replace=TRUE), ]
    return(newData)
  }
  # Return the generator
  return(cov.generator.kernel)
}

### Generator more covariantes from a given set. Add some pertubation.
cov.generator2 <- function(xData)
{
  # Define the pertubation
  nRowX <- nrow(xData)
  nColX <- ncol(xData)
  pert <- c(1, -1, 0, 0, 0)
  # Generator samples
  cov.generator.kernel <- function(nSample)
  {
    newData <- xData[sample(nRowX, nSample, replace=TRUE), ]
    newData <- newData + matrix(sample(pert, nSample*nColX, replace=TRUE), nrow=nSample)
    return(newData)
  }
  # Return the generator
  return(cov.generator.kernel)
}

### Generator more covariantes from a given set. 
ccg.f <- function(xData, path=NULL)
{
  # Define the pertubation
  nRowX <- nrow(xData)
  nColX <- ncol(xData)
  pert <- c(1, -1, 0, 0, 0, 0, 0, 0, 0, 0)
  # Generator samples
  cov.generator.kernel <- function(nSample)
  {
    newData <- xData[sample(nRowX, nSample, replace=TRUE), ]
    newData <- newData + matrix(sample(pert, nSample*nColX, replace=TRUE), nrow=nSample)
    if(!is.null(path))
    {
      for(j in 1:nrow(path))
        if(path[j,3]==0) 
          newData[newData[,path[j,1]]>floor(path[j,2]-0.01), path[j,1]] <- floor(path[j,2]-0.01)
      else 
        newData[newData[,path[j,1]]<ceiling(path[j,2]), path[j,1]] <- ceiling(path[j,2])
    }
    return(newData)
  }
  # Return the generator
  return(cov.generator.kernel)
}

###Compute the Gini coefficients for all
get.Gini <- function(X, Y, split)
{
  #   print(dim(X))
  #   print(dim(Y))
  nSplit <- nrow(split)
  nSample <- nrow(X)
  nClass <- ncol(Y)
  Gini <- c()
  for(i in 1:nSplit)
  {
    #     splitCov <- split[i, 1]
    #     splitVal <- split[i, 2]
    tempFlag <- X[, split[i, 1]] < split[i, 2]
    tempLeft <- rbind(Y[tempFlag,])
    tempRight <- rbind(Y[!tempFlag,])
    #     print(split[i,])
    #     print(tempLeft)
    #     print(tempRight)
    ws <- sum(tempFlag)
    if(ws==0) 
    {
      tempLeft <- rbind(rep(0, nClass))
      tempRight <- Y
    }
    else if(ws==nSample) 
    {
      tempRight <- rbind(rep(0, nClass))
      tempLeft <- Y
    }
    else
    {
      tempLeft <- rbind(Y[tempFlag,])
      tempRight <- rbind(Y[!tempFlag,])
    }
    #     print(tempLeft)
    #     print(tempRight)
    tempLeft <- apply(tempLeft, c(2), mean)
    tempRight <- apply(tempRight, c(2), mean)
    #     print(tempLeft)
    #     print(tempRight)
    Gini <- c(Gini, 1-ws/nSample*t(tempLeft)%*%tempLeft
              -(nSample-ws)/nSample*t(tempRight)%*%tempRight)
    
  }
  return(Gini)
}

### Providing the prior distribution of coviriates. Here goes uniform distribution.
uniform.unit <- function(split)
{
  if(is.vector(split)) return(abs(split[3]-split[2]))
  ans <- 1
  for(i in 1:dim(split)[1])
    ans <- ans * abs(split[i, 3]-split[i, 2])
  return(ans)
}

### Estimate the ratio of each class on both sides of one split
single.count <- function(X, Y, split)
{
  ans <- array(0, dim=c(2, dim(Y)[2]))
  #   print("single.count")
  #   print(split)
  for(i in 1:dim(X)[1])
    if(X[i, split[1]]<split[2]) ans[1, ] <- ans[1, ]+Y[i, ] else ans[2, ] <- ans[2, ]+Y[i, ]
  return(ans/dim(X)[1])
}

### Estimate the ratio of each class after a cross of two splits
cross.count <- function(X, Y, split)
{
  ans <- array(0, dim=c(4, dim(Y)[2]))
  #   print("cross.count")
  for(i in 1:dim(X)[1])
  {
    if(X[i, split[1,1]]<split[1,2])
    {
      if(X[i, split[2,1]]<split[2,2]) ans[1, ] <- ans[1, ]+Y[i, ] else ans[2, ] <- ans[2, ]+Y[i, ]
    }
    else
    {
      if(X[i, split[2,1]]<split[2,2]) ans[3, ] <- ans[3, ]+Y[i, ] else ans[4, ] <- ans[4, ]+Y[i, ]
    }
  }
  return(ans/dim(X)[1])
}

### Calculate the variance
calculate.variance <- function(nClass, priorProb, priorCrossProb, theta, thetaCross)
{
  ans <- 0
  for(p in 1:2)
    for(q in 1:2)
    {
      temp <- sum(theta[p, q, ]^3*(1-theta[p, q, ]))
      for(i in 1:(nClass-1))
        for(j in (i+1):nClass)
          temp <- temp - 2*(theta[p, q, i]*theta[p, q,j])^2
      ans <- ans + temp*priorProb[p, q]^2
    }
  for(p in 1:2)
    for(q in 1:2)
    {
      pq <- (p-1)*2+q
      temp <- sum(theta[1, p, ]*theta[2,q, ]*thetaCross[pq, ]*(1-thetaCross[pq, ]))
      for(i in 1:nClass)
        for(j in 1:nClass)
          if(i!=j) temp <- temp - 2*theta[1, p, i]*theta[2, q, j]*thetaCross[pq, i]*thetaCross[pq, j]
      ans <- ans - priorCrossProb[pq]*sqrt(priorProb[1, p]*priorProb[2, q])*temp
    }
  return(4*ans)
}

### Give confidence assertion about the best split among all
compare.split <- function(X, Y, split, densityFun=uniform.unit)
{
  nX <- dim(X)[1]
  Gini <- get.Gini(X, Y, split)
  #   print(Gini)
  orderGini <- order(Gini, decreasing=TRUE)
  #   print(orderGini)
  top <- orderGini[1]
  nClass <- dim(Y)[2]
  nSplit <- dim(split)[1]
  
  # Calculate the prior distribution (priorBlah)
  priorProb <- matrix(0, nrow=nSplit, ncol=2)
  priorCrossProb <- matrix(0, nrow=nSplit, ncol=4) #11 12 21 22
  priorProb[1, 1] <- densityFun(c(split[top,1], 0, split[top,2]))
  priorProb[1, 2] <- 1-priorProb[1, 1]
  
  theta <- array(0, dim=c(nSplit, 2, nClass))
  crossTheta <- array(0, dim=c(nSplit, 4, nClass))
  #   print(single.count(X, Y, split[top, ]))
  theta[1, , ] <- single.count(X, Y, split[top, ])
  alpha <- c()
  sigma <- c()
  for(i in 2:nSplit)
  {
    bot <- orderGini[i]
    delta <- Gini[top] - Gini[bot]
    
    priorProb[i, 1] <- densityFun(c(split[bot,1], 0, split[bot,2]))
    priorProb[i, 2] <- 1 - priorProb[i, 1]
    priorCrossProb[i, 1] <- densityFun(rbind(c(split[bot,1], 0, split[bot,2]),
                                             c(split[top,1], 0, split[top,2])))
    priorCrossProb[i, 2] <- densityFun(rbind(c(split[bot,1], 1, split[bot,2]),
                                             c(split[top,1], 0, split[top,2])))
    priorCrossProb[i, 3] <- densityFun(rbind(c(split[bot,1], 0, split[bot,2]),
                                             c(split[top,1], 1, split[top,2])))
    priorCrossProb[i, 4] <- 1 - priorCrossProb[i, 1] - priorCrossProb[i, 2] - priorCrossProb[i, 3]
    
    theta[i, ,] <- single.count(X, Y, split[bot, ])
    crossTheta[i, ,] <- cross.count(X, Y, rbind(split[top, ], split[bot, ]))
    
    sigma2Hat <- calculate.variance(nClass, priorProb[c(1,i), ], priorCrossProb[i, ], theta[c(1,i), , ], crossTheta[i, , ])
    print(c(i, sigma2Hat))
    sigma2Hat <- sigma2Hat/nX
    #     print(sigma2Hat)
    sigma <- c(sigma, sigma2Hat)
    alpha <- c(alpha, 1-pnorm(delta, mean=0, sd=sqrt(sigma2Hat)))
  }
  return(list(orderGini, alpha, sum(alpha), sigma))
}

### Previous one fails in the scenario if samples are densities
### Use full estimation here. No need to specify prior density
calculate.variance2 <- function(X, Y, split)
{
  nX <- nrow(X)
  nClass <- ncol(Y)
  ll <- (X[, split[1,1]]<split[1,2]) & (X[,split[2,1]]<split[2,2])  ##left left
  lr <- (X[, split[1,1]]<split[1,2]) & (X[,split[2,1]]>=split[2,2]) ##left right
  rl <- (X[, split[1,1]]>=split[1,2]) & (X[,split[2,1]]<split[2,2])
  rr <- (X[, split[1,1]]>=split[1,2]) & (X[,split[2,1]]>=split[2,2])
  
  d <- ll|lr
  if(sum(d)>1) Theta <- -apply(as.matrix(Y[d, ]), c(2), sum)/sum(d)
  else if(sum(d)==1) Theta <- -Y[d, ] 
  else Theta <- rep(0, nClass)
  
  d <- rl|rr
  if(sum(d)>1) Theta <- c(Theta, -apply(as.matrix(Y[d, ]), c(2), sum)/sum(d))
  else if(sum(d)==1) Theta <- c(Theta, -Y[d, ])
  else Theta <- c(Theta, rep(0, nClass))
  
  d <- ll|rl
  if(sum(d)>1) Theta <- c(Theta, apply(as.matrix(Y[d, ]), c(2), sum)/sum(d))
  else if(sum(d)==1) Theta <- c(Theta, Y[d, ])
  else Theta <- c(Theta, rep(0, nClass))
  
  d <- lr|rr
  if(sum(d)>1) Theta <- c(Theta, apply(as.matrix(Y[d, ]), c(2), sum)/sum(d))
  else if(sum(d)==1) Theta <- c(Theta, Y[d, ])
  else Theta <- c(Theta, rep(0, nClass))
  
  Theta <- Theta*2
  
  bigY <- matrix(0, nrow = nX, ncol=4*nClass)
  for(i in 1:nX)
  {
    if(ll[i] || lr[i]) bigY[i, 1:nClass] <- Y[i, ]
    if(rl[i] || rr[i]) bigY[i, nClass+1:nClass] <- Y[i, ]
    if(ll[i] || rl[i]) bigY[i, 2*nClass+1:nClass] <- Y[i, ]
    if(lr[i] || rr[i]) bigY[i, 3*nClass+1:nClass] <- Y[i, ]
  }

  Sigma <- cov(bigY)
 
  return(t(Theta)%*%Sigma%*%Theta/nX)
}

#### Maybe a much better version
calculate.variance.v0 <- function(X, Y, split)
{
  nX <- nrow(X)
  nClass <- ncol(Y)
  l1 <- (X[, split[1,1]]<split[1,2])  ##left left
  r1 <- (X[, split[1,1]]>=split[1,2]) ##left right
  l2 <- (X[,split[2,1]]<split[2,2])
  r2 <- (X[,split[2,1]]>=split[2,2])
  
  d <- l1
  if(sum(d)>1) Theta <- -apply(as.matrix(Y[d, ]), c(2), sum)/sum(d)
  else if(sum(d)==1) Theta <- -Y[d, ] 
  else Theta <- rep(0, nClass)
  
  d <- r1
  if(sum(d)>1) Theta <- c(Theta, -apply(as.matrix(Y[d, ]), c(2), sum)/sum(d))
  else if(sum(d)==1) Theta <- c(Theta, -Y[d, ])
  else Theta <- c(Theta, rep(0, nClass))
  
  d <- l2
  if(sum(d)>1) Theta <- c(Theta, apply(as.matrix(Y[d, ]), c(2), sum)/sum(d))
  else if(sum(d)==1) Theta <- c(Theta, Y[d, ])
  else Theta <- c(Theta, rep(0, nClass))
  
  d <- r2
  if(sum(d)>1) Theta <- c(Theta, apply(as.matrix(Y[d, ]), c(2), sum)/sum(d))
  else if(sum(d)==1) Theta <- c(Theta, Y[d, ])
  else Theta <- c(Theta, rep(0, nClass))
  
  Theta <- Theta*2
  
  bigY <- cbind(Y*l1, Y*r1, Y*l2, Y*r2)
  
  Sigma <- cov(bigY)
  
  return(t(Theta)%*%Sigma%*%Theta/nX)
}

compare.split2 <- function(X, Y, split, maxRow=10000)
{  
  nX <- nrow(X)
  #   print("Getting Ginis")
  Gini <- get.Gini(X, Y, split)
  #   print(Gini)
  orderGini <- order(Gini, decreasing=FALSE)
  #   print(orderGini)
  top <- orderGini[1]
  nClass <- dim(Y)[2]
  nSplit <- dim(split)[1]  
  
  alpha <- c()
  sigma <- c()
  for(i in 2:nSplit)
  {
    bot <- orderGini[i]
    delta <- Gini[top] - Gini[bot]
    #     print("Calculating Variance")
    if(nX>maxRow) 
    {
      tempL <- sample(1:nX, maxRow) 
      sigma2Hat <- calculate.variance.v0(X[tempL, ], Y[tempL, ], rbind(split[top,], split[bot,]))
    }
    else sigma2Hat <- calculate.variance.v0(X, Y, rbind(split[top,], split[bot,]))
    
    if(sigma2Hat<0) sigma2Hat <- 0
    sigma <- c(sigma, sigma2Hat)
    alpha <- c(alpha, pnorm(delta, mean=0, sd=sqrt(sigma2Hat)))
  }
  return(list(orderGini=orderGini, alpha=alpha, sumAlpha=sum(alpha), sigma=sigma, gini=Gini[orderGini]))
}

### choose one of the splits at the a given level
run.test.v0 <- function(covGen, forest, split, alpha=0.05, maxcut=10000, maxinc=1000, mininc=10)
{
  nSample <- 1000
  alpha50 <- alpha/50
  newX <- c()
  newY <- c()
  X <- covGen(nSample)
  Y <- predict(forest, X, type="Prob")[]
  #   Y <- model.matrix(~predict(forest, X,)-1)
  print("begin_cs")
  result <- compare.split2(X, Y, split)
  print("end_cs")
  while(nSample<maxcut && result$sumAlpha>alpha)
  {
    # Throw some splits 
    tempSplit <- split[result$orderGini[1], ]
    for(i in 2:nrow(split))
      if(result$alpha[i-1]>alpha50) tempSplit <- rbind(tempSplit, split[result$orderGini[i], ])
    split <- tempSplit
    alphaTop <- alpha/result$sumAlpha*result$alpha[1]
    thresh <-  (qnorm(alphaTop)/(qnorm(result$alpha[1])))^2
    inc <- trunc(nSample*(thresh-1))
    if(inc>maxinc) inc <- maxinc
    if(inc<mininc) inc <- mininc
    print(inc)
    newX <- covGen(inc)
    #     print("Predicting")
    newY <- predict(forest, newX, type="Prob")[]
    #     newY <- model.matrix(~predict(forest, newX)-1)
    X <- rbind(X, newX)
    Y <- rbind(Y, newY)
    #     print("Comparing Splits")
    cat("begin_cs_a\n")
    result <- compare.split2(X, Y, split)
    cat("end_cs_a\n")
    nSample <- nSample + inc
    #     print(nSample)
  }
  
  # If we need some mechanic to manually choose between several
  suc <- FALSE
  if(result$sumAlpha>alpha)
  {
    label <- 1
    tempSplit <- rbind(split[result$orderGini[1], ])
    for(i in 2:length(result$orderGini)) 
      if (result$alpha[i-1]>alpha)
      {
        label <- label + 1
        tempSplit <- rbind(tempSplit, split[result$orderGini[i], ])
      }
    #     print(tempSplit)
    best <- which(tempSplit[1:label, 1]==min(tempSplit[1:label, 1]))
    best <- which(tempSplit[1:label, 1]==min(tempSplit[1:label, 1])&tempSplit[1:label, 2]==min(tempSplit[best, 2]))
    best <- tempSplit[best, ]
    #     print(split[best,])
  }
  else
  {
    suc <- TRUE
    best <- split[result$orderGini[1], ]
  }
  return(list(best=best, success=suc, nSample=nSample, orderGini=result[[1]], alpha=result[[2]], sigma=result[[4]], Gini=result[[5]]))
}

### choose one of the splits at the a given level
run.test.v1 <- function(covGen, forest, split, alpha=0.05, maxcut=10000, maxinc=1000, mininc=10)
{
  nSample <- 1000
  alpha50 <- alpha/50
  newX <- c()
  newY <- c()
  X <- covGen(nSample)
  Y <- predict(forest, X, type="Prob")[]
  #   Y <- model.matrix(~predict(forest, X,)-1)
  result <- compare.split2(X, Y, split)
  while(nSample<maxcut && result$sumAlpha>alpha)
  {
    # Throw some splits 
    tempSplit <- rbind(split[result$orderGini[1], ])
    for(i in 2:nrow(split))
      if(result$alpha[i-1]>alpha50) tempSplit <- rbind(tempSplit, split[result$orderGini[i], ])
    split <- tempSplit
    alphaTop <- alpha/result$sumAlpha*result$alpha[1]
    thresh <-  (qnorm(alphaTop)/(qnorm(result$alpha[1])))^2
    inc <- trunc(nSample*(thresh-1))
    if(inc>maxinc) inc <- maxinc
    if(inc<mininc) inc <- mininc
    #     print(inc)
    newX <- covGen(inc)
    #     print("Predicting")
    newY <- predict(forest, newX, type="Prob")[]
    #     newY <- model.matrix(~predict(forest, newX)-1)
    X <- rbind(X, newX)
    Y <- rbind(Y, newY)
    #     print("Comparing Splits")
    result <- compare.split2(X, Y, split)
    nSample <- nSample + inc
    #     print(nSample)
  }
  
  # If we need some mechanic to manually choose between several
  suc <- FALSE
  if(result$sumAlpha>alpha)
  {
    label <- 1
    tempSplit <- t(as.matrix(split[result$orderGini[1], ])) 
    for(i in 2:length(result$orderGini)) 
      if (result$alpha[i-1]>alpha)
      {
        label <- label + 1
        tempSplit <- rbind(tempSplit, split[result$orderGini[i], ])
      }
    #     print(tempSplit)
    best <- which(tempSplit[1:label, 1]==min(tempSplit[1:label, 1]))
    best <- which(tempSplit[1:label, 1]==min(tempSplit[1:label, 1])&tempSplit[1:label, 2]==min(tempSplit[best, 2]))
    best <- tempSplit[best, ]
    #     print(split[best,])
  }
  else
  {
    suc <- TRUE
    best <- split[result$orderGini[1], ]
  }
  return(list(best=best, success=suc, nSample=nSample, orderGini=result[[1]], alpha=result[[2]], sigma=result[[4]], Gini=result[[5]]))
}

run.trail2 <- function(covGen, forest, split, alpha=0.05, maxcut=10000, maxinc=1000, mininc=10)
{
  nSample <- 1000
  newX <- c()
  newY <- c()
  X <- covGen(nSample)
  Y <- predict(forest, X, type="Prob")[]
  result <- compare.split2(X, Y, split)
  while(nSample<maxcut && result$sumAlpha>alpha)
  {
    # Throw some splits 
    tempSplit <- split[result$orderGini[1], ]
    for(i in 2:nrow(split))
      if(result$alpha[i-1]>1e-4) tempSplit <- rbind(tempSplit, split[result$orderGini[i], ])
    split <- tempSplit
    alphaTop <- alpha/result$sumAlpha*result$alpha[1]
    thresh <-  (qnorm(alphaTop)/(qnorm(result$alpha[1])))^2
    inc <- trunc(nSample*(thresh-1))
    if(inc>maxinc) inc <- maxinc
    if(inc<mininc) inc <- mininc
    #     print(inc)
    newX <- covGen(inc)
    #     print("Predicting")
    newY <- predict(forest, newX, type="Prob")[]
    X <- rbind(X, newX)
    Y <- rbind(Y, newY)
    cat("Comparing Splits\n")
    result <- compare.split2(X, Y, split)
    nSample <- nSample + inc
    #     print(nSample)
  }
  
  # If we need some mechanic to manually choose between several
  suc <- FALSE
  if(result$sumAlpha>alpha)
  {
    best <- split[result$orderGini[1], ]
  }
  else
  {
    suc <- TRUE
    best <- split[result$orderGini[1], ]
  }
  return(list(best=best, success=suc, nSample=nSample, orderGini=result[[1]], alpha=result[[2]], sigma=result[[4]], Gini=result[[5]]))
}

run.trail3 <- function(covGen, forest, split, alpha=0.05, maxcut=10000, maxinc=1000, mininc=10)
{
  nSample <- 1000
  alpha50 <- alpha/50
  newX <- c()
  newY <- c()
  X <- covGen(nSample)
  Y <- forest(X)
  #   Y <- model.matrix(~predict(forest, X,)-1)
  result <- compare.split2(X, Y, split)
  while(nSample<maxcut && result$sumAlpha>alpha)
  {
    # Throw some splits 
    tempSplit <- split[result$orderGini[1], ]
    for(i in 2:nrow(split))
      if(result$alpha[i-1]>alpha50) tempSplit <- rbind(tempSplit, split[result$orderGini[i], ])
    split <- tempSplit
    alphaTop <- alpha/result$sumAlpha*result$alpha[1]
    thresh <-  (qnorm(alphaTop)/(qnorm(result$alpha[1])))^2
    inc <- trunc(nSample*(thresh-1))
    if(inc>maxinc) inc <- maxinc
    if(inc<mininc) inc <- mininc
    #     print(inc)
    newX <- covGen(inc)
    #     print("Predicting")
    newY <- forest(newX)
    #     newY <- model.matrix(~predict(forest, newX)-1)
    X <- rbind(X, newX)
    Y <- rbind(Y, newY)
    #     print("Comparing Splits")
    result <- compare.split2(X, Y, split)
    nSample <- nSample + inc
    #     print(nSample)
  }
  
  # If we need some mechanic to manually choose between several
  suc <- FALSE
  if(result$sumAlpha>alpha)
  {
    label <- 1
    tempSplit <- t(as.matrix(split[result$orderGini[1], ])) 
    for(i in 2:length(result$orderGini)) 
      if (result$alpha[i-1]>alpha)
      {
        label <- label + 1
        tempSplit <- rbind(tempSplit, split[result$orderGini[i], ])
      }
    print(tempSplit)
    best <- which(tempSplit[1:label, 1]==min(tempSplit[1:label, 1]))
    best <- which(tempSplit[1:label, 1]==min(tempSplit[1:label, 1])&tempSplit[1:label, 2]==min(tempSplit[best, 2]))
    best <- tempSplit[best, ]
    #     print(split[best,])
  }
  else
  {
    suc <- TRUE
    best <- split[result$orderGini[1], ]
  }
  return(list(best=best, success=suc, nSample=nSample, orderGini=result[[1]], alpha=result[[2]], sigma=result[[4]], Gini=result[[5]]))
}

split.tree <- function(xData, yData, forest, ccg, ncand=5, alpha = 0.05, nodeNumber=1, tree.depth=4, max.leaf=10, maxcut=10000, maxinc=1000, mininc=10, char="", path=NULL)
{
#   print(char)
#   print(dim(xData))
#   print(dim(yData))
  node <- NULL
  # Find all potential candidate splits
  nRowData <- nrow(xData)
  nCovData <- ncol(xData)
  cat("Gini\n")
  bGini <- get.Gini(xData, yData, rbind(c(1,-1000)))
  cat("Gini_end\n")
  split <- NULL
  for(i in 1:nCovData)
  {
    temp <- sort(unique(xData[, i]))
    ltemp <- length(temp) - 1
    if(ltemp>0)
    {
      split$cov <- c(split$cov, rep(i, ltemp))
      split$val <- c(split$val, temp[1:ltemp]+diff(temp)/2)
    }  
  }
  split <- cbind(split$cov, split$val)
  # Get the Gini coefs and choose the ncand cadidates

  # Split the tree 
  if(nodeNumber<2^(tree.depth-1) && nrow(xData)>max.leaf && bGini > 1e-3 && !is.null(split))
  {
    gini <- get.Gini(xData, yData, split)
    if(ncand < nrow(split)) split <- split[order(gini)[1:ncand], ]
    # Find the best split, need to change the function accordingly
    bsplit <- run.test.v0(ccg(xData, path), forest, split, alpha=alpha, maxcut=maxcut, maxinc=maxinc, mininc=mininc)
    # Zhengze's testing function
    pvalues <- c()
    for (TTime in 1:10) {
      testPoints <- ccg(xData, path)(10)
      pvalues <- c(pvalues, ptestu(testPoints, forest))
    }
    node$pvalues <- pvalues
    nSample <- bsplit$nSample
    bsplit <- bsplit$best
    node$split <- c(bsplit, nodeNumber, nRowData, nSample, pvalues)
    node$leaf <- FALSE

    left <- xData[, bsplit[1]] < bsplit[2]
    node$id <- paste(bsplit[1], round(bsplit[2],2))
    node$lnode <- split.tree(rbind(xData[left, ]), rbind(yData[left, ]), forest, ccg=ccg,
                             ncand=ncand, alpha = alpha, nodeNumber = nodeNumber*2, tree.depth=tree.depth, max.leaf=max.leaf, 
                             maxcut=maxcut, maxinc=maxinc, mininc=mininc, 
                             char=paste(char, 'L', sep=""), 
                             path=rbind(path, c(bsplit[1], bsplit[2], 0)))
    node$id <- paste(node$id, 'L', node$lnode$id)
    node$rnode <- split.tree(rbind(xData[!left, ]), rbind(yData[!left, ]), forest, ccg=ccg,
                             ncand=ncand, alpha = alpha, nodeNumber = nodeNumber*2+1, tree.depth=tree.depth, max.leaf=max.leaf, 
                             maxcut=maxcut, maxinc=maxinc, mininc=mininc, 
                             char=paste(char, 'R', sep=""),
                             path=rbind(path, c(bsplit[1], bsplit[2], 1)))
    node$id <- paste(node$id, 'R', node$rnode$id)
    node$splits <- as.matrix(rbind(node$split, node$lnode$splits, node$rnode$splits))
  }
  else
  {
    node$leaf <- TRUE
    node$lnode <- NULL
    node$rnode <- NULL
    node$value <- apply(yData, c(2), mean)
#     print(yData)
#     print(node$value)
  }
  return(node)
}

mimic.forest <- function(xData, yData, ccg=ccg.f, ncand=5, tree.depth=4, max.leaf=10, alpha=0.05, forest=NULL, maxcut=10000, maxinc=1000, mininc=10)
{
  # Train a random forest 
  if(is.null(forest)) forest <- randomForest(x=xData, y=yData, ntree=100)
  
  # Get the root of the mimic tree
#   yData <- model.matrix(~yData-1)
  yData <- predict(forest, xData, type="prob")[]
  root <- split.tree(xData, yData, forest, ccg=ccg, ncand=ncand, alpha = alpha, tree.depth=tree.depth, max.leaf=max.leaf,
                     maxcut=maxcut, maxinc=maxinc, mininc=mininc)
  colnames(root$splits) <- c("_Covariate", "_Value", "_Node No", "_Sample No", "_Pseudo Sample No", paste("p", 1:10, sep=""))
#   print(root$splits)
  return(root)
}
  
mimic.tree.predict <- function(tree, xData)
{
  yData <- NULL
  if(ncol(xData)==1) xData <- t(as.matrix(xData))
  for(i in 1:nrow(xData))
  {
    node <- tree
#     print(node)
    while(!node$leaf)
    {
#       print("==========")
#       print(node$split)
      if(xData[i, node$split[1]]<node$split[2]) node <- node$lnode else node <- node$rnode
    }
#     print(yData)
    yData <- rbind(yData, node$value)
  }
  return(yData)
}

mimic.tree.predict.factor <- function(tree, xData)
{
  result <- mimic.tree.predict(tree, xData)
  return(factor(apply(result, c(1), which.max)-1))
}
