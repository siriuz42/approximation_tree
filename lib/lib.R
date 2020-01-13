library('rpart')
library('randomForest')

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

### Generator more covariantes from a given set. NO pertubation.
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
    tempFlag <- (X[, split[i, 1]] < split[i, 2])
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
  
  return(2*t(Theta)%*%Sigma%*%Theta/nX)
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
run.test.v0 <- function(covGen, forest, split, alpha=0.05, maxcut=10000, maxinc=1000, mininc=10,
                        passX=NULL, passY=NULL)
{
  alpha50 <- alpha/50
  if(is.null(passX)) nSample <- 1000 else nSample <- max(1000-nrow(passX), 0)
  X <- passX
  Y <- passY
  if(nSample>0)
  {
    newX <- covGen(nSample)
    newY <- predict(forest, newX, type="Prob")[]
    X <- rbind(X, newX)
    Y <- rbind(Y, newY)
  }
  nSample <- nrow(X)
  #   Y <- model.matrix(~predict(forest, X,)-1)
#   print("begin_cs")
  result <- compare.split2(X, Y, split)
#   print("end_cs")
  while(nSample<maxcut && result$sumAlpha>alpha)
  {
    # Throw some splits 
    tempSplit <- rbind(split[result$orderGini[1], ])
    for(i in 2:nrow(split))
      if(result$alpha[i-1]>alpha50) tempSplit <- rbind(tempSplit, split[result$orderGini[i], ])
    split <- tempSplit  
    if(result$alpha[1]<1e-5) inc <- mininc 
    else
    {
      alphaTop <- alpha/result$sumAlpha*result$alpha[1]
      thresh <-  (qnorm(alphaTop)/(qnorm(result$alpha[1])))^2
      inc <- trunc(nSample*(thresh-1))
#     cat("result$sumAlpha: ", result$sumAlpha, "\n")
#     cat("result$alpha[1]: ", result$alpha[1], "\n")
#     cat("alphaTop: ", alphaTop, "\n")
#     cat("qnorm(alphaTop): ", qnorm(alphaTop), "\n")
#     cat("qnorm(result$alpha[1]): ", qnorm(result$alpha[1]), "\n")
#     cat("thresh: ",  thresh, "\n")
#     cat("inc: ", inc, "\n")
        
      if(inc>maxinc) inc <- maxinc
      if(inc<mininc) inc <- mininc
    }
#     print(inc)
#     print("gen_b")
    newX <- covGen(inc)
    #     print("Predicting")
#     print("gen_m")
    newY <- predict(forest, newX, type="Prob")[]
#     print("gen_e")
    #     newY <- model.matrix(~predict(forest, newX)-1)
    X <- rbind(X, newX)
    Y <- rbind(Y, newY)
    rm(newX)
    rm(newY)
    gc()
    #     print("Comparing Splits")
#     print("begin_cs_a")
    result <- compare.split2(X, Y, split)
#     print("end_cs_a")
    nSample <- nSample + inc
    #     print(nSample)
  }
  rm(X)
  rm(Y)
  gc()
#   print("finish_b")
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
#   print("finish_s")

  return(list(best=best, success=suc, nSample=nSample))
#               usedX=X, usedY=Y))
}


split.tree.cart <-  function(xData, yData, 
                             nodeNumber=1, tree.depth=4, max.leaf=10, 
                             char="", path=NULL,
                             passX=NULL, passY=NULL)
{
  print(nodeNumber)
  node <- NULL
  # Find all potential candidate splits
  nRowData <- nrow(xData)
  nCovData <- ncol(xData)
  bGini <- get.Gini(xData, yData, rbind(c(1,-1000)))
  split <- NULL
  for(i in 1:nCovData)
  {
    temp <- sort(unique(xData[, i]))
    if (length(temp)>10) temp <- temp[seq(1, length(temp), floor(length(temp)/12))]
    ltemp <- length(temp) - 1
    if(ltemp>0) split <- rbind(split, cbind(rep(i, ltemp), temp[1:ltemp]+diff(temp)/2)) 
  }
  # Get the Gini coefs and choose the ncand cadidates
  
  # Split the tree 
  if(nodeNumber<2^(tree.depth-1) && nrow(xData)>max.leaf && bGini > 1e-3 && !is.null(split))
  {
    if(nrow(split)>1)
    {
      gini <- get.Gini(xData, yData, split)
      bsplit <- list(nSample=nrow(xData), best=split[order(gini)[1], ])
    }
    else bsplit <- list(nSample=nrow(xData), best=split[1,])
    
    nSample <- bsplit$nSample
    bsplit <- bsplit$best
    node$split <- c(bsplit, nodeNumber, nRowData, nSample)
    
    node$leaf <- FALSE
    
    left <- xData[, bsplit[1]] < bsplit[2]
    node$id <- paste(bsplit[1], round(bsplit[2],2))
    node$lnode <- split.tree.cart(rbind(xData[left, ]), rbind(yData[left, ]), 
                                  nodeNumber = nodeNumber*2, tree.depth=tree.depth, max.leaf=max.leaf, 
                                  char=paste(char, 'L', sep=""), 
                                  path=rbind(path, c(bsplit[1], bsplit[2], 0)))
    #                              passX=usedX[left, ], passY=usedY[left, ])
    node$id <- paste(node$id, 'L', node$lnode$id)
    node$rnode <- split.tree.cart(rbind(xData[!left, ]), rbind(yData[!left, ]), 
                                  nodeNumber = nodeNumber*2+1, tree.depth=tree.depth, max.leaf=max.leaf, 
                                  char=paste(char, 'R', sep=""),
                                  path=rbind(path, c(bsplit[1], bsplit[2], 1)))
    #                              passX=usedX[!left, ], passY=usedY[!left, ])
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
  
split.tree <- function(xData, yData, forest, ccg, ncand=5, alpha = 0.05, 
                       nodeNumber=1, tree.depth=4, max.leaf=10, 
                       maxcut=10000, maxinc=1000, mininc=10, char="", path=NULL,
                       passX=NULL, passY=NULL)
{
  print(nodeNumber)
#   print(char)
#   print(dim(xData))
#   print(dim(yData))
  node <- NULL
  # Find all potential candidate splits
  nRowData <- nrow(xData)
  nCovData <- ncol(xData)
#   print("Gini")
  bGini <- get.Gini(xData, yData, rbind(c(1,-1000)))
#  print("Gini_end")
  split <- NULL
  for(i in 1:nCovData)
  {
    temp <- sort(unique(xData[, i]))
    ltemp <- length(temp) - 1
    if(ltemp>0) split <- rbind(split, cbind(rep(i, ltemp), temp[1:ltemp]+diff(temp)/2)) 
  }
  # Get the Gini coefs and choose the ncand cadidates

  # Split the tree 
  if(nodeNumber<2^(tree.depth-1) && nrow(xData)>max.leaf && bGini > 1e-3 && !is.null(split))
  {
    if(nrow(split)>1)
    {
      gini <- get.Gini(xData, yData, split)
      if(ncand < nrow(split)) split <- split[order(gini)[1:ncand], ]
    # Find the best split, need to change the function accordingly
#     print("bsplit_b")
      bsplit <- run.test.v0(ccg(xData, path), 
                            forest, split, alpha=alpha, maxcut=maxcut, maxinc=maxinc, mininc=mininc,
                            passX=passX, passY=passY)
#     usedX <- bsplit$usedX
#     usedY <- bsplit$usedY
    }
    else bsplit <- list(nSample=0, best=split[1,])

    nSample <- bsplit$nSample
    bsplit <- bsplit$best
    node$split <- c(bsplit, nodeNumber, nRowData, nSample)
    
    node$leaf <- FALSE

    left <- xData[, bsplit[1]] < bsplit[2]
    node$id <- paste(bsplit[1], round(bsplit[2],2))
    node$lnode <- split.tree(rbind(xData[left, ]), rbind(yData[left, ]), forest, ccg=ccg,
                             ncand=ncand, alpha = alpha, nodeNumber = nodeNumber*2, tree.depth=tree.depth, max.leaf=max.leaf, 
                             maxcut=maxcut, maxinc=maxinc, mininc=mininc, 
                             char=paste(char, 'L', sep=""), 
                             path=rbind(path, c(bsplit[1], bsplit[2], 0)))
#                              passX=usedX[left, ], passY=usedY[left, ])
    node$id <- paste(node$id, 'L', node$lnode$id)
    node$rnode <- split.tree(rbind(xData[!left, ]), rbind(yData[!left, ]), forest, ccg=ccg,
                             ncand=ncand, alpha = alpha, nodeNumber = nodeNumber*2+1, tree.depth=tree.depth, max.leaf=max.leaf, 
                             maxcut=maxcut, maxinc=maxinc, mininc=mininc, 
                             char=paste(char, 'R', sep=""),
                             path=rbind(path, c(bsplit[1], bsplit[2], 1)))
#                              passX=usedX[!left, ], passY=usedY[!left, ])
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
  colnames(root$splits) <- c("_Covariate", "_Value", "_Node No", "_Sample No", "_Pseudo Sample No")
#   print(root$splits)
  return(root)
}

mimic.forest.cart <- function(xData, yData, xDataBuild, tree.depth=4, max.leaf=10, alpha=0.05, forest=NULL, maxcut=10000, maxinc=1000, mininc=10)
{
  # Train a random forest 
  if(is.null(forest)) forest <- randomForest(x=xData, y=yData, ntree=100)
  
  # Get the root of the mimic tree
  #   yData <- model.matrix(~yData-1)
  yDataBuild <- predict(forest, xDataBuild, type="prob")[]
  root <- split.tree.cart(xDataBuild, yDataBuild, tree.depth=tree.depth, max.leaf=max.leaf)
  colnames(root$splits) <- c("_Covariate", "_Value", "_Node No", "_Sample No", "_Pseudo Sample No")
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

compare.tree.and.forest <- function(n, covGenerator, tree, forest, eval=NULL)
{
  testData <- covGenerator(n)
  treeResult <- mimic.tree.predict(tree, testData)
  forestResult <- predict(forest, testData, type="prob")[]
  #   print(treeResult)
  #   print(forestResult)
  diff = 0.0
  for(i in 1:n)
    if(is.null(eval)) diff <- diff + max(abs(treeResult[i, ]-forestResult[i, ])) else diff <- diff + eval(treeResult[i,], forestResult[i])
  return(list(score = diff/n, treeResult = treeResult, forestResult = forestResult))
}

compare.tree.and.forest.class <- function(n, covGenerator, tree, forest, eval=NULL)
{
  testData <- covGenerator(n)
  treeResult <- mimic.tree.predict(tree, testData)
  forestResult <- predict(forest, testData, type="prob")[]
  treeClass <- apply(treeResult, c(1), which.max)
  forestClass <- apply(forestResult, c(1), which.max)
  #   print(treeResult)
  #   print(forestResult)
  return(list(score = sum(treeClass==forestClass)/n, treeResult = treeResult, forestResult = forestResult))
}

compare.tree.and.tree <- function(n, covGenerator, tree, forest, eval=NULL)
{
  testData <- covGenerator(n)
  treeResult <- mimic.tree.predict(tree, testData)
  forestResult <- mimic.tree.predict(forest, testData)
  
  diff = 0.0
  for(i in 1:n)
    if(is.null(eval)) diff <- diff + max(abs(treeResult[i, ]-forestResult[i, ])) else diff <- diff + eval(treeResult[i,], forestResult[i])
  return(list(score = diff/n, treeResult = treeResult, forestResult = forestResult))
}

compare.forest.and.forest <- function(n, covGenerator, tree, forest, eval=NULL)
{
  testData <- covGenerator(n)
  treeResult <- predict(tree, testData, type="prob")
  forestResult <- predict(forest, testData, type="prob")
  diff = 0.0
  for(i in 1:n)
    if(is.null(eval)) diff <- diff + max(abs(treeResult[i, ]-forestResult[i, ])) else diff <- diff + eval(treeResult[i,], forestResult[i])
  return(list(score = diff/n, treeResult = treeResult, forestResult = forestResult))
}

lumbermill <- function(forest)
{
  ntree <- forest$ntree
  warehouse <- array(list(), ntree)
  ncat <- 0
  for(i in 1:ntree)
  {
    tempTree <- getTree(forest, i)
    ncat <- max(c(ncat, tempTree[, 6]))
  }
  #   print(ncat)
  for(i in 1:ntree)
  {
    tempTree <- getTree(forest, i)
    canvas <- tempTree[,1:5]
    for(j in 1:ncat)
      canvas <- cbind(canvas, tempTree[,6]==j)
    #     canvas <- array(list(), nrow(tempTree))
    # #     print("hello0")
    #     for(j in nrow(tempTree):1)
    #       if(tempTree[j, 5] == -1)  #leaf
    #       {
    #         canvas[[j]]$leaf <- TRUE
    #         canvas[[j]]$lnode <- NULL
    #         canvas[[j]]$rnode <- NULL
    #         canvas[[j]]$value <- rep(0, ncat)
    #         canvas[[j]]$value[tempTree[j,6]] <- 1
    #       }
    #       else
    #       {
    # #         print("x")
    #         canvas[[j]]$leaf <- FALSE
    #         canvas[[j]]$split <- c(tempTree[j, 3], tempTree[j, 4])
    #         canvas[[j]]$lnode <- canvas[[tempTree[j,1]]]
    #         canvas[[j]]$rnode <- canvas[[tempTree[j,2]]]
    # #         print("xx")
    #       }
    # #     print("hello1")
    warehouse[[i]] <- canvas
  }
  return(warehouse)
}

forest.split <- function(forest)
{
  split <- NULL
  for(i in 1:forest$ntree)
  {
    tmp <- getTree(forest, i)
    split <- rbind(split, tmp[tmp[,5]==1, 3:4])
  }
  return(unique(split))
}


subforest <- function(forest, nrep)
{
  trees <- lumbermill(forest)
  ntree <- forest$ntree
  nCol <- ncol(trees[[1]])
  subforest.kernel <- function(xData)
  {
    yData <- c()
    for(i in 1:nrow(xData))
    {
      temp <- 0
      index <- sample(ntree, nrep)
      for(j in index)
      {
        node <- 1
        while(trees[[j]][node, 5]==1)
        {
          if(xData[i, trees[[j]][node, 3]]<trees[[j]][node, 4]) node <- trees[[j]][node, 1] else node <- trees[[j]][node, 2] 
        }
        #     print(yData)
        temp <- temp + trees[[j]][node, 6:nCol] 
      }
      yData <- rbind(yData, temp)
    }
    return(yData/nrep)
  }
  return(subforest.kernel)
}

treeID <- function(node, depth)
{
  if(depth==0) return(NULL)
  if(is.null(node$split)) 
    return(c(treeID(node$lnode, depth-1), 0, 0, treeID(node$rnode, depth-1)))
  else
    return(c(treeID(node$lnode, depth-1), node$split[1:2], treeID(node$rnode, depth-1)))
}
