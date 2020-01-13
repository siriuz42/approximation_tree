source('treeAndForest.R')

timeStampStart <- proc.time()

synDataGen <- function(n)
{
  x1 <- runif(n)
  x2 <- runif(n)
  x3 <- runif(n)
  x4 <- runif(n)
  x5 <- runif(n)
  x6 <- runif(n)
  x7 <- runif(n)
  x8 <- runif(n)

  p <- x1+2*x2-3*x3+sin(x4*x5*6.28)+x6*x7+x7*x8
  p <- exp(p)/(1+exp(p))
  x <- cbind(x1, x2, x3, x4, x5, x6, x7, x8)
  y <- rbinom(n, 1, p)
  return(list(x=x, y=y))
}

synCCG <- function(xData, path=NULL)
{
  # Define the pertubation
  nRowX <- nrow(xData)
  nColX <- ncol(xData)
  # Generator samples
  cov.generator.kernel <- function(nSample)
  {
    newData <- xData[sample(nRowX, nSample, replace=TRUE), ]
    newData <- newData + matrix(rnorm(nSample*nColX, mean=0, sd=0.05), nrow=nSample)
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

result2 <- c()
result2class <- c()
result.ff <- c()
tempCovGen <- synCCG(xRaw)
result1 <- c()
result2cart <- c()
strID <- c()
result2rf <- c()
result2apt <- c()

for(TIME in 1:100)
{
  print(TIME)
  dataRaw <- synDataGen(1000)
  xRaw <- dataRaw$x
  yRaw <- as.factor(dataRaw$y)
  forest <- randomForest(x=xRaw, y=yRaw, ntree=100)
  
  cart <- rpart(yRaw ~ xRaw, method="class", control=rpart.control(maxdepth=6, minsplit=5))
  tree <- mimic.forest(xRaw, yRaw, ccg=synCCG, ncand=20, alpha = 0.1, tree.depth=7, max.leaf=5, forest=forest,
                       maxcut=100000, maxinc=100000, mininc=1000)
  print(tree)
  strID <- c(strID, tree$id)
  result1 <- rbind(result1, tree$splits)
  result2 <- c(result2, compare.tree.and.forest(1000, tempCovGen, tree, forest)$score)
  result2class <- c(result2class, compare.tree.and.forest.class(1000, tempCovGen, tree, forest)$score)
  testData <- synDataGen(1000)
  result2cart <- c(result2cart, sum(predict(cart, newData=testData$x, type="class")==factor(testData$y))/1000)
  result2rf <- c(result2rf, sum(predict(forest, newData=testData$x, type="class")==factor(testData$y))/1000)
  result2apt <- c(result2apt, sum(mimic.tree.predict.factor(tree, testData$x)==factor(testData$y))/1000)
  #   forest2 <- randomForest(x=xRaw, y=yRaw, ntree=1000)
  #   result.ff <- c(result.ff, compare.forest.and.forest(1000, tempCovGen, forest, forest2)$score)
  
}

print(result2)
hist(result2)
summary(result2)

print(result2class)
hist(result2class)
summary(result2class)

timeStampEnd <- proc.time()
print(timeStampEnd - timeStampStart)
# 
# result <- result1
# tempM <- cbind(t(matrix(c(result[,1]), nrow=15)), t(matrix(c(result[,2]), nrow=15)))
# tempU <- unique(tempM)
# ans <- c()
# for(i in 1:nrow(tempU))
# {
#   # print(tempU[i, ])
#   count <- 0
#   for(j in 1:nrow(tempM))
#     if(sum(tempM[j, ]!=tempU[i, ])==0) count <- count+1
#   ans <- c(ans, count)
# }
# print(tempU)
# print(ans)

print(cbind(sort(unique(strID))))
save.image("synacc0401.RData")
