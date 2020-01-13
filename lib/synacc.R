source('treeAndForest.R')

timeStampStart <- proc.time()

synDataGen <- function(n)
{
  x1 <- runif(n)
  x2 <- runif(n)
  x3 <- runif(n)
  x4 <- runif(n)
  x5 <- runif(n)
  p <- (x1>0.5)*(3*(x2>0.7)+-3*(x2<0.7 & x2>0.2)-4*(x2<0.2)) + 
    ((x1<0.5) & (x5<0.5))*(3*((x3+x4)>1.4)+2*(((x3+x4)<1.4) & ((x3+x4)>0.5))-2*((x3+x4)<0.5)) +
    ((x1<0.5) & (x5>0.5))*2
  p <- exp(p)/(1+exp(p))
  x <- cbind(x1, x2, x3, x4, x5)
  colnames(x)=NULL
  y <- rbinom(n, 1, p)
  return(list(x=x, y=y, df=data.frame(x,y)))
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
    return(newData)
  }
  # Return the generator
  return(cov.generator.kernel)
}

result2 <- c()
result2class <- c()
result2cart <- c()
result2rf <- c()
result2apt <- c()
treels <- list()

for(TIME in 1:100)
{
  print(TIME)
  dataRaw <- synDataGen(1000)
  xRaw <- dataRaw$x
  yRaw <- as.factor(dataRaw$y)
  forest <- randomForest(x=xRaw, y=yRaw, ntree=100)
  
  cart <- rpart(y~., data=dataRaw$df, method="class", control=rpart.control(maxdepth=4, minsplit=5))
  tree <- mimic.forest(xRaw, yRaw, ccg=synCCG, ncand=10, alpha = 0.1, tree.depth=5, max.leaf=5, forest=forest,
                       maxcut=1000000, maxinc=1000000, mininc=10000)
  treels[[TIME]] <- tree
  result2 <- c(result2, compare.tree.and.forest(1000, synCCG(xRaw), tree, forest)$score)
  result2class <- c(result2class, compare.tree.and.forest.class(1000, synCCG(xRaw), tree, forest)$score)
  testData <- synDataGen(2000)
  result2cart <- c(result2cart, sum(predict(cart, newdata=testData$df, type="class")==factor(testData$y))/2000)
  result2rf <- c(result2rf, sum(predict(forest, newdata=testData$x, type="class")==factor(testData$y))/2000)
  result2apt <- c(result2apt, sum(mimic.tree.predict.factor(tree, testData$x)==factor(testData$y))/2000)
  #   forest2 <- randomForest(x=xRaw, y=yRaw, ntree=1000)
  #   result.ff <- c(result.ff, compare.forest.and.forest(1000, tempCovGen, forest, forest2)$score)
  
}

summary(result2cart)
summary(result2rf)
summary(result2apt)

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

sfSource('lib.R')

######## HERE BEGINS THE SIMULATION ########

synDataGen <- function(n)
{
  x1 <- runif(n)
  x2 <- runif(n)
  x3 <- runif(n)
  x4 <- runif(n)
  x5 <- runif(n)
  p <- (x1>0.5)*(3*(x2>0.7)+-3*(x2<0.7 & x2>0.2)-4*(x2<0.2)) + 
    ((x1<0.5) & (x5<0.5))*(3*((x3+x4)>1.4)+2*(((x3+x4)<1.4) & ((x3+x4)>0.5))-2*((x3+x4)<0.5)) +
    ((x1<0.5) & (x5>0.5))*2
  p <- exp(p)/(1+exp(p))
  x <- cbind(x1, x2, x3, x4, x5)
  colnames(x)=NULL
  y <- factor(rbinom(n, 1, p))
  return(list(x=x, y=y,df=data.frame(x=x, y=y)))
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
    newData <- newData + matrix(rnorm(nSample*nColX, mean=0, sd=0.01), nrow=nSample)
    return(newData)
  }
  # Return the generator
  return(cov.generator.kernel)
}

simulation <- function(dummy)
{
  dataRaw <- synDataGen(1000)
  xRaw <- dataRaw$x
  yRaw <- dataRaw$y
  print("learn forest")
  forest <- randomForest(x=xRaw, y=yRaw, data=dataRaw$df, ntree=100)
  
  print("learn tree")
  cart <- rpart(y~., data=dataRaw$df, method="class", control=rpart.control(maxdepth=4, minsplit=5))
  
  print("learn mimic tree")
  tree <- mimic.forest(xRaw, yRaw, ccg=synCCG, ncand=20, alpha = 0.1, tree.depth=5, max.leaf=5, forest=forest,
                       maxcut=10000, maxinc=10000, mininc=1000)
  print("comparison")
  result2 <- compare.tree.and.forest(1000, synCCG(xRaw), tree, forest)
  result2class <- compare.tree.and.forest.class(1000, synCCG(xRaw), tree, forest)
  print("new prediction")
  testData <- synDataGen(2000)
  result2cart <- sum(predict(cart, newdata=testData$df, type="class")==testData$y)/1000
  result2rf <- sum(predict(forest, newdata=testData$x, type="class")==testData$y)/1000
  result2apt <- sum(mimic.tree.predict.factor(tree, testData$x)==testData$y)/1000
  return(list(accCART=result2cart, accRF=result2rf, accSTA=result2apt, con=result2, conclass=result2class,
              tree=tree, rf=forest, cart=cart))
}

sfExportAll()

bigResult <- list()
dir.create("0416")
for(TIME in 1:5)
{
  bigResult <- c(bigResult, sfLapply(1:20, simulation))
  save(bigResult, file="0416//bigResult_syn_1.RData")
}

sfStop()

save.image("synacc0413_1.RData")

boxplot(cbind(CART= result2cart, RF=result2rf, STA=result2apt), xlab=expression('N'[ps]*'=100,000'), ylab="Predictive Accuracy")
boxplot(cbind(PROB= result2, CLASS=1-result2class), xlab=expression('N'[ps]*'=100,000'), ylab="Mimicking Accuracy")

par(mfrow=c(1,2))

par(mfrow=c(2,2))

for(depth in 1:4)
{
  count <- c()
  board <- c()
  for(i in 1:100)
  {
    treei <- treeID(treels[[i]], depth)
    flag <- FALSE
    if(length(count)>0)
      for(j in 1:length(count))
      {     
        if(max(abs(treei-board[j,]))<0.1) 
        {
          count[j] <- count[j]+1
          flag <- TRUE
          break
        }
      }
    if(!flag)
    {
      count <- c(count, 1)
      board <- rbind(board, treei)
    }
  }
  print(board)
  print(count)
}