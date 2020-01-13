source('lib.R')

#### 161228 and 170111 ####
data <- read.csv("edstars.csv")
dataX <- data[, 2:87]
dataY <- data[["y"]]
xRow <- nrow(dataX)
load("161228//bigResult.RData")
#load("170111//bigResult.RData")

# Plot ROC 
png("170111//auc.png", width=480, height=480)
forest <- randomForest(x=dataX, y=as.factor(dataY), ntree=100)
resForest <- predict(forest, x=dataX, type="prob")

x <- c()
y <- c()
for (thresh in seq(0, 1, by=0.02)) {
  treeRes <- (resForest[, 2] >= thresh)*1
  x <- c(x, sum(treeRes*dataY)/sum(dataY))
  y <- c(y, sum(treeRes*(1-dataY)/sum(1-dataY)))
}
plot(y, x, xlim=c(0,1), ylim=c(0,1), xlab="False Positive", ylab="True Positive", type="l", col="blue")


for (i in 1:100) {
  thisTree <- bigResult[[i]]$tree
  res <- mimic.tree.predict(thisTree, dataX)
  x <- c()
  y <- c()
  for (thresh in seq(0, 1, by=0.02)) {
    treeRes <- (res[, 2] >= thresh)*1
    x <- c(x, sum(treeRes*dataY)/sum(dataY))
    y <- c(y, sum(treeRes*(1-dataY)/sum(1-dataY)))
  }
  points(y, x, type="l", col="red")
}
legend("bottomright", legend=c("Forest", "Tree"), col=c("blue", "red"), lty=c(1,1))
dev.off()

# Plot Structures of the Tree
treeID <- function(node, depth) {
  if(depth==0) return(NULL)
  if(is.null(node$split)) 
    return(c(treeID(node$lnode, depth-1), 0, 0, treeID(node$rnode, depth-1)))
  else
    return(c(treeID(node$lnode, depth-1), node$split[1:2], treeID(node$rnode, depth-1)))
}

par(mfrow=c(1,4))
par(mar=c(5.1,4.1,2.1,2.1))
treels <- list()

for(i in 1:100) {
  treels[[i]] <- bigResult[[i]]$tree
}
  
for(depth in 1:4) {
  count <- c()
  board <- c()
  for(i in 1:100) {
    treei <- treeID(treels[[i]], depth)
    len <- length(treei)
    tl <- seq(from=1, to=len, by=2)
    flag <- FALSE
    if(length(count)>0)
      for(j in 1:length(count)) {     
        if(max(abs(treei[tl]-board[j,tl]))<0.1) {
          count[j] <- count[j]+1
          flag <- TRUE
          break
        }
      }
    if(!flag) {
      count <- c(count, 1)
      board <- rbind(board, treei)
    }
  }
  print(board)
  print(count)
  count <- sort(count, decreasing=TRUE)
  names(count) <- c(1:length(count))
  barplot(count, ylim=c(0,100), xlab=paste("Depth:", depth), horiz=FALSE)
}

#### 170128 and 170201 ####

print.tree <- function(node, indent.str=NULL, name=NULL) {
  if (node$leaf) {
    cat(paste(c(" return 0/1: ", 
                round(node$value[1], digits=3), "/", round(node$value[2], digits=3)), 
              collapse="", sep=""))
  } 
  else {
    var <- as.integer(node$split[1])
    if (is.null(name)) var <- paste("x", var, collapse="", sep="") else var <- name[var]
    cat(paste(c("\n", indent.str, "if ", var, " <= ", round(node$split[2],digits=1), ":"), collapse="", sep=""))
    print.tree(node$lnode, indent.str=paste(c(indent.str, " |  "), collapse="", sep=""), name=name)
    cat(paste(c("\n", indent.str, "if ", var, " >  ", round(node$split[2],digits=1), ":"), collapse="", sep=""))
    print.tree(node$rnode, indent.str=paste(c(indent.str, "    "), collapse="", sep=""), name=name)
  }
}

part1 <- 1:467
part2 <- 468:934
part3 <- 935:1401

data <- read.csv("edstars.csv")
dataX <- data[, 2:87]
dataY <- data[["y"]]
xRow <- nrow(dataX)
xnames <- names(dataX)

train1 <- data[c(part2, part3), ]
train1.X <- train1[, 2:87]
train1.Y <- train1[["y"]]
test1 <- data[part1, ]
test1.X <- test1[, 2:87]
test1.Y <- test1[["y"]]

train2 <- data[c(part1, part3), ]
train2.X <- train2[, 2:87]
train2.Y <- train2[["y"]]
test2 <- data[part2, ]
test2.X <- test2[, 2:87]
test2.Y <- test2[["y"]]

train3 <- data[c(part1, part2), ]
train3.X <- train3[, 2:87]
train3.Y <- train3[["y"]]
test3 <- data[part3, ]
test3.X <- test3[, 2:87]
test3.Y <- test3[["y"]]

# Build 3 RFs
load("170130//forests.RData")
forest1 <- forests[[1]]
forest1 <- randomForest(x=train1.X, y=as.factor(train1.Y), ntree=100)
forest2 <- forests[[2]]
forest2 <- randomForest(x=train2.X, y=as.factor(train2.Y), ntree=100)
forest3 <- forests[[3]]
forest3 <- randomForest(x=train3.X, y=as.factor(train3.Y), ntree=100)

load("170130//d4.p1.RData")
load("170130//d4.p2.RData")
load("170130//d4.p3.RData")
load("170130//d5.p1.RData")
load("170130//d5.p2.RData")
load("170130//d5.p3.RData")
load("170130//d6.p1.RData")
load("170130//d6.p2.RData")
load("170130//d6.p3.RData")

treeSets <- list(bigResult.d4.p1, bigResult.d4.p2, bigResult.d4.p3,
                 bigResult.d5.p1, bigResult.d5.p2, bigResult.d5.p3,
                 bigResult.d6.p1, bigResult.d6.p2, bigResult.d6.p3)


## Plot the ROC regarding forests
png("170130//3FCVI2Train.png", width=480, height=480)
resForest <- predict(forest2, newdata=train2.X, type="prob")

x <- c()
y <- c()
for (thresh in seq(0, 1, by=0.02)) {
  treeRes <- (resForest[, 2] >= thresh)*1
  x <- c(x, sum(treeRes*train2.Y)/sum(train2.Y))
  y <- c(y, sum(treeRes*(1-train2.Y)/sum(1-train2.Y)))
}
plot(y, x, xlim=c(0,1), ylim=c(0,1), 
     xlab="False Positive", ylab="True Positive", type="l", col="black", lwd=3,
     main="3-Fold CV, Iteration 2 (Train)")
abline(v=(0:10)/10, h=(0:10)/10, col="grey", lty=3)

res <- list()
colSets = rep(c("blue", "red", "darkorchid"), each=3)
for (j in c(2,5,8)) {
  bigResult <- treeSets[[j]]
  for (i in 1:100) {
    thisTree <- bigResult[[i]]$tree
    res[[i]] <- mimic.tree.predict(thisTree, train2.X)
  }
  x <- c()
  y <- c()
  for (thresh in seq(0, 1, by=0.01)) {
    wx <- 0
    wy <- 0
    for (i in 1:100) {
      thisRes <- res[[i]][,2]
      treeRes <- (thisRes >= thresh)*1
      wx <- wx + sum(treeRes*train2.Y)/sum(train2.Y)
      wy <- wy + sum(treeRes*(1-train2.Y))/sum(1-train2.Y)
    }
    x <- c(x, wx/100)
    y <- c(y, wy/100)
  }
  points(y, x, type="l", lty=1, col=colSets[j], lwd=2)
}
legend("bottomright", legend=c("Forest", "Tree Depth=4", "Tree Depth=5", "Tree Depth=6"), 
       col=c("black", "blue", "red", "darkorchid"), lty=1)
dev.off()

count <- c()
board <- c()
bigResult <- treeSets[[9]]
for(i in 1:100) {
  treei <- treeID(bigResult[[i]]$tree, 5)
  len <- length(treei)
  tl <- seq(from=1, to=len, by=2)
  flag <- FALSE
  if(length(count)>0)
    for(j in 1:length(count)) {     
      if(max(abs(treei[tl]-board[j,tl]))<0.1) {
        count[j] <- count[j]+1
        flag <- TRUE
        break
      }
    }
  if(!flag) {
    count <- c(count, 1)
    board <- rbind(board, treei)
  }
}
print(board)
print(count)

treeID(bigResult[[1]]$tree, depth=5) == board[1, ]
print.tree(bigResult[[1]]$tree, name=names(dataX))

#### 170201 new ####
