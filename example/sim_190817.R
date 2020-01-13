source("mimicTreeF.R")

sigmoid <- function(x) return(1/(1+exp(-x)))

#### PART 1 : Synthetic Data 1 ####
set.seed(42)
x.generator <- function(n, d = 3) {
    return(matrix(runif(n*d, -1, 1), nrow=n))
}

y.generator <- function(x) {
    return(as.numeric(
        runif(nrow(x)) < sigmoid(-3*(x[, 1]<0) + 3*x[, 2] * (x[, 1] > 0))
    ))
}

ccg <- function(xData, path=NULL) {
    nRowX <- nrow(xData)
    nColX <- ncol(xData)
    cov.generator.kernel <- function(nSample) {
        newData <- xData[sample(nRowX, nSample, replace=TRUE), ]
        newData <- newData + matrix(rnorm(nSample*nColX, 0, 0.05), nrow=nSample)
        if(!is.null(path)) {
            for(j in 1:nrow(path))
                if(path[j,3]==0) 
                    newData[newData[,path[j,1]]>path[j,2], path[j,1]] <- path[j,2]
                else 
                    newData[newData[,path[j,1]]<path[j,2], path[j,1]] <- path[j,2]
        }
        return(newData)
    }
    return(cov.generator.kernel)
}

X <- x.generator(10000)
y <- as.factor(y.generator(X))

forest <- randomForest(X, y, sampsize = 1000, replace = FALSE, keep.inbag = TRUE, ntree = 10000)
mimic_tree <- mimic.forest(xData=X, yData=y, ccg=ccg, tree.depth=4, forest=forest, maxcut=500000, maxinc=300000, mininc=100000)
save(mimic_tree, file="PART1.RData")

mimic_tree_sameforest = list()
for (i in 1:100) {
  mimic_tree_sameforest[[i]] <- mimic.forest(xData=X, yData=y, ccg=ccg, tree.depth=4, forest=forest, maxcut=500000, maxinc=300000, mininc=100000)
}
save(mimic_tree_sameforest, file="PART1_1.RData")


#### PART 2 : Synthetic Data 2 ####
set.seed(42)
x.generator <- function(n, d = 5) {
    return(matrix(runif(n*d, -1, 1), nrow=n))
}

y.generator <- function(x) {
    return(as.numeric(
        runif(nrow(x)) < sigmoid(-3 * (x[,3] * x[, 4]) * (x[, 1]<0) + 3 * x[, 2] * (x[, 1] > 0))
    ))
}

ccg <- function(xData, path=NULL) {
    nRowX <- nrow(xData)
    nColX <- ncol(xData)
    cov.generator.kernel <- function(nSample) {
        newData <- xData[sample(nRowX, nSample, replace=TRUE), ]
        newData <- newData + matrix(rnorm(nSample*nColX, 0, 0.05), nrow=nSample)
        if(!is.null(path)) {
            for(j in 1:nrow(path))
                if(path[j,3]==0) 
                    newData[newData[,path[j,1]]>path[j,2], path[j,1]] <- path[j,2]
                else 
                    newData[newData[,path[j,1]]<path[j,2], path[j,1]] <- path[j,2]
        }
        return(newData)
    }
    return(cov.generator.kernel)
}

X <- x.generator(10000)
y <- as.factor(y.generator(X))

forest <- randomForest(X, y, sampsize = 1000, replace = FALSE, keep.inbag = TRUE, ntree = 10000)
mimic_tree <- mimic.forest(xData=X, yData=y, ccg=ccg, tree.depth=4, forest=forest, maxcut=500000, maxinc=300000, mininc=100000)
save(mimic_tree, file="PART2.RData")

#### PART 3 : Breast Cancer ####
set.seed(42)
data <- read.csv("breast-cancer-wisconsin.data", header=FALSE)
# 1. Sample code number: id number 
# 2. Clump Thickness: 1 - 10 
# 3. Uniformity of Cell Size: 1 - 10 
# 4. Uniformity of Cell Shape: 1 - 10 
# 5. Marginal Adhesion: 1 - 10 
# 6. Single Epithelial Cell Size: 1 - 10 
# 7. Bare Nuclei: 1 - 10 
# 8. Bland Chromatin: 1 - 10 
# 9. Normal Nucleoli: 1 - 10 
# 10. Mitoses: 1 - 10 
# 11. Class: (2 for benign, 4 for malignant)
X <- cbind(data$V2, data$V3, data$V4, data$V5, data$V6, data$V7, data$V8, data$V9, data$V10)
y <- as.factor(as.numeric(data$V11 == 2))

forest <- randomForest(X, y, sampsize = 200, replace = FALSE, keep.inbag = TRUE, ntree = 10000)
mimic_tree <- mimic.forest(xData=X, yData=y, ccg=ccg.f, tree.depth=4, forest=forest, maxcut=500000, maxinc=300000, mininc=100000)
save(mimic_tree, file="PART3.RData")

#### PART 4 : Wine ####
set.seed(42)
data <- read.table("winequality-white.csv", header=TRUE, sep=";")
predictors <- c("fixed.acidity", "volatile.acidity", "citric.acid",
                "residual.sugar", "chlorides", "free.sulfur.dioxide",
                "total.sulfur.dioxide", "density", "pH",
                "sulphates", "alcohol")
X <- data[predictors]
y <- as.factor(as.numeric(data[["quality"]] > 5))
ccg <- function(xData, path=NULL)
{
    nRowX <- nrow(xData)
    nColX <- ncol(xData)
    scale <- diag(c(0.0843868228, 0.0100794548, 0.0121019804,
                    0.5072057784, 0.0021847968, 1.7007137325,
                    4.2498064554, 0.0002990907, 0.0151000600,
                    0.0114125834, 0.1230620568))
    cov.generator.kernel <- function(nSample)
    {
        newData <- xData[sample(nRowX, nSample, replace=TRUE), ]
        newData <- newData + matrix(rnorm(nSample*nColX), nrow=nSample)%*%scale
        if(!is.null(path))
        {
            for(j in 1:nrow(path))
                if(path[j,3]==0) 
                    newData[newData[,path[j,1]]>path[j,2], path[j,1]] <- path[j,2]
                else 
                    newData[newData[,path[j,1]]<path[j,2], path[j,1]] <- path[j,2]
        }
        return(newData)
    }
    return(cov.generator.kernel)
}

forest <- randomForest(X, y, sampsize = 1000, replace = FALSE, keep.inbag = TRUE, ntree = 10000)
mimic_tree <- mimic.forest(xData=X, yData=y, ccg=ccg, tree.depth=4, forest=forest, maxcut=500000, maxinc=300000, mininc=100000)
save(mimic_tree, file="PART4.RData")


#### PART 5 : Draw Trees ####
treeID <- function(node, depth)
{
  if(depth==0) return(NULL)
  if(is.null(node$split)) 
    return(c(treeID(node$lnode, depth-1), 0, 0, treeID(node$rnode, depth-1)))
  else
    return(c(treeID(node$lnode, depth-1), node$split[1:2], treeID(node$rnode, depth-1)))
}
library(ggplot2)
sheet <- data.frame()

depth = 3
count <- c()
board <- c()
for(i in 1:100) {
  treei <- treeID(treels[[i]], depth)
  len <- length(treei)
  tl <- seq(from=1, to=len, by=1)
  flag <- FALSE
  if(length(count)>0)
    for(j in 1:length(count)) {     
      if(max(abs(treei[tl]-board[j,tl]))<0.05) {
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
for(i in 1:length(count)) {
  sheet <- rbind(sheet, data.frame(name=titles[TIME], depth=depth, class=i, count=count[i]))
}



tempPlot <- ggplot(sheet, aes(depth, y=count)) + theme(legend.position="none") + facet_wrap(~name)
tempPlot + geom_bar(stat="identity", position="stack", colour="white") + scale_fill_brewer() + xlab("") + ylab("Unique Structure Counts") + scale_x_continuous(breaks=c(3,4), labels=c("4-layer","5-layer")) 


