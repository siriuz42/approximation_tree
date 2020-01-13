
library(snowfall)
# library(snow)
pbsnodefile = Sys.getenv("PBS_NODEFILE")
machines <- scan(pbsnodefile, what="")
machines
nmach = length(machines)

# Initializing the nodes
sfInit(parallel=TRUE,type='SOCK',cpus=nmach,socketHosts=machines)


sfSource('lib.R')

data <- read.table("mdd-final.csv", sep=',', head=FALSE)
yRaw <- data[1:(nrow(data)-1), 2]
for(i in 1:length(yRaw)) if(yRaw[i]>0) yRaw[i] <- 1
yRaw <- as.factor(yRaw)
xRaw <- data[1:(nrow(data)-1), 3:90]

forest <- randomForest(x=xRaw, y=yRaw, ntree=100)


evaluate <- function(TIME)
{
  message(TIME)
  tree <- mimic.forest(xRaw, yRaw, ncand=20, alpha = 0.10, tree.depth=5, forest=forest,
                       maxcut=1000000, maxinc=1000000, mininc=100000)
  message(DONE)
  result2 <- compare.tree.and.forest(10000, ccg.f(xRaw), tree, forest)$score
  result2class <- compare.tree.and.forest.class(10000, ccg.f(xRaw), tree, forest)$score
  return(list(tree=tree, res=result2, resclass=result2class))
}

sfExportAll()
bigResult <- sfLapply(c(1:10), evaluate)
save(bigResult, file="bigResult0414.RData")

# save.image("mddTree.RData")
# 
# save.image("mddTree0410_2.RData")
sfStop()