source('mimicTreeF.R')

data <- read.table("mdd-final.csv", sep=',', head=FALSE)
yRaw <- as.factor(data[1:(nrow(data)-1), 2])
xRaw <- data[1:(nrow(data)-1), 3:ncol(data)]

result2 <- c()
forest <- randomForest(x=xRaw, y=yRaw, ntree=100)

split <- matrix(c(59, 40, 1.5, 1.5), nrow=2)
yRaw <- predict(forest, xRaw, type="Prob")
dataGen <- cov.generator(xRaw)

# svar <- c()
# for(nSample in seq(1000,20000, by=1000))
# {
#   print(nSample)
#   xx <- dataGen(nSample)
#   yy <- predict(forest, xx, type="Prob")
#   svar <- c(svar, calculate.variance2(xx, yy, split)*nrow(xx))
# }
# plot(svar)

chart1 <- c()
chart2 <- c()
for(nSample in seq(1000, 10000, by=1000))
{
  print(nSample)
  temp <- c()
  for(TIME in 1:100)
  {
    print(TIME)
    xTemp <- dataGen(nSample)
    yTemp <- predict(forest, xTemp, type="Prob")
    giniTemp <- get.Gini(xTemp, yTemp, split)
    print(giniTemp)
    temp <- c(temp, diff(giniTemp))
  }
  chart1 <- c(chart1, var(temp)*nSample)
  chart2 <- c(chart2, calculate.variance2(xTemp, yTemp, split))
}

plot(chart1, ylim=c(0, 0.15), col="red")
points(chart2*seq(1000, 10000, by=1000), col="blue")