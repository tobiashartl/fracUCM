# out of sample forecast experiment


#Packages
gc()
rm(list = ls())
library(fUCpack)
library(dplyr)
library(parallel)
library(CFFpack)
library(xtable)


# Wd, etc
setwd("/Users/tobias/Dokumente/Projekte/filtering unknown persistence/R/code")
source("./help functions/fUC_arma_approx.R")

# load data
load(file = "./app2/mRVOL.RData")

y <- data.m$lRVOL.close_m
y0 <- y[1]
y <- y-y0
plot(y, type="l")
corr <- F
START = 2
diffuse = T
det <- "none"
p   <- 1
h   <- 120


# load the forecasts
y.pred <- x.pred <- matrix(NA,ncol = h, nrow =length(y)+h)
y.arima.111 <- y.arima.213 <- y.arima.101 <- y.arima.201 <- 
    y.ari.210 <- y.ari.410 <- y.ari.200 <- y.ari.300 <-
    y.arfima.0d0 <- y.arfima.6d6 <- y.arfima.0d10 <- y.arfima.4d13 <- 
    y.RW <- matrix(NA,ncol = h, nrow =length(y)+h)
for(j in 372:length(y)){
    pred <- readRDS(paste("/Users/tobias/Dokumente/Projekte/Filtering unknown persistence/R/code/app2/iterative/par/pred_t",
                       j, ".RDS", sep=""))
    # store the predictions for y and for x
    diag(y.pred[(j+1):(j+h),1:h]) <- pred$y.hat[(j+1):(j+h)]
    diag(x.pred[(j+1):(j+h),1:h]) <- pred$x.hat[(j+1):(j+h)]
    
    
    
    # benchmarks
    pred.ben <- readRDS(paste("/Users/tobias/Dokumente/Projekte/Filtering unknown persistence/R/code/app2/iterative/ben/pred_t",
                          j, ".RDS", sep=""))
    diag(y.arima.111[(j+1):(j+h),1:h]) <- pred.ben$pred.arima.111[(j+1):(j+h)]
    diag(y.arima.213[(j+1):(j+h),1:h]) <- pred.ben$pred.arima.213[(j+1):(j+h)]
    diag(y.arima.101[(j+1):(j+h),1:h]) <- pred.ben$pred.arima.101[(j+1):(j+h)]
    diag(y.arima.201[(j+1):(j+h),1:h]) <- pred.ben$pred.arima.201[(j+1):(j+h)]
    
    diag(y.ari.210[(j+1):(j+h),1:h]) <- pred.ben$pred.ari.210[(j+1):(j+h)]
    diag(y.ari.410[(j+1):(j+h),1:h]) <- pred.ben$pred.ari.410[(j+1):(j+h)]
    diag(y.ari.200[(j+1):(j+h),1:h]) <- pred.ben$pred.ari.200[(j+1):(j+h)]
    diag(y.ari.300[(j+1):(j+h),1:h]) <- pred.ben$pred.ari.300[(j+1):(j+h)]
    
    diag(y.arfima.0d0[(j+1):(j+h),1:h]) <- pred.ben$pred.arfima.0d0[(j+1):(j+h)]
    diag(y.arfima.6d6[(j+1):(j+h),1:h]) <- pred.ben$pred.arfima.6d6[(j+1):(j+h)]
    diag(y.arfima.0d10[(j+1):(j+h),1:h]) <- pred.ben$pred.arfima.0d10[(j+1):(j+h)]
    diag(y.arfima.4d13[(j+1):(j+h),1:h]) <- pred.ben$pred.arfima.4d13[(j+1):(j+h)]
    
    
    
}
y.RW <- embed0(y, 121)[,-1] 
y.RW[1:372, ] <- NA
y.RW <- rbind(y.RW, matrix(NA, nrow = 120, ncol = 120))


###Benchmarks: AR, random walk, arfima
plot(y, type="l")
lines(y.pred[,1], col=5)
lines(x.pred[,1], col=2)


lines(y.pred[,12], col=2)
lines(y.pred[,24], col=3)
lines(y.pred[,60], col=4)
dim(y.pred)
# note: y.pred is now a data frame consisting of 120 forecasts (1 month up to 10 years=


plot(y, type="l")
hh <- 48
lines(y.pred[,hh], col=2)
lines(y.arima.111[,hh], col=3)
lines(y.arima.213[,hh], col=4)
lines(y.arima.101[,hh], col=5)
lines(y.arima.201[,hh], col=6)

plot(y, type="l")
lines(y.pred[,hh], col=2)
lines(y.ari.210[,hh], col=3)
lines(y.ari.410[,hh], col=4)
lines(y.ari.200[,hh], col=5)
lines(y.ari.300[,hh], col=6)

plot(y, type="l")
lines(y.pred[,hh], col=2)
lines(y.arfima.0d0[,hh], col=3)
lines(y.arfima.6d6[,hh], col=4)
lines(y.arfima.0d10[,hh], col=5)
lines(y.arfima.4d13[,hh], col=6)

hh=1
plot(y-y,ylim = c(150, -150), type="l")
lines(y-y.pred[,hh], col=2)
lines(y.RW[, hh], col=3)
lines(y.arfima.0d0[, hh], col=4)


calcMSE <- function(pred, y){
    sqrt(mean(((pred) - (y))^2))
}
calcASE <- function(pred, y){
    (mean(abs((pred) - (y))))
}
calcLp <- function(pred, y, p=3){
    (sum((abs(pred - y))^p))^(1/p)
}
calcQlike <- function(pred, y){
    pred <- (pred )/100
    y    <- (y ) / 100
    mean( exp(y)/exp(pred) - (y - pred) -1 , na.rm = TRUE)
}


# all models
h <- c(1, 3, 6, 12, 24, 48, 60, 90, 120)
RESULTSmse <- RESULTSabs <- RESULTSpmse <- RESULTSpabs <- matrix(NA, nrow = length(h), ncol = 8)
k=0
for(hh in h){
    k <- k+1
    data.check <- data.frame(
        y.pred = y.pred[, hh], 
        #y.LR   = x.pred[, hh],
        y.arima.111 = y.arima.111[, hh], 
        y.arima.213 = y.arima.213[, hh], 
        #y.arima.101 = y.arima.101[, hh], 
        #y.arima.201 = y.arima.201[, hh], 
        y.ari.210 = y.ari.210[,hh],
        y.ari.410 = y.ari.410[,hh],
        #y.ari.200 = y.ari.200[,hh],
        #y.ari.300 = y.ari.300[,hh],
        #y.arfima.0d0 = y.arfima.0d0[,hh],
        #y.arfima.6d6 = y.arfima.6d6[,hh],
        y.arfima.0d10 = y.arfima.0d10[,hh],
        y.arfima.4d13 = y.arfima.4d13[,hh],
        y.RW = y.RW[, hh]
    )
    
    results <- apply(data.check, 2, function(x) calcMSE(x[(372 + hh):length(y)], y[(372 + hh):length(y)]))
    results <- results/results[1]
                     
    cat("Horizon ", hh, "\n")
    print(results)
    RESULTSmse[k, ] <- results
    results2 <- apply(data.check, 2, function(x) calcLp(x[(372 + hh):length(y)], y[(372 + hh):length(y)], p=4))
    results2 <- results2/results2[1]
    
    RESULTSabs[k, ] <- results2
    
    
    
    
    dm <- apply(data.check[,-1], 2, function(x) forecast::dm.test(data.check[(372 + hh):length(y),1] - y[(372 + hh):length(y)],
                                                             x[(372 + hh):length(y)] - y[(372 + hh):length(y)],
                                                             alternative = "two.sided", 
                                                             h=hh,
                                                             power = 2)$p.value)
    dm2 <- apply(data.check[,-1], 2, function(x) forecast::dm.test(data.check[(372 + hh):length(y),1] - y[(372 + hh):length(y)],
                                                                  x[(372 + hh):length(y)] - y[(372 + hh):length(y)],
                                                                  alternative = "two.sided", 
                                                                  h=hh,
                                                                  power = 4)$p.value)
    RESULTSpmse[k, ] <- c(NA, dm)
    RESULTSpabs[k, ] <- c(NA, dm2)
    
}

colnames(RESULTSmse) <- colnames(RESULTSabs)  <- c("fUC", "ARIMA.bic", "ARIMA.aic", "ARI.bic", 
                                                   "ARI.aic", "ARFIMA.bic", "ARFIMA.aic", "RW")


cbind(h, apply(RESULTSpmse, c(1, 2), function(x) round(x, digits = 3)) < 0.05)
cbind(h, apply(RESULTSpabs, c(1, 2), function(x) round(x, digits = 3)) < 0.05)



xtable::xtable(RESULTSmse, digits = c(0, 0, rep(2, 7)))
xtable::xtable(RESULTSabs, digits = c(0, 0, rep(2, 7)))

rownames(RESULTSmse) <- rownames(RESULTSabs) <- h
RESULTS <- cbind(RESULTSmse[,1],
      RESULTSmse[,2], RESULTSabs[,2],
      RESULTSmse[,3], RESULTSabs[,3],
      RESULTSmse[,4], RESULTSabs[,4],
      RESULTSmse[,5], RESULTSabs[,5],
      RESULTSmse[,6], RESULTSabs[,6],
      RESULTSmse[,7], RESULTSabs[,7],
      RESULTSmse[,8], RESULTSabs[,8])

xtable::xtable(RESULTS, digits = c(0, 0, rep(2, 14)))

calcMSE(y.pred[(372 + hh):length(y), hh], y[(372 + hh):length(y)])
calcMSE(y.arima.111[(372 + hh):length(y), hh], y[(372 + hh):length(y)])
calcMSE(y.arima.111[(372 + hh):length(y), hh], y[(372 + hh):length(y)])

# What benchmarks would we like to have?
# obvious ones: auto AR, auto ARMA, auto ARIMA, RW
# 
# get the orders
