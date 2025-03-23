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



library(forecast)
# arma
auto.arima(y[1:372], max.p=12, max.q=12, seasonal=F, ic = "bic") # p=q=1
auto.arima(y[1:372], max.p=12, max.q=12, seasonal=F, ic = "aic") # p=2, q=3
auto.arima(y[1:372], max.p=12, max.q=12, seasonal=F, ic = "bic", d=0) # p=q=1
auto.arima(y[1:372], max.p=12, max.q=12, seasonal=F, ic = "aic", d=0) # p=2, q=1
# ar
auto.arima(y[1:372], max.p=12, max.q=0, seasonal=F, ic = "bic") # p=2, d=1
auto.arima(y[1:372], max.p=12, max.q=0, seasonal=F, ic = "aic") # p=4 d=1
auto.arima(y[1:372], max.p=12, max.q=0, seasonal=F, ic = "bic", d=0) # p=2, d=0
auto.arima(y[1:372], max.p=12, max.q=0, seasonal=F, ic = "aic", d=0) # p=3 d=0

# arfima
arfimabic <- function(a, b, AIC = F, d=0){
    mod <- arfima::arfima(y[1:372], order = c(a, d, b))
    if(AIC){
        return(as.numeric(AIC(mod)))
    }else{
        return(as.numeric(BIC(mod)))
    }
}
grid <- expand.grid(0:12, 0:12)
res1 <- cbind(grid, 
      sapply(1:nrow(grid), function(j) arfimabic(grid[j,1], grid[j, 2])),
      sapply(1:nrow(grid), function(j) arfimabic(grid[j,1], grid[j, 2], AIC =TRUE))
)

grid2 <- expand.grid(0:12, 0:12)
res2 <- cbind(grid2, 
      sapply(1:nrow(grid2), function(j) arfimabic(grid2[j,1], grid2[j, 2], d=1)),
      sapply(1:nrow(grid2), function(j) arfimabic(grid2[j,1], grid2[j, 2], d=1, AIC =TRUE))
)


# d = 0
### BIC: (0, 0); 3089
### AIC: (6, 6); 3075
### d= 1
### BIC: (0, 0); 3079 x
### AIC: (4, 3); 3065 x





# models to be estimated:
# ARIMA(1, 1, 1) -> BIC for ARIMA
# ARIMA(2,1,3) -> AIC for ARIMA
# ARIMA(1, 0, 1) -> BIC for ARIMA with d=0
# ARIMA(2, 0, 1) -> AIC for ARIMA with d =0
# AR(2, 1, 0) -> BIC for ARI
# AR(4, 1, 0) -> AIC for ARI
# AR(2, 0, 0) -> BIC for ARI with d =0
# AR(3, 0, 0) -> AIC for ARI with d = 0
# ARFIMA(0, d, 0) -> BIC for ARIMA with integer-diff = 0
# ARFIMA(6, d, 6) -> AIC for ARIMA with integer-diff = 0
# ARFIMA(0, d, 0) -> BIC for ARIMA with integer-diff = 1
# ARFIMA(4, d, 3) -> AIC for ARIMA with integer-diff = 1

re_est_seq <- seq(372, length(y), by = 1)
for(j in 372:length(y)){
    cat("Iteration ", j, " of ", length(y), "\n")
    if(file.exists(file = paste("./app2/iterative/ben/pred_t", j, ".RDS", sep=""))) next 

    y_tilde <- y[1:j]
    
    # Estimate all models
    ### ARIMA
    arima.111 <- arima(y_tilde, order = c(1, 1, 1))
    arima.213 <- arima(y_tilde, order = c(2, 1, 3))
    arima.101 <- arima(y_tilde, order = c(1, 0, 1))
    arima.201 <- arima(y_tilde, order = c(2, 0, 1))
    
    ### ARI
    ari.210 <- arima(y_tilde, order = c(2, 1, 0))
    ari.410 <- arima(y_tilde, order = c(4, 1, 0))
    ari.200 <- arima(y_tilde, order = c(2, 0, 0))
    ari.300 <- arima(y_tilde, order = c(3, 0, 0))
    
    ### ARFIMA
    sink("/dev/null") 
    arfima.0d0 <- arfima::arfima(y_tilde, order = c(0, 0, 0))
    arfima.6d6 <- arfima::arfima(y_tilde, order = c(6, 0, 6))
    arfima.0d10 <- arfima::arfima(y_tilde, order = c(0, 1, 0))
    arfima.4d13 <- arfima::arfima(y_tilde, order = c(4, 1, 3))
    sink() 
    # do the predictions
    pred.arima.111 <- c(y_tilde, as.numeric(predict(arima.111, h)$pred))
    pred.arima.213 <- c(y_tilde, as.numeric(predict(arima.213, h)$pred))
    pred.arima.101 <- c(y_tilde, as.numeric(predict(arima.101, h)$pred))
    pred.arima.201 <- c(y_tilde, as.numeric(predict(arima.201, h)$pred))
    pred.ari.210 <- c(y_tilde, as.numeric(predict(ari.210, h)$pred))
    pred.ari.410 <- c(y_tilde, as.numeric(predict(ari.410, h)$pred))
    pred.ari.200 <- c(y_tilde, as.numeric(predict(ari.200, h)$pred))
    pred.ari.300 <- c(y_tilde, as.numeric(predict(ari.300, h)$pred))
    pred.arfima.0d0 <- c(y_tilde, as.numeric(predict(arfima.0d0, h)[[1]]$Forecast))
    pred.arfima.6d6 <- c(y_tilde, as.numeric(predict(arfima.6d6, h)[[1]]$Forecast))
    pred.arfima.0d10 <- c(y_tilde, as.numeric(predict(arfima.0d10, h)[[1]]$Forecast))
    pred.arfima.4d13 <- c(y_tilde, as.numeric(predict(arfima.4d13, h)[[1]]$Forecast))
    
   
    
    results <- data.frame(pred.arima.111 = pred.arima.111, 
                          pred.arima.213 = pred.arima.213, 
                          pred.arima.101 = pred.arima.101, 
                          pred.arima.201 = pred.arima.201,
                          pred.ari.210 = pred.ari.210,
                          pred.ari.410 = pred.ari.410,
                          pred.ari.200 = pred.ari.200,
                          pred.ari.300 = pred.ari.300,
                          pred.arfima.0d0 = pred.arfima.0d0,
                          pred.arfima.6d6 = pred.arfima.6d6,
                          pred.arfima.0d10 = pred.arfima.0d10,
                          pred.arfima.4d13 = pred.arfima.4d13,
                          y.true = y[1:length(pred.arfima.4d13)]
                          )
    saveRDS(results, file = paste("./app2/iterative/ben/pred_t", j, ".RDS", sep=""))
}


