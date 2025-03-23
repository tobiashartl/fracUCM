
#Packages
gc()
rm(list = ls())
library(fUCpack)
library(dplyr)
library(parallel)
library(CFFpack)
library(tidyquant)
#devtools::install_github('tobiashartl/fracUCM/packages_newfortran/fUCpack')
#install_github("https://github.com/tobiashartl/CFFpack")
data.name <- "RVOLm"

# Wd, etc
setwd("~/R/")
source("./help functions/fUC_arma_approx.R")
source("./help functions/fUC_predict.R")

load(file = "./Applications/app2/mRVOL.RData")


data <- data.m
R <- 1000
corr <- FALSE
START = 2
diffuse = T
data(list=paste("d2arma",3,3,"_n",1000,sep =""))
det <- c("none")
dum <- c(FALSE)
p   <- 1 # according to BIC
grid <- expand.grid(p, det, dum)
y <- data.m$lRVOL.close_m
y0 <- y[1]
y <- y-y0
n <- length(y)
irregular = FALSE
eta <- NULL
plot(y, type="l")

# grid starting values
load(file = paste("./Applications/app2/iterative/start/1_start_", det, "_p", p, "_corr0.RData", sep=""))
RESULTS <- RESULTS[order(RESULTS[,1]),] %>%
    as.data.frame %>%
    subset(!is.na(V1)) 
RESULTS.ar <- RESULTS[, 5:(4+p), drop=F]
#stable <- apply(RESULTS.ar, 1, function(x) all(abs(toComp(x)$eigv) < 0.99))
#RESULTS <- RESULTS[stable, ]
theta <- RESULTS[1,-1]

### start a loop with optimization in each step
### re-estimation sequence: every 12 observations
h=120
re_est_seq <- seq(372, length(y), by = 1)
for(j in 372:length(y)){
    y_tilde <- y[1:j]
    
    # estimate the TC model
    if(j %in% re_est_seq){
        theta_new <- optim(theta, fn = fUC_opt_ML, 
                           method = "BFGS", y=y_tilde, nulim = c(0, Inf),
                           corr = F, deterministics = F, quiet=T,
                           START = START, d.int = c(0, 5/2),
                           diffuse = F, nu.opt = F, eta = NULL, Q.trans = "mlv")
        theta <- theta_new$par
        # save it
        saveRDS(theta, file = paste("./Applications/app2/iterative/par/theta_t", j, ".RDS", sep=""))
        #theta <- readRDS(file = paste("./Applications/app2/iterative/par/theta_t", j, ".RDS", sep=""))
        cat("Iteration ", j, ": Theta = ", theta, "\n")
    }
    
    # generate the forecasts
    Q <- diag(exp(theta[2:3]))
    TC <- fUC_smooth(y_tilde, theta[1], Q, theta[-(1:3)], corr = T)
    
    
    # do the forecasts
    x.hat <- trend_predict(TC$x, theta[1], h=h)
    c.hat <- cycle_predict(TC$c, theta[-(1:3)], h=h)
    plot(y)
    lines(x.hat + c.hat, col=2)
    results <- data.frame(x.hat = x.hat, 
                          c.hat = c.hat, 
                          y.hat = x.hat + c.hat, 
                          y.true = y[1:length(x.hat)])
    saveRDS(results, file = paste("./Applications/app2/iterative/par/pred_t", j, ".RDS", sep=""))
}

trend_predict <- function(x, d, h){
    n       <- length(x)
    ar_coef <- -frac_diff(c(1, rep(0, n+h-1)), d)[-1]
    x.pred <- c(x, rep(NA, h))
    for(ii in 1:h){
        x.pred[n+ii] <- x.pred[1:(n+ii-1)]%*%(ar_coef[1:(n+ii-1)][(n+ii-1):1])
    }
    return(x.pred)
}

cycle_predict <- function(x, ar, h){
    n       <- length(x)
    ar_coef <- c(-ar, rep(0, n+h-length(ar)-1))
    x.pred <- c(x, rep(NA, h))
    for(ii in 1:h){
        x.pred[n+ii] <- x.pred[1:(n+ii-1)]%*%(ar_coef[1:(n+ii-1)][(n+ii-1):1])
    }
    return(x.pred)
}



# load the different d estimates
dload <- function(j){
    theta <- load(paste("./Applications/app2/iterative/par/theta_t", j, ".RDS" , sep =""))
    return(theta[1])
}
sapply(372:length(y), dload)



