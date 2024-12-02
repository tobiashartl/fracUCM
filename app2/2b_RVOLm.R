#Packages
gc()
rm(list = ls())
library(fUCpack)
library(dplyr)
library(parallel)
library(CFFpack)
library(tidyquant)

#install_github("https://github.com/tobiashartl/CFFpack")
data.name <- "SPSPREAD"

getSymbols("^GSPC", from = '1928-01-01',
           to = "2024-10-31",warnings = FALSE,
           auto.assign = TRUE)

# Wd, etc
setwd("~/R/")
source("./help functions/fUC_arma_approx.R")

# Calculate monthly volatility:
DATE = (index(GSPC))
data <- GSPC %>%
    as.data.frame %>%
    mutate(DATE = DATE) %>%
    mutate(DATE.m = as.yearmon(DATE))
    

data.m <- data %>%
    group_by(DATE.m) %>%
    mutate(MVOL.open.adj = ((GSPC.Open - mean(GSPC.Open))/mean(GSPC.Open))^2,
           MVOL.close.adj = ((GSPC.Close - mean(GSPC.Close))/mean(GSPC.Close))^2,
           MVOL.open = (GSPC.Open - mean(GSPC.Open))^2,
           MVOL.close = (GSPC.Close - mean(GSPC.Close))^2,
           MVOL.OvC = (GSPC.Open - GSPC.Close)^2,
           MVOL.HvL = (GSPC.High - GSPC.Low)^2) %>%
    summarize(MVOL.open_m = mean(MVOL.open),
              MVOL.close_m = mean(MVOL.close),
              MVOL.open_m.adj = mean(MVOL.open.adj),
              MVOL.close_m.adj = mean(MVOL.close.adj),
              MVOL.hvl_m = mean(MVOL.HvL),
              MVOL.ovc_m = mean(MVOL.OvC)
              ) %>%
    mutate(lMVOL.o = log(MVOL.open_m)*100,
           lMVOL.c = log(MVOL.close_m)*100,
           lMVOL.o.adj = log(MVOL.open_m.adj)*100,
           lMVOL.c.adj = log(MVOL.close_m.adj)*100,)

plot(data.m$lMVOL.c, type="l")
plot(data.m$lMVOL.c.adj, type="l")

save(data.m, file = "mRVOL.RData")


R <- 1000
corr <- TRUE
START = 2
diffuse = T
data(list=paste("d2arma",3,3,"_n",2082,sep =""))
det <- c("const", "trend", "frac")
dum <- c(FALSE)
p   <- 1:10
grid <- expand.grid(p, det, dum)
y <- data$lRV

EW(y, alpha = .7, type = 2, interval = c(0, 2))

c(log(var(frac_diff(y, 2))), log(var(frac_diff(y, 1.5))), log(var(frac_diff(y, 1))), 
  # variance in fractional differences: about -12
  log(var(frac_diff(y, 0.5))), log(var(y)))  # 0: maximum variance to sample from 
eta <- NULL

for(i in 1:nrow(grid)){
    p <- grid[i, 1]
    det <- grid[i, 2]
    dummies <- grid[i, 3]
    if(det == "none"){
        det.type = FALSE
    }else{
        det.type = det
    }
    ### get starting values
    #if(file.exists(paste("./Applications/app2/", data.name, "/1_start_", det, "_p", p, ".RData", sep=""))) next
    
    
    set.seed(42)
    d  <- runif(R, min = 1/2, max = 3/2)
    
    # Q
    if(corr){
        s1s <- runif(R,-10, 10)
        s2s <- runif(R, -5, 10)
        rs  <- runif(R, 0, 0)
        Sigvec <- sapply(1:R, function(x){
            Sig <- diag(sqrt(exp(c(s1s[x],s2s[x]))))
            Cor <- matrix(c(1,rs[x],rs[x],1),2,2)
            return(mlogvech(Sig%*%Cor%*%Sig))
        }) 
        
        
        
    }else{
        Sigma <- cbind(runif(R, -10, 10),
                       runif(R, -5, 10))
        Sigvec <- t(log(Sigma))
    }
    
    
    if(p > 0){
        # Starting values for Phi
        as <- matrix(NA,R,p)
        as[,1] <- runif(R,0.5, 1)
        if (p>1)
        {
            as[,2:p] <- runif(R*(p-1),-0.7, 0.7)
            as <- t(apply(as,1,arfima::PacfToAR))
        }
        A <- -as
    }else{
        A <- NULL
    }
    if(p>0){
        if(!is.null(eta)){
            ev <- apply(-A, 1, function(x) max(abs(toComp(x)$eigv)))
            while(any(ev > 1-2*eta)){
                
                A[ev > 1-2*eta, 1] <-  runif(sum(ev > 1-2*eta),0.7, 1)
                if(p>1){
                    A[ev > 1-2*eta, 2:p] <- runif(sum(ev > 1-2*eta)*(p-1),-0.7, 0.7)
                    A[ev > 1-2*eta, ] <- -t(apply(A[ev > 1-2*eta, , drop=F],1,arfima::PacfToAR))
                }
                ev <- apply(-A, 1, function(x) max(abs(toComp(x)$eigv)))
            }
            
        } 
    }
    
    
    
    
    START.val <- cbind(d, t(Sigvec), A)#, Irrv)
    optfn <- function(par0, y){ tryCatch({
        est <- optim(par=par0, method = "BFGS", 
                     fn = fUC_opt_ML_ARMA_approx, deterministics = det.type,  
                     y=y,  corr=corr, START = START, d.int = c(1/2, 2),
                     nulim = c(0, Inf), control = list(trace = 1, REPORT = 1),
                     pq=c(p, 0),penalty.corr=F,
                     diffuse = diffuse, return.det=F, nu.opt=F, Q.trans = "mlv", flip=F, neg=F,
                     seas=dummies, period=12, irregular = F) # custom settings for ML
        Q <- mlogvech2mat(est$par[2:4])
        
        
        cat("ll = ", est$value, ", par = ", est$par, ", corr = ", cov2cor(Q)[2,1],"\n")  
        return(c(est$value, est$par))
    }, 
    error = function(e) return(rep(NA, 1 + ncol(START.val))))
    }
    
    #optfn(START.val[1,], y)
    
    system.time({
        cl <- makeCluster(32)
        clusterExport(cl, ls())
        clusterEvalQ(cl, library(fUCpack))
        clusterEvalQ(cl, library(CFFpack))
        
        RESULTS <- parSapply(cl, 1:R, function(j) optfn(START.val[j,], y))
        stopCluster(cl)
    }
    )
    
    RESULTS <- t(RESULTS)
    RESULTS.win <- RESULTS[which.min(RESULTS[,1]), ]
    
    save(RESULTS, file = paste("./Applications/app2/", data.name, "/1_start_", det, "_p", p, ".RData", sep=""))
    
    cat("Results: det = ", det, ", p = ", p, ", seas = ", dummies, 
        ", diffuse = ", diffuse, "\n", "LL = ", RESULTS.win[1], ", theta = ", 
        RESULTS.win[-1], ", corr = ", cov2cor(mlogvech2mat(RESULTS.win[3:5]))[2,1], "\n")
}


# choose the appropriate specification
for(det in c("const")){
    for(p in c(1:10)){
        cat("Type = ", det, ", p = ", p, "\n")
        load(file = paste("./Applications/app2/", data.name, "/1_start_", det, "_p", p, ".RData", sep=""))
        RESULTS <- RESULTS[order(RESULTS[,1]),] %>%
            as.data.frame %>%
            subset(!is.na(V1)) 
        RESULTS$correlation <- apply(RESULTS[,3:5], 1, function(x) cov2cor(mlogvech2mat(x))[2,1])
        
        RESULTS.ar <- RESULTS[, 6:(5+p), drop=F]
        stable <- apply(RESULTS.ar, 1, function(x) all(abs(toComp(x)$eigv) < 0.99))
        RESULTS <- RESULTS[stable, ]
        
        
        if(det == "none"){
            ndet = 0
        }else{
            if(det %in% c('frac', 'trend', 'const')){
                ndet = 1
            }else{
                ndet = 2
            }
        }
        k <- ncol(RESULTS) - 2 + ndet
        
        RESULTS$BIC         <- k * log(length(y)) + 2*RESULTS[,1]
        RESULTS$AIC         <- 2 * k              + 2*RESULTS[,1]
        
        RESULTS[1:3,] %>% print
    }
}



### pseudo out-of-sample forecast experiment; recursive window
h.ahead <- 30*12 # 30 years ahead
n.start <- n.last
det <- "trend"



p <- 3
load(file = paste("./Applications/app3/", data.name, "/1_start_", det, "_p", p, ".RData", sep=""))
theta <- RESULTS[which.min(RESULTS[,1]),-1]
det <- fUC_opt_ML_ARMA_approx(theta, y, nulim = c(0, Inf), 
                              corr=TRUE, deterministics = det, START=2, diffuse = T, 
                              return.det = TRUE, nu.opt=F, d.int = c(0, 2.5), 
                              eta=NULL, Q.trans = "mlv", flip=F, neg=F, pq=c(p, 0), penalty.corr=F, 
                              seas=F, period=12, irregular=F)
n <- length(y)
Z <- 1:n
dettrend <- cbind(Z) %*% det
plot(y, type="l")
lines(dettrend, col="2")
#det=0

TC <- fUC_KS_ARMA_approx(theta, y - dettrend, nulim = c(0, Inf), quiet=F,
                         corr=TRUE, deterministics = "trend", START=1, diffuse = T, 
                         return.det = F, nu.opt=F, d.int = c(0, 2.5), 
                         eta=NULL, Q.trans = "mlv", flip=F, neg=F, pq=c(p, 0), penalty.corr=F, 
                         seas=F, period=12, irregular=F)
plot(TC$x, type="l")
plot(TC$x + dettrend)
lines(y, col="2")
plot(ts(TC$c, start = c(1850, 1), frequency = 12))
abline(h=0)






### exact likelihood
fUC_opt_ML(theta, y, nulim = c(0, Inf), quiet=F,
           corr = TRUE, deterministics = "frac",
           START = START, d.int = c(0, 2),
           diffuse = diffuse, nu.opt = F, eta = NULL, Q.trans = "mlv")



library("optimParallel")
cl <- makeCluster(6, outfile = "")     # set the number of processor cores
clusterExport(cl, ls())
clusterEvalQ(cl, library(ucminf))
clusterEvalQ(cl, library(fUCpack))
clusterEvalQ(cl, library(CFFpack))
setDefaultCluster(cl=cl) # set 'cl' as default cluster


est.final <- optimParallel(theta, fn = fUC_opt_ML, 
                           method = "BFGS", y=y, quiet=F,
                           nulim = c(0, Inf), corr = TRUE, deterministics = "frac", 
                           START = START, d.int = c(0, 2),
                           diffuse = diffuse, nu.opt = F, eta = NULL, Q.trans = "mlv")



k <- length(est.final$par) + 1
(BIC         <- k * log(length(y)) + 2*est.final$value)
(AIC          <- 2 * k              + 2*est.final$value)




saveRDS(est.final, file = "./app/NOAA_SST_est_frac_p4.RDS")
est.final$Hessian <- optimHess(est.final$par, 
                               fn = fUC_opt_ML, 
                               y=y, quiet=F,
                               nulim = c(0, Inf), corr = TRUE, deterministics = "frac", 
                               START = START, d.int = c(0, 2),
                               diffuse = diffuse, nu.opt = F, eta = NULL, Q.trans = "mlv")
saveRDS(est.final, file = "./app/NOAA_SST_est_frac_p4_hess.RDS")


