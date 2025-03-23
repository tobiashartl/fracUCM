# Note: Looking at the graph, it is highly likely that a purely measurement error
# component will improve the tc decomposition cosiderably; however, this needs to be
# incorporated into the fUCmodel. 
# Note2: well, likelihood does not really improve 




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

load(file = "./Applications/app2/mRVOL.RData")


data <- data.m
# select first 30 years of data -> 360 obs.
data.m <- data.m %>%
    filter(DATE.m <= as.yearmon("1990-12-01")) 


R <- 100
corr <- FALSE
START = 2
diffuse = T
data(list=paste("d2arma",3,3,"_n",1000,sep =""))
det <- c("none", "const", "trend", "frac")
dum <- c(FALSE)
p   <- 1:12
grid <- expand.grid(p, det, dum)
y <- data.m$lRVOL.close_m
y0 <- y[1]
y <- y-y0
n <- length(y)
irregular = FALSE
eta <- NULL
plot(y, type="l")
int <- sapply(seq(0.3, 0.9, by = 0.01), function(j) EW(y, alpha = j, type ="trend", interval = c(0, 2)))
min(int)
max(int)
for(i in 1:nrow(grid)){
    
    p <- grid[i, 1]
    det <- grid[i, 2]
    dummies <- grid[i, 3]
    
    if(file.exists(paste("./Applications/app2/iterative/start/1_start_", det, "_p", p, "_corr0.RData", sep=""))) next
    if(det == "none"){
        det.type = FALSE
    }else{
        det.type = det
    }
    
    set.seed(42)
    d  <- runif(R, min = 0, max = 2)
    
    # Q
    if(corr){
        s1s <- runif(R,-10, 10)
        s2s <- runif(R, -5, 15)
        rs  <- runif(R, 0, 0)
        Sigvec <- sapply(1:R, function(x){
            Sig <- diag(sqrt(exp(c(s1s[x],s2s[x]))))
            Cor <- matrix(c(1,rs[x],rs[x],1),2,2)
            return(mlogvech(Sig%*%Cor%*%Sig))
        }) 
        
        
        
    }else{
        Sigma <- cbind(runif(R, -10, 10),
                       runif(R, -5, 15))
        Sigvec <- t((Sigma))
    }
    
    Qtest <- sapply(1:R, function(j){
        Q <- exp(Sigvec[,j])
        Q[2]/(Q[1]) < 0.1
    })
    Sigvec[ ,Qtest ] <- Sigvec[2:1, Qtest]
    
    
    
    if(p > 0){
        # Starting values for Phi
        as <- matrix(NA,R,p)
        as[,1] <- runif(R,-1, 1)
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
    
    if(irregular){
        Irrv <- runif(R, -0, 20)
    }else{
        Irrv <- NULL
    }
    
    
    
    START.val <- cbind(d, t(Sigvec), A, Irrv)
    optfn <- function(par0, y){ tryCatch({
        est <- optim(par=par0, method = "BFGS", 
                     fn = fUC_opt_ML_ARMA_approx, deterministics = det.type,  
                     y=y,  corr=corr, START = START, d.int = c(0, 5/2),
                     nulim = c(0, Inf), control = list(trace = 1, REPORT = 1),
                     pq=c(p, 0),penalty.corr=F,
                     diffuse = diffuse, return.det=F, nu.opt=F, Q.trans = "mlv", flip=F, neg=F,
                     seas=dummies, period=12, irregular = irregular) # custom settings for ML
        if(corr){
            Q <- mlogvech2mat(est$par[2:4])
        }else{
            Q <- diag(exp(est$par[2:3]))
        }
        
        
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
    
    save(RESULTS, file = paste("./Applications/app2/iterative/start/1_start_", det, "_p", p, "_corr0.RData", sep=""))
    
    cat("Results: det = ", det, ", p = ", p, ", seas = ", dummies, 
        ", diffuse = ", diffuse, "\n", "LL = ", RESULTS.win[1], ", theta = ", 
        RESULTS.win[-1], ", corr = ", 0, "\n")
}


# choose the appropriate specification
for(det in c("none")){
    for(p in c(1:12)){
        cat("Type = ", det, ", p = ", p, "\n")
        load(file = paste("./Applications/app2/iterative/start/1_start_", det, "_p", p, "_corr0.RData", sep=""))
        RESULTS <- RESULTS[order(RESULTS[,1]),] %>%
            as.data.frame %>%
            subset(!is.na(V1)) 
        RESULTS$correlation <- 0
        RESULTS.ar <- RESULTS[, 5:(4+p), drop=F]
        #stable <- apply(RESULTS.ar, 1, function(x) all(abs(toComp(x)$eigv) < 0.99))
        #RESULTS <- RESULTS[stable, ]
        
        
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

### 1 1716.632 0.5027782 7.983144 6.755345 -0.4269204           0 3456.94 3441.264


# const
det <- "none"
p <- 1
load(file = paste("./Applications/app2/iterative/start/1_start_", det, "_p", p, "_corr0.RData", sep=""))
RESULTS <- RESULTS[order(RESULTS[,1]),] %>%
    as.data.frame %>%
    subset(!is.na(V1)) 
RESULTS.ar <- RESULTS[, 5:(4+p), drop=F]
#stable <- apply(RESULTS.ar, 1, function(x) all(abs(toComp(x)$eigv) < 0.99))
#RESULTS <- RESULTS[stable, ]
theta <- as.numeric(RESULTS[1,-1])
# det <- fUC_opt_ML_ARMA_approx(theta, y, nulim = c(0, Inf), 
#                               corr=corr, deterministics = det, START=2, diffuse = T, 
#                               return.det = TRUE, nu.opt=F, d.int = c(0, 2.5), 
#                               eta=NULL, Q.trans = "mlv", flip=F, neg=F, pq=c(p, 0), penalty.corr=F, 
#                               seas=F, period=12, irregular=F)
n <- length(y)
Z <- rep(0, n)
#Z <- 1:n
#Z <- frac_diff(Z, -theta[1])
#dettrend <- cbind(Z) %*% det
plot(y, type="l")
dettrend <- rep(0, n)
lines(dettrend, col="2")
#det=0

TC <- fUC_KS_ARMA_approx(as.numeric(theta), y - dettrend, nulim = c(0.1, Inf), quiet=F,
                         corr=corr, deterministics = det, START=1, diffuse = T, 
                         return.det = F, nu.opt=F, d.int = c(0, 2.5), 
                         eta=NULL, Q.trans = "mlv", flip=F, neg=F, pq=c(p, 0), penalty.corr=F, 
                         seas=F, period=12, irregular=F)

plot(TC$x, type="l")
plot(y)
lines(TC$x + dettrend, col=2)


plot(TC$x)
lines(y, col="2")


plot(ts(TC$c, start = c(1960, 1), frequency = 12), type = "n")
tis::nberShade()
lines(ts(TC$c, start = c(1960, 1), frequency = 12))

### exact likelihood
fUC_opt_ML(theta, y, nulim = c(0, Inf), quiet=F,
           corr = F, deterministics = F,
           START = START, d.int = c(0, 2),
           diffuse = diffuse, nu.opt = F, eta = NULL, Q.trans = "mlv")

fUC_opt_ML_ARMA_approx(theta, y, nulim = c(0, Inf), quiet=F,
                       corr = F, deterministics = F,  pq=c(p, 0),
                       START = START, d.int = c(0, 2),
                       diffuse = diffuse, nu.opt = F, eta = NULL, Q.trans = "mlv")


p=1
load(file = paste("./Applications/app2/iterative/start/1_start_", det, "_p", p, "_corr0.RData", sep=""))

RESULTS <- RESULTS[order(RESULTS[,1]),] %>%
    as.data.frame %>%
    subset(!is.na(V1)) 
RESULTS.ar <- RESULTS[, 5:(4+p), drop=F]
stable <- apply(RESULTS.ar, 1, function(x) all(abs(toComp(x)$eigv) < 0.99))
RESULTS <- RESULTS[stable, ]
theta <- as.numeric(RESULTS[1,-1])
est.final <- optim(theta, fn = fUC_opt_ML, 
                   method = "BFGS", y=y, nulim = c(0, Inf),
                   corr = F, deterministics = "const", quiet=F,
                   START = START, d.int = c(0, 2),
                   diffuse = diffuse, nu.opt = F, eta = NULL, Q.trans = "mlv")


k <- length(est.final$par) + 1
(BIC         <- k * log(length(y)) + 2*est.final$value)
(AIC          <- 2 * k              + 2*est.final$value)




saveRDS(est.final, file = "./Applications/app2/mRVOL_corr0_p1.RDS")
est.final$Hessian <- optimHess(est.final$par, 
                               fn = fUC_opt_ML, 
                               y=y, nulim = c(0, Inf),
                               corr = F, deterministics = F, quiet=T,
                               START = START, d.int = c(0, 2),
                               diffuse = diffuse, nu.opt = F, eta = NULL, Q.trans = "mlv")
solve(est.final$Hessian)
saveRDS(est.final, file = "./Applications/app2/mRVOL_corr0_p1.RDS")


