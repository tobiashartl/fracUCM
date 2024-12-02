#Packages
gc()
rm(list = ls())
library(fUCpack)
library(dplyr)
library(parallel)
library(CFFpack)
#install_github("https://github.com/tobiashartl/CFFpack")


# Wd, etc
setwd("~/R/")
source("./help functions/fUC_arma_approx.R")

data.name <- "PCEPI"

# Load inflation data
data<- read.csv(file = paste("./Applications/app3/", data.name, ".csv", sep=""), skip = 0)[, ] %>%
    as.data.frame() 
data <- data %>%
    mutate(DATE = as.Date(DATE)) %>%
    mutate(ly = log(data[, data.name]))  %>%
    filter(DATE >= as.Date("1960-01-01"))


# set a starting date for the inflation application: We choose 1989-12-01
n.last <- as.Date("1989-12-01")
data <- data %>% filter(DATE <= n.last)
plot(data$ly ~ data$DATE, type="l")
plot(diff(data$ly) ~ data$DATE[-1], type="l")

y <- data$ly
(n <- length(y))

R <- 1000
corr <- TRUE
START = 2
diffuse = T
data(list=paste("d2arma",3,3,"_n",300,sep =""))
det <- c("trend")
dum <- c(FALSE)
p   <- 1:12
grid <- expand.grid(p, det, dum)

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
    
    #if(file.exists(paste("./Applications/app3/", data.name, "/1_start_", det, "_p", p, ".RData", sep=""))) next
    
    
    set.seed(42)
    d  <- runif(R, min = 1, max = 3/2)
    
    # Q
    if(corr){
        s1s <- runif(R,-25, -3)
        s2s <- runif(R, -15, 0)
        rs  <- runif(R, 0, 0)
        Sigvec <- sapply(1:R, function(x){
            Sig <- diag(sqrt(exp(c(s1s[x],s2s[x]))))
            Cor <- matrix(c(1,rs[x],rs[x],1),2,2)
            return(mlogvech(Sig%*%Cor%*%Sig))
        }) 
        
        
        
    }else{
        Sigma <- cbind(runif(R, -25, -3),
                       runif(R, -15, 0))
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
                     y=y,  corr=corr, START = START, d.int = c(1/2, 5/2),
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
    
    save(RESULTS, file = paste("./Applications/app3/", data.name, "/1_start_", det, "_p", p, ".RData", sep=""))
    
    cat("Results: det = ", det, ", p = ", p, ", seas = ", dummies, 
        ", diffuse = ", diffuse, "\n", "LL = ", RESULTS.win[1], ", theta = ", 
        RESULTS.win[-1], ", corr = ", cov2cor(mlogvech2mat(RESULTS.win[3:5]))[2,1], "\n")
}


# choose the appropriate specification
for(det in c("trend")){
    for(p in c(1:12)){
        cat("Type = ", det, ", p = ", p, "\n")
        load(file = paste("./Applications/app3/", data.name, "/1_start_", det, "_p", p, ".RData", sep=""))
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


