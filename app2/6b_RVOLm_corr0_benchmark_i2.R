#Packages
gc()
rm(list = ls())
library(fUCpack)
library(dplyr)
library(parallel)
library(CFFpack)


# Wd, etc
setwd("/Users/tobias/Dokumente/Projekte/filtering unknown persistence/R/code")
source("./help functions/fUC_arma_approx.R")
source("./help functions/UC_i1.R")


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
R <- 100
START = 2
diffuse = T
seas = F
det = "none"
eta <- NULL

for(p in c(1:12)){
    
    ### get starting values
    set.seed(42)
    
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
        Sigma <- cbind(runif(R, -30, 10),
                       runif(R, -30, 10))
        Sigvec <- t((Sigma))
    }
    
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
    
    
    
    
    START.val <- cbind(t(Sigvec), A)#, Irrv)
    RESULTS <- matrix(NA, nrow = R, ncol = ncol(START.val)+1)
    optfn <- function(par0, y){ tryCatch({
        est <- optim(par=par0, method = "BFGS", 
                     fn = UC_i2_opt_ML, deterministics = FALSE,  
                     y=y,  corr=corr, START = START,
                     nulim = c(0, Inf), quiet=T,
                     pq=c(p, 0), penalty.corr=F,
                     diffuse = diffuse, return.det=F, nu.opt=F, Q.trans = "mlv", flip=F, neg=F,
                     seas=seas, period=12, irregular = F) # custom settings for ML
        Q <- diag(exp(est$par[1:2]))
        
        
        cat("ll = ", est$value, ", par = ", est$par, ", corr = ", cov2cor(Q)[2,1],"\n")  
        return(c(est$value, est$par))
    }, 
    error = function(e) return(rep(NA, 1 + ncol(START.val))))
    }
    
    #optfn(START.val[1,], y)
    
    system.time({
        cl <- makeCluster(6)
        clusterExport(cl, ls())
        clusterEvalQ(cl, library(ucminf))
        clusterEvalQ(cl, library(fUCpack))
        clusterEvalQ(cl, library(CFFpack))
        
        RESULTS <- parSapply(cl, 1:R, function(j) optfn(START.val[j,], y))
        stopCluster(cl)
    }
    )
    
    
    RESULTS <- t(RESULTS)
    save(RESULTS, file = paste("./app2/Benchmark_ML_i2_p", p, "_none_corr0.RData", sep=""))
    RESULTS.win <- RESULTS[which.min(RESULTS[,1]), ]
    
    
    cat("Results: det = ", det, ", p = ", p, "\n", "LL = ", RESULTS.win[1], ", theta = ", 
        RESULTS.win[2], ", corr = ", 0, "\n")
}



for(p in c(1:12)){
    cat("p = ", p, "\n")
    load(file = paste("./app2/Benchmark_ML_i2_p", p, "_none_corr0.RData", sep=""))
    RESULTS <- RESULTS[order(RESULTS[,1]),] %>%
        as.data.frame %>%
        subset(!is.na(V1)) 
    RESULTS$correlation <- 0
    
    RESULTS.ar <- RESULTS[, 4:(3+p), drop=F]
    #stable <- apply(RESULTS.ar, 1, function(x) all(abs(toComp(-x)$eigv) < 0.99))
    #RESULTS <- RESULTS[stable, ]
    
    
    
    ndet = 0
    
    k <- ncol(RESULTS) - 2 + ndet
    
    RESULTS$BIC         <- k * log(length(y)) + 2*RESULTS[,1]
    RESULTS$AIC         <- 2 * k              + 2*RESULTS[,1]
    
    RESULTS[1:1,] %>% print
    
}

p <- 5
load(file = paste("./app2/Benchmark_ML_i2_p", p, "_none_corr0.RData", sep=""))
theta <- RESULTS[which.min(RESULTS[,1]),-1]


# confindence bands
hessian <- optimHess(theta, 
                     fn = UC_i2_opt_ML, deterministics = FALSE,  
                     y=y,  corr=corr, START = START,
                     nulim = c(0, Inf), quiet=T,
                     pq=c(p, 0), penalty.corr=F,
                     diffuse = diffuse, return.det=F, nu.opt=F, Q.trans = "mlv", flip=F, neg=F,
                     seas=seas, period=12, irregular = F)

J <- numDeriv::jacobian(exp, theta[1:2])
cov0 <- solve(hessian)
Cov <- bdiag(J[,], diag(p))%*%cov0%*%t(bdiag(J[, ], diag(p)))
theta.trans <- c((exp(theta[1:2]))[], theta[-(1:2)])
se.trans <- sqrt(abs(diag(Cov)))
EST <- cbind(c(NA, theta.trans), c(NA, se.trans)) 
colnames(EST) <- c("par", "se") 
rownames(EST) <- c("d", "Q11", "Q22", paste("ar_", 1:p, sep=""))
EST
EST[,1] / EST[,2]





TC <- fUC_KS_ARMA_approx(c(2, theta), y , nulim = c(0, Inf), 
                         corr=F, deterministics = "none", START=1, diffuse = T, 
                         return.det = F, nu.opt=F, d.int = c(0, 2.5), 
                         eta=NULL, Q.trans = "mlv", flip=F, neg=F, pq=c(p, 0), penalty.corr=F, 
                         seas=F, period=12, irregular=F)
plot(TC$x, type="l")
#plot(TC$x + dettrend)
lines(y, col="2")
plot(ts(TC$c, start = c(1960, 1), frequency = 12))
abline(h=0)

ll <- fUC_opt_ML(theta = c(2, theta), 
                 y=y, nulim =c(0, Inf), quiet = TRUE, corr=F, deterministics = F, 
                 START=2, diffuse =TRUE, return.det=F, nu.opt = F, d.int = c(0, 2.5), 
                 eta = NULL, Q.trans = "mlv", flip =F, neg=F)



k <- length(theta)+1
(BIC          <- k * log(length(y)) + 2*ll)
(AIC          <- 2 * k              + 2*ll)
Q <- diag(exp(theta[2:3]))
SSR <- (fUC_comp((y), 1, Q, theta[-(1:2)], corr=T)$v[-1])^2 %>% sum(.)
add <- c(Q[2,2]/Q[1,1], Q[2,1]/Q[1,1], cov2cor(Q)[2,1], ll, SSR, AIC, BIC)

EST <- rbind(EST, cbind(add, NA))
rownames(EST) <- c("d", "Q11", "Q22", paste("ar_", 1:p, sep=""),
                   "nu_1", "nu_2", "rho", "ll", "CSS", "AIC", "BIC")
colnames(EST)<- c("par", "se")
saveRDS(EST, file = "./app2/Benchmark_i2_ML_Results.RDS")
data.plot <- data.frame(
    y=y,
    trend = TC$x, 
    x = TC$x, 
    cycle = TC$c,
    time = seq(from = as.Date("1960-01-01"), to = as.Date("2024-10-01"),  by = "month"))
saveRDS(data.plot, file = "./app2/Benchmark_i2_ML_TC.RDS")
plot(data.plot$y)
lines(data.plot$x)
