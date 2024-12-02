# Simulation:
# Performance of fUC for parameter estimation
#   - High and normal signal to noise ratios
#   - AR plus non-AR
#   - compare d with other nonparametric estimators
gc()
rm(list = ls())
library(fUCpack)
library(dplyr)
library(parallel)
library(CFFpack)

# Wd, etc
setwd("~/R/")
source("./help functions/fUC_arma_approx.R")


# Settings
n <- c(100, 200, 300)
d <- c(0.75, 1, 1.75)
ratio <- c(1, 10, 30, 0)
R <- 1000
nulim <- c(1/1000, 1e+06)

# Strong persistence
ar <- c(-1.6, 0.8)

# Long-run variance part
cvar <- function(ar){
    p <- length(ar)
    A <- toComp(-ar)$CompMat
    return(matrix(solve(diag(p^2)-(A%x%A))%*%c(1,rep(0,p^2-1)),p,p)[1,1])
}
xvar <- function(n, d){
    pisq <- frac_diff(c(1, rep(0, n-1)), -d)^2
    var  <- cumsum(pisq)
    return(mean(var))
}

setups <- expand.grid(n, d, ratio)
colnames(setups) <- c("n", "d", "ratio")
setups <- setups[order(setups[, "n"]), ]


# calculate the true nu:
nucalc <- function(n, d, ratio, ar){
    if(ratio == 0) return(1)
    var_x <- xvar(n, d)
    var_c <- cvar(ar)
    nu    <- var_x/var_c * (1 / ratio)
    return(nu)
}

nuvec <- sapply(1:nrow(setups), function(j) nucalc(setups[j, 1], setups[j, 2], setups[j, 3], ar))
Q21vec <- -0.8*sqrt(nuvec)
# matrices of Qs
Qmatvech <- (lapply(1:NROW(setups), function(j) 
    mlogvech(matrix(c(1, Q21vec[j], Q21vec[j], nuvec[j]), 2, 2))))

for ( i in 1:NROW(setups)){
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    d <- as.numeric(setup[2])
    nu <- nucalc(n, d, setups[i, 3], ar)
    r <- setups[i, 3]
    cat("Iteration ", i, "\n")
    if(file.exists(file = paste("./MC/MC_2/ARMA/Sim_R", R, "_n", n, "_d", d, "_r", r, "_corr.RData", sep=""))) next
    set.seed(42)
    
    # Generate the data
    Q   <- matrix(c(1, -0.8*sqrt(nu), -0.8*sqrt(nu), nu), 2, 2)
    ETA <- EPS <- matrix(NA, nrow = n, ncol=R)
    
    for(j in 1:R){
        ERR <-   mvtnorm::rmvnorm(n, mean=c(0, 0), sigma = Q)
        ETA[, j] <- ERR[,1]
        EPS[, j] <- ERR[,2]
    }
    
    x <- frac_diff_multi(ETA, -d)
    c <- apply(EPS, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- yy<-x + c
    
    # here comes the grid of starting values
    # here comes the grid of starting values
    d.grid <- seq(0.75, 1.75, 0.25)
    # get a reasonable grid of starting values for sigma:
    
    
    #sig1.grid <- c(-0.5, -0.75, -1)
    #sig2.grid <- c(-0.5, 3, 10)
    #sig12.grid <- c(-1, -0.5, 0)
    ar.1.grid <- c(-1.8,  -1.6, -1.4)
    ar2.grid  <- c(0.6,  0.8,  1)
    
    
    gr.start <- expand.grid(d.grid, Qmatvech, ar.1.grid, ar2.grid) 
    gr.start <- cbind(gr.start[,1], gr.start[,2] %>% unlist() %>% matrix(ncol = 3, byrow = TRUE), gr.start[, -(1:2)]) %>% as.matrix()
    gr.start <- gr.start[sapply(1:nrow(gr.start), function(i) toComp(-c(gr.start[i, 5:6]))$stable) ,]
    
    
    optfn <- function(n, y, x){
        tryCatch({
            # Estimate
            # check grid
            grid.st <- cbind(gr.start, sapply(1:nrow(gr.start), function(i) fUC_opt_ML_ARMA_approx(gr.start[i, ], corr=TRUE,
                                                                                                   y=y, START = 2, pq=c(2,0), 
                                                                                                   penalty.corr=FALSE,
                                                                                                   nulim = nulim, nu.opt=F,
                                                                                                   deterministics = FALSE)))
            
            par0 <- grid.st[which.min(grid.st[,7]),-7]
            
            est <- tryCatch({
                (optim(par=par0, control = list(), corr=TRUE, pq=c(2,0), penalty.corr=FALSE,
                       fn = fUC_opt_ML_ARMA_approx, quiet = TRUE, START = 2, nu.opt=F,
                       y=y, method = "BFGS", nulim = nulim, d.int = c(0.5, 2.5),
                       deterministics = FALSE))
            }, error = function(e) return(NA)
            )
            
            
            if(is.na(est[1])) est <- list(par = NA)
            
            # check for corner solutions: d->0, d->2, nu->1/100, nu->100
            j=1
            while(is.na(est$par[1])  |  est$par[1] < .55 | est$par[1] > 2.45){
                if(j > 10) break
                j=j+1
                par0 <- grid.st[order(grid.st[,7]),-7][j , ]
                est <- tryCatch({
                    (optim(par=par0, control = list(), corr=TRUE, pq=c(2,0), penalty.corr=FALSE,
                           fn = fUC_opt_ML_ARMA_approx, quiet = TRUE, START = 2, nu.opt=F,
                           y=y, method = "BFGS", nulim = nulim, d.int = c(0.5, 2.5),
                           deterministics = FALSE))
                }, error = function(e) return(NA)
                )
                if(is.na(est[1])) est <- list(par = NA)
            }
            
            
            if(is.na(est$par[1])  |  est$par[1] < .55 | est$par[1] > 2.45){
                # primitive starting value: set 1
                par0 <- c(2, grid.st[which.min(grid.st[,7]),2:6])
                est <- tryCatch({
                    (optim(par=par0, control = list(), corr=TRUE, pq=c(2,0), penalty.corr=FALSE,
                           fn = fUC_opt_ML_ARMA_approx, quiet = TRUE, START = 2, nu.opt=F,
                           y=y, method = "BFGS", nulim = nulim, d.int = c(0.5, 2.5),
                           deterministics = FALSE))
                }, error = function(e) return(NA)
                )
                if(is.na(est$par[1]) |  est$par[1] < .55 | est$par[1] > 2.45){
                    # primitive starting value: set 2
                    par0 <- c(1, grid.st[which.min(grid.st[,7]),2:6])
                    est <- tryCatch({
                        (optim(par=par0, control = list(), corr=TRUE, pq=c(2,0), penalty.corr=FALSE,
                               fn = fUC_opt_ML_ARMA_approx, quiet = TRUE, START = 2, nu.opt=F,
                               y=y, method = "BFGS", nulim = nulim, d.int = c(0.5, 2.5),
                               deterministics = FALSE))
                    }, error = function(e) return(NA)
                    )
                }
                if(is.na(est[1])) est <- list(par = NA)
            }
            
            
            
            
            par <- est$par
            Q <- mlogvech2mat(est$par[2:4])
            ar <- est$par[-(1:4)]
            #KF <- fUC_comp(y, par[1], nu, ar)
            #KF2 <- fUC_comp((1:length(y))^par[1], par[1], nu, ar)
            #mu  <- summary(lm(KF$v ~ -1 + KF2$v))$coefficients[1,1]
            #Q <- matrix(c(1, corr, corr, nu), 2, 2)
            
            KS <- fUC_smooth(y, par[1], Q, ar = ar, corr = TRUE)
            Rsq <- summary(lm(x ~ KS$x))$r.squared
            
            
            # Return results
            results <- c(est$value, est$par, Rsq)
            names(results) <- c("ll", "d", "Q11", "Q21", "Q22", "ar_1", "ar_2", "Rsq")
            return(results)
        }, error = function(e) return(rep(NA, 8))
        )
        
    }
    #optfn(n, y[,1], x[,1])
    #optfn(n, y[,917], x[,917])
    
    ### fUC part
    cl <- makeCluster(32)
    clusterExport(cl, c("optfn", "fUC_opt_ML_ARMA_approx", "y", "fUC_KS_ARMA_approx", "ma_inf", "nulim",
                        "embed0", "fUC_KF_approx", "fUC_KS_approx", "x", "frac_diff", "lm", "n", "gr.start"))
    clusterEvalQ(cl, library(fUCpack))
    clusterEvalQ(cl, library(CFFpack))
    
    RESULTS <- parSapply(cl, 1:R, function(j) optfn(n, y[, j], x[, j]))
    stopCluster(cl)
    
    dEW_45 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2.5), alpha = .45))
    dEW_50 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2.5), alpha = .50))
    dEW_55 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2.5), alpha = .55))
    dEW_60 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2.5), alpha = .60))
    dEW_65 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2.5), alpha = .65))
    dEW_70 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2.5), alpha = .70))
    
    RESULTS <- rbind(RESULTS, c(dEW_45), c(dEW_50), c(dEW_55), c(dEW_60), c(dEW_65), c(dEW_70))
    rownames(RESULTS) <- c("ll", "d", "Q11", "Q21", "Q22", "ar_1", "ar_2", "Rsq", 
                           "d_45", "d_50", "d_55", "d_60", "d_65", "d_70")
    save(RESULTS, file = paste("./MC/MC_2/ARMA/Sim_R", R, "_n", n, "_d", d, "_r", r, "_corr.RData", sep=""))
}

