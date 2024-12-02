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

# Wd, etc
setwd("~/R/")


# Settings
n <- c(100, 200, 300)
d <- c(0.75, 1, 1.75)
ratio <- c(1, 10, 30, 0)
R <- 1000
nulim = c(1/1000, 1e+06)

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
# calculate variance ratio for grid


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


# corr_grid <- c(tcrossprod(-c(0.7, 0.8, 0.9), sqrt(r.grid)))

for ( i in 1:NROW(setups)){
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    d <- as.numeric(setup[2])
    nu <- nucalc(n, d, setups[i, 3], ar)
    r <- setups[i, 3]
    cat("Iteration ", i, "\n")
    if(file.exists(file = paste("./MC/MC_2/CSS/Sim_R", R, "_n", n, "_d", d, "_r", r, "_corr.RData", sep=""))) next
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
    d.grid <- seq(0.75, 1.75, 0.25)
    nu.grid <- (c(0.02, 1, 100, 1000))
    ar.1.grid <- c(-1.8,  -1.6, -1.4)
    ar2.grid  <- c(0.6,  0.8,  1)
    corr_grid <- c(tcrossprod(-c(0.9, 0.8, 0.7), sqrt(nu.grid)))
    
    gr.start <- expand.grid(d.grid, nu.grid, corr_grid, ar.1.grid, ar2.grid) %>% as.matrix()
    gr.start <- gr.start[sapply(1:nrow(gr.start), function(i) toComp(-c(gr.start[i, 4:5]))$stable) ,]
    
    
    optfn <- function(n, y, x){
        tryCatch({
            # Estimate
            # check grid
            grid.st <- cbind(gr.start, sapply(1:nrow(gr.start), function(i) fUC_opt(gr.start[i, ], corr=TRUE,
                                                                                    y=y, START = 1, 
                                                                                    nulim = c(1/1000, 1e+06),
                                                                                    deterministics = FALSE)))
            
            par0 <- grid.st[which.min(grid.st[,6]),-6]
            
            (est <- optim(par=par0, control = list(), corr=TRUE,
                          fn = fUC_opt, quiet = TRUE, START = 1,
                          y=y, method = "BFGS", nulim = c(1/100, 1e+06), d.int = c(0.5, 2.5),
                          deterministics = FALSE))
            
            # check for corner solutions: d->0, d->2, nu->1/100, nu->100
            j=1
            while(est$par[1] < .55 | est$par[1] > 2.45 | exp(est$par[2]) < .01 | exp(est$par[2]) > 250000){
                if(j > 10) break
                j=j+1
                par0 <- grid.st[order(grid.st[,6]),-6][j , ]
                (est <- optim(par=par0, corr=TRUE,
                              fn = fUC_opt, quiet = TRUE, d.int = c(0.5, 2.5),
                              y=y, method = "BFGS", nulim = c(1/100, 250000),
                              deterministics = FALSE, START = 1))
                
                
            }
            
            if(est$par[1] < .55 | est$par[1] > 2.45 | exp(est$par[2]) < .01 | exp(est$par[2]) > 250000){
                # try primitive starting values
                par0 <- c(2, grid.st[which.min(grid.st[,6]),2:5])
                (est <- optim(par=par0, corr=TRUE,
                              fn = fUC_opt, quiet = TRUE, d.int = c(0.5, 2.5),
                              y=y, method = "BFGS", nulim = c(1/100, 250000),
                              deterministics = FALSE, START = 1))
                
                if(est$par[1] < .55 | est$par[1] > 2.45 | exp(est$par[2]) < .01 | exp(est$par[2]) > 250000){
                    # try primitive starting values
                    par0 <- c(1, grid.st[which.min(grid.st[,6]),2:5])
                    (est <- optim(par=par0, corr=TRUE,
                                  fn = fUC_opt, quiet = TRUE, d.int = c(0.5, 2.5),
                                  y=y, method = "BFGS", nulim = c(1/100, 250000),
                                  deterministics = FALSE, START = 1))
                    
                    
                }
            }
            
            
            par <- est$par
            nu  <- (par[2])
            corr <- par[3]
            ar <- est$par[-(1:3)]
            Q <- matrix(c(1, corr, corr, nu), 2, 2)
            
            KS <- fUC_smooth(y, par[1], Q, ar = ar, corr = TRUE)
            Rsq <- summary(lm(x ~ KS$x))$r.squared
            
            
            # Return results
            results <- c(par[1], (par[2]), par[3], ar, Rsq)
            names(results) <- c("d", "nu", "corr", "ar_1", "ar_2", "Rsq")
            return(results)
        }, error = function(e) return(rep(NA, 6))
        )
        
    }
    
    #optfn(n, y[,3], x[,3])
    
    ### fUC part
    cl <- makeCluster(32)
    clusterExport(cl, c("optfn", "fUC_opt", "y", "fUC_comp", "ma_inf",
                        "embed0", "fUC_smooth", "x", "frac_diff", "lm", "n", "gr.start"))
    clusterEvalQ(cl, library(fUCpack))
    RESULTS <- parSapply(cl, 1:R, function(j){ 
        optfn(n, y[, j], x[, j])})
    stopCluster(cl)
    
    dEW_45 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2.5), alpha = .45))
    dEW_50 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2.5), alpha = .50))
    dEW_55 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2.5), alpha = .55))
    dEW_60 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2.5), alpha = .60))
    dEW_65 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2.5), alpha = .65))
    dEW_70 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2.5), alpha = .70))
    
    RESULTS <- rbind(RESULTS, c(dEW_45), c(dEW_50), c(dEW_55), c(dEW_60), c(dEW_65), c(dEW_70))
    rownames(RESULTS) <- c("d", "nu", "corr", "ar_1", "ar_2", "Rsq", 
                           "d_45", "d_50", "d_55", "d_60", "d_65", "d_70")
    save(RESULTS, file = paste("./MC/MC_2/CSS/Sim_R", R, "_n", n, "_d", d, "_r", r, "_corr.RData", sep=""))
}

