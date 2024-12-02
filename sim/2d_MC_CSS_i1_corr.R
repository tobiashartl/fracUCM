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
library(KFAS)

# Wd, etc
setwd("/Users/tobias/Dokumente/Projekte/filtering unknown persistence/R/code")
source("./help functions/KF.R")



# Settings
n <- c(100, 200, 300)
d <- c(0.75, 1, 1.75)
ratio <- c(1, 10, 30, 0, 0.1)
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


for ( i in 1:NROW(setups)){
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    d <- as.numeric(setup[2])
    nu <- nucalc(n, d, setups[i, 3], ar)
    r <- setups[i, 3]
    
    set.seed(42)
    
    rho <- -0.8
    Q <- matrix(c(1, rho*sqrt(1)*sqrt(nu), rho*sqrt(1)*sqrt(nu), nu), 2, 2)
    
    # Generate the data
    U <- sapply(1:R, function(x) mvtnorm::rmvnorm(n, mean=c(0, 0), sigma = Q))
    x <- frac_diff_multi(U[1:n,], d=-d)
    u <- U[-(1:n),]
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- x + c
    
    nu.grid <- (c(0.01, 1, 100, 1000))
    ar.1.grid <- c(-1.8,  -1.6, -1.4)
    ar2.grid  <- c(0.6,  0.8,  1)
    corr_grid <- c(tcrossprod(-c(0.7, 0.8, 0.9), sqrt(nu.grid)))
    
    gr.start <- expand.grid(nu.grid, corr_grid, ar.1.grid, ar2.grid) %>% as.matrix()
    gr.start <- gr.start[sapply(1:nrow(gr.start), function(i) toComp(-c(gr.start[i, 3:4]))$stable) ,]
    
    
    optfn <- function(n, y, x){
        # Estimate
        tryCatch({
            grid.st <- cbind(gr.start, sapply(1:nrow(gr.start), function(i) UC_opt_KF_i1(gr.start[i, ], corr=TRUE,
                                                                                         y=y, START = 1, ll=FALSE,
                                                                                         nulim = nulim,
                                                                                         deterministics = FALSE)))
            
            par0 <- grid.st[which.min(grid.st[,5]),-5]
            
            
            est <- optim(par=par0, START=1,
                         fn = UC_opt_KF_i1, ll = FALSE, corr = TRUE,
                         y=y, method = "BFGS", nulim = nulim)
            
            
            par <- est$par
            nu  <- matrix(c(1, est$par[2], est$par[2], est$par[1]), 2, 2)
            ar  <- par[-(1:2)]
            KS <- fUC_smooth(y, 1, nu, ar = -ar, corr=TRUE)
            KF <- fUC_comp(y, 1, nu, ar=-ar, corr=TRUE)
            
            # Calculate Rsq, etc
            SSR <- mean((x - KS$x)^2)
            SST <- mean((x - mean(x))^2)
            Rsq <- summary(lm(x ~ KS$x))$r.squared
            
            # Return results
            results <- c(nu[2,2], nu[1,2], ar, SSR, Rsq)
            names(results) <- c("sigma_eps", "corr", "ar_1", "ar_2","SSR", "Rsq")
            
            
            return(results)
        }, error = function(e) return(rep(NA, 6))
        )
        
    }
    ### fUC part
    cl <- makeCluster(32)
    clusterExport(cl, c("optfn", "UC_opt_KF_i1", "y", "fUC_comp", "ma_inf", "nulim",
                        "embed0", "fUC_smooth", "x", "frac_diff", "gr.start"))
    clusterEvalQ(cl, library(fUCpack))
    clusterEvalQ(cl, library(KFAS))
    
    RESULTS <- parSapply(cl, 1:R, function(j) optfn(n, y[, j], x[, j]))
    stopCluster(cl)
    
    
    
    rownames(RESULTS) <- c("sigma_eps", "corr", "ar_1", "ar_2", "SSR", "Rsq")
    save(RESULTS, file = paste("./MC/MC_2/CSS_INTEGER/Sim_R", R, "_n", n, "_d", d, "_r", r, "_corr.RData", sep=""))
}
