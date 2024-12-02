# Simulation:
# Performance of fUC for parameter estimation
#   - High and normal signal to noise ratios
#   - AR plus non-AR
#   - compare d with other nonparametric estimators
library(devtools)
install_github("https://github.com/tobiashartl/CFFpack")
gc()
rm(list = ls())
library(fUCpack)
library(dplyr)
library(parallel)
library(CFFpack)
# Wd, etc
setwd("~/R/")
source("./help functions/fUC_arma_approx.R")

#setwd("/Users/tobias/Dokumente/Projekte/filtering unknown persistence/R/")
#source("./help functions/help_functions.R")
source("./help functions/KF.R")

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

# get nu
nuc <- sapply(1:nrow(setups), function(x) nucalc(setups[x, 1], setups[x, 2], setups[x, 3], ar))

#

for ( i in 1:nrow(setups)){
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    d <- as.numeric(setup[2])
    nu <- nucalc(n, d, setups[i, 3], ar)
    r <- setups[i, 3]
    cat("Iteration ", i, "\n")
    if(file.exists(file = paste("./MC/MC_1/ARMA/Sim_R", R, "_n", n, "_d", d, "_r", r, "_corr0.RData", sep=""))) next
    set.seed(42)
    
    
    # load the ARMA coefficients
    data(list = paste("d2arma",3,3,"_n",n, sep =""))
    
    # Generate the data
    x <- frac_diff_multi(matrix(rnorm(n*R, mean=0, sd=1), n, R), -d)
    u <- matrix(rnorm(n*R, mean=0, sd=sqrt(nu)), n, R)
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- yy<-x + c
    
    
    # grid of starting values
    d.grid <- seq(0.75, 1.75, 0.25)
    nu_grid <- c(0.01, 1, 100, 1000) %>% log()
    ar.1.grid <- c(-1.8,  -1.6, -1.4)
    ar2.grid  <- c(0.6,  0.8,  1)
    
    
    gr.start <- expand.grid(d.grid, nu_grid, ar.1.grid, ar2.grid) %>% as.matrix()
    gr.start <- gr.start[sapply(1:nrow(gr.start), function(i) toComp(-c(gr.start[i, 3:4]))$stable) ,]
    
    
    
    optfn <- function(n, y, x){
        tryCatch({
            # Estimate
            # check grid
            grid.st <- cbind(gr.start, 
                             sapply(1:nrow(gr.start), 
                                    function(i) fUC_opt_ML_ARMA_approx(gr.start[i, ], nu.opt=TRUE,
                                                                       y=y, START = 2, corr =FALSE,
                                                                       pq=c(2, 0),penalty.corr=FALSE,
                                                                       nulim = nulim,  
                                                                       deterministics = FALSE)))
            
            par0 <- grid.st[which.min(grid.st[,5]),-5]
            
            (est <- optim(par=par0, control = list(), 
                          fn = fUC_opt_ML_ARMA_approx, nu.opt=TRUE,
                          y=y, START = 2, corr =FALSE,
                          pq=c(2, 0),penalty.corr=FALSE,
                          nulim = nulim, d.int = c(0.5, 2.5),
                          deterministics = FALSE))
            
            # check for corner solutions: 
            j=1
            while(est$par[1] < .55 | est$par[1] > 2.45 | exp(est$par[2]) < .001 | exp(est$par[2]) > 250000){
                if(j > 10) break
                j=j+1
                par0 <- grid.st[order(grid.st[,5]),-5][j , ]
                (est <- optim(par=par0, 
                              fn = fUC_opt_ML_ARMA_approx, nu.opt=TRUE,
                              y=y, START = 2, corr =FALSE,
                              pq=c(2, 0),penalty.corr=FALSE,
                              nulim = nulim, d.int = c(0.5, 2.5),
                              deterministics = FALSE))
                
            }
            
            if(is.na(est$par[1])  |  est$par[1] < .55 | est$par[1] > 2.45){
                # primitive starting value: set 1
                par0 <- c(2, grid.st[which.min(grid.st[,5]),2:4])
                est <- tryCatch({
                    (optim(par=par0, 
                           fn = fUC_opt_ML_ARMA_approx, nu.opt=TRUE,
                           y=y, START = 2, corr =FALSE, d.int = c(0.5, 2.5),
                           pq=c(2, 0),penalty.corr=FALSE,
                           nulim = nulim,
                           deterministics = FALSE))
                }, error = function(e) return(NA)
                )
                if(is.na(est$par[1]) |  est$par[1] < .55 | est$par[1] > 2.45){
                    # primitive starting value: set 2
                    par0 <- c(1, grid.st[which.min(grid.st[,5]),2:4])
                    est <- tryCatch({
                        (optim(par=par0, 
                               fn = fUC_opt_ML_ARMA_approx, nu.opt=TRUE,
                               y=y, START = 2, corr =FALSE,
                               pq=c(2, 0),penalty.corr=FALSE,
                               nulim = nulim, d.int = c(0.5, 2.5),
                               deterministics = FALSE))
                    }, error = function(e) return(NA)
                    )
                }
                if(is.na(est[1])) est <- list(par = NA)
            }
            
            
            
            
            par <- est$par
            nu  <- exp(par[2])
            ar <- est$par[-(1:2)]
            
            KS <- fUC_KS_ARMA_approx(est$par, d.int = c(0.5, 2.5),
                                     y=y, START = 1, corr =FALSE,
                                     pq=c(2, 0),penalty.corr=FALSE,
                                     nulim = nulim,
                                     deterministics = FALSE)
            
            Rsq <- summary(lm(x ~ KS$x))$r.squared
            
            
            # Return results
            results <- c(par[1], exp(par[2]), ar, Rsq)
            names(results) <- c("d", "nu", "ar_1", "ar_2", "Rsq")
            return(results)
        }, error = function(e) return(rep(NA, 5))
        )
        
    }
    
    #optfn(n, y[,1], x[,1])
    ### fUC part
    cl <- makeCluster(32)
    clusterExport(cl, c("optfn", "fUC_opt_ML_ARMA_approx", "y", "fUC_KS_ARMA_approx", "ma_inf", "nulim",
                        "embed0", "fUC_KF_approx", "fUC_KS_approx", "x", "frac_diff", "lm", "n", "gr.start"))
    clusterEvalQ(cl, library(fUCpack))
    clusterEvalQ(cl, library(CFFpack))
    
    RESULTS <- parSapply(cl, 1:R, function(j) optfn(n, y[, j], x[, j]))
    stopCluster(cl)
    
    
    RESULTS <- rbind(RESULTS, c(dEW_45), c(dEW_50), c(dEW_55), c(dEW_60), c(dEW_65), c(dEW_70))
    rownames(RESULTS) <- c("d", "nu", "ar_1", "ar_2", "Rsq", 
                           "d_45", "d_50", "d_55", "d_60", "d_65", "d_70")
    save(RESULTS, file = paste("./MC/MC_1/ARMA/Sim_R", R, "_n", n, "_d", d, "_r", r, "_corr0.RData", sep=""))
}

