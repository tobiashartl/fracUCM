gc()
rm(list = ls())
library(fUCpack)
library(dplyr)

# Wd, etc
setwd("~/R/")
source("./help functions/KF.R")
# Settings
n <- c(100, 200, 300)
d <- c(0.75, 1, 1.75)
ratio <- c(1, 10, 30, 0)
R <- 1000
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


results.list <- list()
X.res <- list()
C.res <- list()
T.res <- list()
X.true <- list()
C.true <- list()
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    d <- as.numeric(setup[2])
    nu <- nucalc(n, d, setups[i, 3], ar)
    r <- setups[i, 3]
    # check if simulation has already been done
    set.seed(42)
    
    # Generate the data
    x <- frac_diff_multi(matrix(rnorm(n*R, mean=0, sd=1), n, R), -d)
    u <- matrix(rnorm(n*R, mean=0, sd=sqrt(nu)), n, R)
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- yy<-x + c
    
    trend <- (1:n) 
    y <- y + trend
    
    X.mat <- x
    C.mat <- c
    X.true[[i]] <- X.mat
    C.true[[i]] <- C.mat
    try({
        load(paste("./MC/MC_6/CSS/Sim_R1000_n", true.par[1], "_d", true.par[2], "_r", true.par[3],
                   "_corr0_trend.RData", sep=""))
        results.list[[i]] <- RESULTS
        
        filterfun <- function(theta, y){
            tryCatch({
                det      <- fUC_opt(c(theta[1], log(theta[2]), theta[-(1:2)]), y, 
                                    nulim = c(0, Inf), d.int = c(0, 3),
                                    corr=FALSE, deterministics = "trend", START = 1, 
                                    return.det = TRUE)
                KS <- fUC_smooth(y - det*(1:n), theta[1], theta[2], theta[-(1:2)], corr=FALSE)
                return(c(KS$x, KS$c))
            }, error = function(e) return(rep(NA, 3*length(y))))
        }
        
        cl <- makeCluster(6)
        clusterExport(cl, c("fUC_opt", "y", "fUC_comp", "ma_inf", "true.par", "RESULTS",
                            "embed0", "fUC_smooth", "x", "frac_diff", "lm", "n", "filterfun"))
        clusterEvalQ(cl, library(fUCpack))
        XC <- parSapply(cl, 1:R, function(j) filterfun(RESULTS[1:4,j], y[,j]))
        stopCluster(cl)
        X.mat <- XC[1:n, ]
        C.mat <- XC[(n+1):(2*n), ]
    })
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}



save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/MC_6/CSS/results.RData")


# ML
# Load esimates
R <- 1000
ar <- c(-1.6, 0.8)
results.list <- list()
X.res <- list()
C.res <- list()
X.true <- list()
C.true <- list()
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    d <- as.numeric(setup[2])
    nu <- nucalc(n, d, setups[i, 3], ar)
    r <- setups[i, 3]
    # check if simulation has already been done
    set.seed(42)
    
    # Generate the data
    x <- frac_diff_multi(matrix(rnorm(n*R, mean=0, sd=1), n, R), -d)
    u <- matrix(rnorm(n*R, mean=0, sd=sqrt(nu)), n, R)
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- yy<-x + c
    
    
    trend <-  (1:n) 
    
    y <- y + trend
    
    X.mat <- x
    C.mat <- c
    X.true[[i]] <- X.mat
    C.true[[i]] <- C.mat
    try({
        load(paste("./MC/MC_6/ML/Sim_R1000_n", true.par[1], "_d", true.par[2], "_r", true.par[3],
                   "_corr0_trend.RData", sep=""))
        results.list[[i]] <- RESULTS
        
        
        filterfun <- function(theta, y){
            tryCatch({
                det      <- fUC_opt_ML(c(theta[1], log(theta[2]), theta[-(1:2)]), y, 
                                       nulim = c(0, Inf), nu.opt = TRUE, d.int = c(0, 3),
                                       corr=FALSE, deterministics = "trend", START = 1, 
                                       return.det = TRUE)
                KS <- fUC_smooth(y - det*(1:n), theta[1], theta[2], theta[-(1:2)], corr=FALSE)
                return(c(KS$x, KS$c))
            }, error = function(e) return(rep(NA, 2*length(y))))
        }
        filterfun(RESULTS[1:4,1], y[,1])[1:100] %>% plot
        lines(-x[,1], type="l")
        
        cl <- makeCluster(6)
        clusterExport(cl, c("fUC_opt_ML", "y", "fUC_comp", "ma_inf", "true.par", "RESULTS",
                            "embed0", "fUC_smooth", "x", "frac_diff", "lm", "n", "filterfun"))
        clusterEvalQ(cl, library(fUCpack))
        XC <- parSapply(cl, 1:R, function(j) filterfun(RESULTS[1:4,j], y[,j]))
        stopCluster(cl)
        X.mat <- XC[1:n, ]
        C.mat <- XC[-(1:n), ]
    })
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}


save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/MC_6/ML/results_ML.RData")


# I(1) CSS
# Load esimates
R <- 1000
ar <- c(-1.6, 0.8)
results.list <- list()
X.res <- list()
C.res <- list()
X.true <- list()
C.true <- list()
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    d <- as.numeric(setup[2])
    nu <- nucalc(n, d, setups[i, 3], ar)
    r <- setups[i, 3]
    # check if simulation has already been done
    set.seed(42)
    
    # Generate the data
    x <- frac_diff_multi(matrix(rnorm(n*R, mean=0, sd=1), n, R), -d)
    u <- matrix(rnorm(n*R, mean=0, sd=sqrt(nu)), n, R)
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- yy<-x + c
    
    trend <- (1:n) 
    
    y <- y + trend
    
    X.mat <- x
    C.mat <- c
    X.true[[i]] <- X.mat
    C.true[[i]] <- C.mat
    try({
        load(paste("./MC/MC_6/CSS_INTEGER/Sim_R1000_n", true.par[1], "_d", true.par[2], "_r", true.par[3],
                   "_trend.RData", sep=""))
        results.list[[i]] <- RESULTS
        
        filterfun <- function(theta, y){
            tryCatch({
                det      <- fUC_opt(c(1, log(theta[1]), theta[-1]), y, nulim = c(0, Inf), d.int = c(0, 3),
                                    corr=FALSE, deterministics = "trend", START = 1, 
                                    return.det = TRUE)
                
                KS <- fUC_smooth(y-det*(1:n), 1, theta[1], theta[-(1:1)], corr=FALSE)
                return(c(KS$x, KS$c))
            }, error = function(e) return(rep(NA, 2*length(y))))
        }
        
        cl <- makeCluster(6)
        clusterExport(cl, c("fUC_opt", "y", "fUC_comp", "ma_inf", "true.par", "RESULTS",
                            "embed0", "fUC_smooth", "x", "frac_diff", "lm", "n", "filterfun"))
        clusterEvalQ(cl, library(fUCpack))
        XC <- parSapply(cl, 1:R, function(j) filterfun(RESULTS[1:3,j], y[,j]))
        stopCluster(cl)
        X.mat <- XC[1:n, ]
        C.mat <- XC[-(1:n), ]
    })
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}

save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/MC_6/CSS_INTEGER/results_i1.RData")


# I(1) ML
# Load esimates
R <- 1000
ar <- c(-1.6, 0.8)
results.list <- list()
X.res <- list()
C.res <- list()
X.true <- list()
C.true <- list()
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    d <- as.numeric(setup[2])
    nu <- nucalc(n, d, setups[i, 3], ar)
    r <- setups[i, 3]
    # check if simulation has already been done
    set.seed(42)
    
    # Generate the data
    x <- frac_diff_multi(matrix(rnorm(n*R, mean=0, sd=1), n, R), -d)
    u <- matrix(rnorm(n*R, mean=0, sd=sqrt(nu)), n, R)
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- yy<-x + c
    
    trend <- (1/1) * (1:n) 
    y <- y + trend
    
    X.mat <- x
    C.mat <- c
    X.true[[i]] <- X.mat
    C.true[[i]] <- C.mat
    try({
        load(paste("./MC/MC_6/ML_INTEGER/Sim_R1000_n", true.par[1], "_d", true.par[2], "_r", true.par[3],
                   "_trend.RData", sep=""))
        results.list[[i]] <- RESULTS
        
        filterfun <- function(theta, y){
            tryCatch({
                det      <- fUC_opt_ML(c(1, log(theta[1]), theta[-1]), y, 
                                       nulim = c(0, Inf), nu.opt = TRUE, d.int = c(0, 3),
                                       corr=FALSE, deterministics = "trend", START = 1, 
                                       return.det = TRUE)
                
                KS <- fUC_smooth(y - det*(1:n), theta[1], theta[2], theta[-(1:2)], corr=FALSE)
                return(c(KS$x, KS$c))
            }, error = function(e) return(rep(NA, 2*length(y))))
        }
        
        cl <- makeCluster(6)
        clusterExport(cl, c("fUC_opt_ML", "y", "fUC_comp", "ma_inf", "true.par", "RESULTS",
                            "embed0", "fUC_smooth", "x", "frac_diff", "lm", "n", "filterfun"))
        clusterEvalQ(cl, library(fUCpack))
        XC <- parSapply(cl, 1:R, function(j) filterfun(RESULTS[1:3,j], y[,j]))
        stopCluster(cl)
        X.mat <- XC[1:n, ]
        C.mat <- XC[-(1:n), ]
    })
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}

save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/MC_6/ML_Integer/results_i1_ML.RData")





# ARMA
source("./help functions/fUC_arma_approx.R")
library(CFFpack)
# Load esimates
R <- 1000
ar <- c(-1.6, 0.8)
results.list <- list()
X.res <- list()
C.res <- list()
X.true <- list()
C.true <- list()
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    d <- as.numeric(setup[2])
    nu <- nucalc(n, d, setups[i, 3], ar)
    r <- setups[i, 3]
    # check if simulation has already been done
    set.seed(42)
    
    # Generate the data
    x <- frac_diff_multi(matrix(rnorm(n*R, mean=0, sd=1), n, R), -d)
    u <- matrix(rnorm(n*R, mean=0, sd=sqrt(nu)), n, R)
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- yy<-x + c
    
    trend <- (1:n) 
    
    y <- y + trend
    
    X.mat <- x
    C.mat <- c
    X.true[[i]] <- X.mat
    C.true[[i]] <- C.mat
    try({
        load(paste("./MC/MC_6/ARMA/Sim_R1000_n", true.par[1], "_d", true.par[2], "_r", true.par[3],
                   "_corr0_trend.RData", sep=""))
        results.list[[i]] <- RESULTS[, colSums(is.na(RESULTS))==0 & RESULTS[1, ] < 1e+06]
        
        
        filterfun <- function(theta, y){
            tryCatch({
                trend <- fUC_opt_ML_ARMA_approx(theta,
                                                y=y, START = 2, corr =FALSE, d.int = c(0, 3),
                                                pq=c(2, 0),penalty.corr=FALSE,
                                                nulim = c(1/1000, 1e+07),
                                                deterministics = "trend", return.det = TRUE)
                KS <- fUC_KS_ARMA_approx(theta,
                                         y=y - trend*(1:n), START = 2, corr =FALSE, d.int = c(0, 3),
                                         pq=c(2, 0),penalty.corr=FALSE,
                                         nulim = c(1/1000, 1e+07))
                
                #KS <- fUC_smooth(y, theta[1], theta[2], theta[-(1:2)], corr=FALSE)
                return(c(KS$x, KS$c))
            }, error = function(e) return(rep(NA, 2*length(y))))
        }
        cl <- makeCluster(cl)
        clusterExport(cl, c("fUC_opt_ML_ARMA_approx", "fUC_KS_ARMA_approx", "y", "fUC_comp", "ma_inf", "true.par", "RESULTS",
                            "embed0", "fUC_smooth", "x", "frac_diff", "lm", "n", "filterfun"))
        clusterEvalQ(cl, library(fUCpack))
        clusterEvalQ(cl, library(CFFpack))
        XC <- parSapply(cl, 1:R, function(j) filterfun(c(RESULTS[1,j], log(RESULTS[2,j]), RESULTS[3:4, j]), y[,j]))
        stopCluster(cl)
        X.mat <- XC[1:n, ]
        C.mat <- XC[-(1:n), ]
    })
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}

save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/MC_6/ARMA/results_ARMA.RData")






