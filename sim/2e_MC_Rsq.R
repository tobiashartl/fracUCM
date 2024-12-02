gc()
rm(list = ls())
library(fUCpack)
library(dplyr)

# Wd, etc
setwd("/Users/tobias/Dokumente/Projekte/filtering unknown persistence/R/code")
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
X.true <- list()
C.true <- list()
for(i in 1:nrow(setups)){
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
    
    try({
        load(paste("./MC/MC_2/CSS/Sim_R1000_n", n, "_d", d, "_r", r,
                   "_corr.RData", sep=""))
        results.list[[i]] <- RESULTS
        
        filterfun <- function(theta, y){
            nu  <- matrix(c(1, theta[3], theta[3], theta[2]), 2, 2)
            KS <- fUC_smooth(y, theta[1], nu, theta[-(1:3)], corr=TRUE)
            return(c(KS$x, KS$c))
        }
        
        XC <- sapply(1:R, function(j) filterfun(RESULTS[1:5,j], y[,j]))
        X.mat <- XC[1:n, ]
        C.mat <- XC[-(1:n), ]
    })
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}

save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/MC_2/CSS/results.RData")


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
    
    set.seed(42)
    
    rho <- -0.8
    Q <- matrix(c(1, rho*sqrt(1)*sqrt(nu), rho*sqrt(1)*sqrt(nu), nu), 2, 2)
    
    # Generate the data
    U <- sapply(1:R, function(x) mvtnorm::rmvnorm(n, mean=c(0, 0), sigma = Q))
    x <- frac_diff_multi(U[1:n,], d=-d)
    u <- U[-(1:n),]
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- x + c
    try({
        load(paste("./MC/MC_2/ML/Sim_R1000_n", n, "_d", d, "_r", r,
                   "_corr.RData", sep=""))
        results.list[[i]] <- RESULTS[, colSums(is.na(RESULTS))==0 & RESULTS[1, ] < 1e+06]
        
        
        filterfun <- function(theta, y){
            nu  <- mlogvech2mat(theta[2:4])
            KS <- fUC_smooth(y, theta[1], nu, theta[-(1:4)], corr=TRUE)
            return(c(KS$x, KS$c))
        }
        
        XC <- sapply(1:R, function(j) filterfun(RESULTS[2:7,j], y[,j]))
        X.mat <- XC[1:n, ]
        C.mat <- XC[-(1:n), ]
    })
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}

save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/MC_2/ML/results_ML.RData")


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
    
    set.seed(42)
    rho <- -0.8
    Q <- matrix(c(1, rho*sqrt(1)*sqrt(nu), rho*sqrt(1)*sqrt(nu), nu), 2, 2)
    
    # Generate the data
    U <- sapply(1:R, function(x) mvtnorm::rmvnorm(n, mean=c(0, 0), sigma = Q))
    x <- frac_diff_multi(U[1:n,], d=-d)
    u <- U[-(1:n),]
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- x + c
    try({
        load(paste("./MC/MC_2/CSS_INTEGER/Sim_R1000_n", n, "_d", d, "_r", r,
                   "_corr.RData", sep=""))
        results.list[[i]] <- RESULTS
        
        
        filterfun <- function(theta, y){
            Q <- matrix(c(1, theta[2], theta[2], theta[1]), 2, 2)
            KS <- fUC_smooth(y, 1, Q, -theta[-(1:2)], corr=TRUE)
            return(c(KS$x, KS$c))
        }
        
        XC <- sapply(1:R, function(j) filterfun(RESULTS[1:4,j], y[,j]))
        X.mat <- XC[1:n, ]
        C.mat <- XC[-(1:n), ]
    })
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}

save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/MC_2/CSS_INTEGER/results_i1.RData")


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
    
    set.seed(42)
    
    rho <- -0.8
    Q <- matrix(c(1, rho*sqrt(1)*sqrt(nu), rho*sqrt(1)*sqrt(nu), nu), 2, 2)
    
    # Generate the data
    U <- sapply(1:R, function(x) mvtnorm::rmvnorm(n, mean=c(0, 0), sigma = Q))
    x <- frac_diff_multi(U[1:n,], d=-d)
    u <- U[-(1:n),]
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- x + c
    try({
        load(paste("./MC/MC_2/ML_INTEGER/Sim_R1000_n", n, "_d", d, "_r", r,
                   "_corr.RData", sep=""))
        results.list[[i]] <- RESULTS[, colSums(is.na(RESULTS))==0]
        
        
        filterfun <- function(theta, y){
            nu  <- matrix(c(theta[1], theta[3], theta[3], theta[2]), 2, 2)
            KS <- fUC_smooth(y, 1, nu, theta[-(1:3)], corr=TRUE)
            
            
            
            return(c(KS$x, KS$c))
        }
        
        XC <- sapply(1:R, function(j) tryCatch(filterfun(RESULTS[1:5,j], y[,j]),
                                               error = function(e) return(rep(NA, 2*n))
        )
        )
        
        X.mat <- XC[1:n, ]
        C.mat <- XC[-(1:n), ]
    })
    
    
    calcRSQ <- function(x.true, x.est){
        Rsq <- mean(sapply(1:ncol(x.true), function(j) tryCatch(summary(lm(x.true[, j] ~ x.est[,j]))$r.squared,
                                                                error = function(e) return(NA))), na.rm = TRUE)
    }
    calcRSQ(x, X.mat)-> rsq
    
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}

save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/MC_2/ML_INTEGER/results_i1_ML.RData")








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
    
    set.seed(42)
    
    rho <- -0.8
    Q <- matrix(c(1, rho*sqrt(1)*sqrt(nu), rho*sqrt(1)*sqrt(nu), nu), 2, 2)
    
    # Generate the data
    U <- sapply(1:R, function(x) mvtnorm::rmvnorm(n, mean=c(0, 0), sigma = Q))
    x <- frac_diff_multi(U[1:n,], d=-d)
    u <- U[-(1:n),]
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- x + c
    try({
        load(paste("./MC/MC_2/ARMA/Sim_R1000_n", true.par[1], "_d", true.par[2], "_r", true.par[3],
                   "_corr.RData", sep=""))
        results.list[[i]] <- RESULTS[, colSums(is.na(RESULTS))==0]
        
        
        
        
        filterfun <- function(theta, y){
            tryCatch({
                KS <- fUC_KS_ARMA_approx(theta, d.int = c(0, 3),
                                         y=y, START = 2, corr =TRUE,
                                         pq=c(2, 0),penalty.corr=FALSE,
                                         nulim = c(1/100, 1e+07),
                                         deterministics = FALSE)
                
                #KS <- fUC_smooth(y, theta[1], theta[2], theta[-(1:2)], corr=FALSE)
                return(c(KS$x, KS$c))
            }, error = function(e) return(rep(NA, 2*length(y))))
        }
        XC <- sapply(1:R, function(j) filterfun(c(RESULTS[2:7,j]), y[,j]))
        X.mat <- XC[1:n, ]
        C.mat <- XC[-(1:n), ]
    })
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}

save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/MC_2/ARMA/results_ARMA.RData")









