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
ratio <- c(0, 30, 10, 1)
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


# Load CSS esimates
results.list <- list()
statistics <- statisticsb <- matrix(NA, nrow(setups), 12)
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    try({
        load(paste("./MC/MC_5/CSS/Sim_R1000_n", true.par[1], "_d", true.par[2], "_r", true.par[3],
                   "_outlier.RData", sep=""))
        results.list[[i]] <- RESULTS
        ### get the true nu
        nu <- nucalc(true.par$n, true.par$d, true.par$ratio, ar)
        
        MSE_d <- rowMeans((results.list[[i]][c("d", "d_45", "d_50", "d_55", "d_60", "d_65", "d_70"), ] - true.par$d)^2) %>%
            sqrt(.)
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8))^2)%>%
            sqrt(.)
        #SSR   <- mean(results.list[[i]][c("SSR"), ])
        Rsq   <- mean(results.list[[i]][c("Rsq"), ])
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - nu)^2) %>%
            sqrt(.)
        statistics[i, ] <- c(nu, MSE_d, MSE_nu, Rsq, MSE_ar)
        
        MSE_d <- apply((results.list[[i]][c("d", "d_45", "d_50", "d_55", "d_60", "d_65", "d_70"), ] - true.par$d), 1, mean) 
        MSE_ar <- -apply((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8)), 1, mean)
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - nu)) 
        statisticsb[i, ] <- c(nu, MSE_d, MSE_nu, Rsq, MSE_ar)
        
    })
    
}

rownames(statistics) <- rownames(statisticsb) <- apply(setups[,], 1, function(x) paste(x, collapse="_"))
colnames(statistics) <- colnames(statisticsb) <- c("nu0", "d",
                                                   "d_45", "d_50", "d_55", "d_60", "d_65", "d_70",
                                                   "nu", "Rsq", "ar1", "ar2" )


# Load ML esimates
results.list <- list()
statistics2 <- statistics2b <- matrix(NA, nrow(setups), 6)
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    try({
        load(paste("./MC/MC_5/ML/Sim_R1000_n", true.par[1], "_d", true.par[2], "_r", true.par[3],
                   "_outlier.RData", sep=""))
        results.list[[i]] <- RESULTS
        nu <- nucalc(true.par$n, true.par$d, true.par$ratio, ar)
        
        
        MSE_d <- mean((results.list[[i]][c("d"), ] - true.par$d)^2) %>%
            sqrt(.)
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8))^2)%>%
            sqrt(.)
        #SSR   <- mean(results.list[[i]][c("SSR"), ])
        Rsq   <- mean(results.list[[i]][c("Rsq"), ])
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - nu)^2) %>%
            sqrt(.)
        statistics2[i, ] <- c(nu, MSE_d, MSE_nu, Rsq, MSE_ar)
        
        
        MSE_d <- mean((results.list[[i]][c("d"), ] - true.par$d)) 
        MSE_ar <- -apply((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8)), 1, mean)
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - nu)) 
        statistics2b[i, ] <- c(nu, MSE_d, MSE_nu, Rsq, MSE_ar)
    })
    
}

rownames(statistics2) <- rownames(statistics2b) <- apply(setups[,], 1, function(x) paste(x, collapse="_"))
colnames(statistics2) <- colnames(statistics2b) <- c("nu0","d_ML",
                                                     "nu_ML", "Rsq_ML", "ar1_ML", "ar2_ML" )


# Load I(1) CSS estimates
results.list <- list()
statistics3 <- statistics3b <- matrix(NA, nrow(setups), 5)
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    try({
        load(paste("./MC/MC_5/CSS_INTEGER/Sim_R1000_n", true.par[1], "_d", true.par[2], "_r", true.par[3],
                   "_outlier.RData", sep=""))
        results.list[[i]] <- RESULTS
        nu <- nucalc(true.par$n, true.par$d, true.par$ratio, ar)
        
        SSR   <- mean(results.list[[i]][c("SSR"), ])
        Rsq   <- mean(results.list[[i]][c("Rsq"), ])
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8))^2)%>%
            sqrt(.)
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - nu)^2) %>% sqrt(.)
        statistics3[i, ] <- c(MSE_nu, SSR, Rsq, MSE_ar)
        
        MSE_ar <- apply((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8)), 1, mean)
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - nu)) 
        statistics3b[i, ] <- c(MSE_nu, SSR, Rsq, MSE_ar)
        
    })
    
}
rownames(statistics3) <- rownames(statistics3b) <- apply(setups[1:nrow(setups),], 1, function(x) paste(x, collapse="_"))
colnames(statistics3) <- colnames(statistics3b) <- c("nu_i1", "SSR_i1", "Rsq_i1", "ar1_i1", "ar2_i1")


# Load I(1) ML estimates
results.list <- list()
statistics4 <- statistics4b <- matrix(NA, nrow(setups), 5)
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    try({
        load(paste("./MC/MC_5/ML_INTEGER/Sim_R1000_n", true.par[1], "_d", true.par[2], "_r", true.par[3],
                   "_outlier.RData", sep=""))
        results.list[[i]] <- RESULTS
        nu <- nucalc(true.par$n, true.par$d, true.par$ratio, ar)
        SSR   <- mean(results.list[[i]][c("SSR"), ])
        Rsq   <- mean(results.list[[i]][c("Rsq"), ])
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8))^2)%>%
            sqrt(.)
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - nu)^2) %>% sqrt(.)
        statistics4[i, ] <- c(MSE_nu, SSR, Rsq, MSE_ar)
        
        MSE_ar <- apply((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8)), 1, mean)
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - nu)) 
        statistics4b[i, ] <- c(MSE_nu, SSR, Rsq, MSE_ar)
        
    })
    
}


rownames(statistics4) <- rownames(statistics4b) <- apply(setups[1:nrow(setups),], 1, function(x) paste(x, collapse="_"))
colnames(statistics4) <- colnames(statistics4b) <- c("nu_i1_ML", "SSR_i1_ML", "Rsq_i1_ML", "ar1_i1_ML", "ar2_i1_ML")




# Load ARMA esimates
results.list <- list()
statistics5 <- statistics5b <- matrix(NA, nrow(setups), 6)
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    try({
        load(paste("./MC/MC_5/ARMA/Sim_R1000_n", true.par[1], "_d", true.par[2], "_r", true.par[3],
                   "_outlier.RData", sep=""))
        results.list[[i]] <- RESULTS
        nu <- nucalc(true.par$n, true.par$d, true.par$ratio, ar)
        
        
        MSE_d <- mean((results.list[[i]][c("d"), ] - true.par$d)^2) %>%
            sqrt(.)
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8))^2)%>%
            sqrt(.)
        #SSR   <- mean(results.list[[i]][c("SSR"), ])
        Rsq   <- mean(results.list[[i]][c("Rsq"), ])
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - nu)^2) %>%
            sqrt(.)
        statistics5[i, ] <- c(nu, MSE_d, MSE_nu, Rsq, MSE_ar)
        
        
        MSE_d <- mean((results.list[[i]][c("d"), ] - true.par$d)) 
        MSE_ar <- -apply((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8)), 1, mean)
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - nu)) 
        statistics5b[i, ] <- c(nu, MSE_d, MSE_nu, Rsq, MSE_ar)
    })
    
}

rownames(statistics5) <- rownames(statistics5b) <- apply(setups[,], 1, function(x) paste(x, collapse="_"))
colnames(statistics5) <- colnames(statistics5b) <- c("nu0","d_ARMA",
                                                     "nu_ARMA", "Rsq_ARMA", "ar1_ARMA", "ar2_ARMA" )






statistics <- statistics %>% as.data.frame%>%tibble::rownames_to_column()
statistics2 <- statistics2 %>% as.data.frame%>%tibble::rownames_to_column()
statistics3 <- statistics3 %>% as.data.frame %>% tibble::rownames_to_column()
statistics4 <- statistics4 %>% as.data.frame %>% tibble::rownames_to_column()
statistics5 <- statistics5 %>% as.data.frame %>% tibble::rownames_to_column()


statistics$n <- setups$n
statistics$d0 <- setups$d
statistics$r   <- setups$r

### Table 1: RMSE d
tab1 <- left_join(statistics, statistics3) %>%
    left_join(statistics2)%>%
    left_join(statistics4)%>%
    left_join(statistics5)%>%
    dplyr::select("n",  "r", "d0", "nu0",
                  "d", "d_ML", "d_ARMA","d_50","d_60","d_70") %>%
    filter(n %in% c(100, 200, 300), # in the end: 100, 200, 300
           r %in% c(0, 30, 10, 1), 
           d0 %in% c(0.75, 1, 1.75))

tab1 <- tab1[with(tab1, order(n, r)),]
# tbd
tab1$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab1$r[!c(1:36)%in%seq(1, 36, by = 3)] <- NA
tab1$r[tab1$r==0.0 & !is.na(tab1$r)] <- NA

table <- tab1%>% 
    xtable::xtable(digits = c(0, 0, 0, 2, -2, rep(2, 6)), include.rownames=FALSE)
print(table, include.rownames = FALSE)



### Table 2: Bias d
statisticsb <- statisticsb %>% as.data.frame%>%tibble::rownames_to_column()
statistics2b <- statistics2b %>% as.data.frame%>%tibble::rownames_to_column()
statistics3b <- statistics3b %>% as.data.frame %>% tibble::rownames_to_column()
statistics4b <- statistics4b %>% as.data.frame %>% tibble::rownames_to_column()
statistics5b <- statistics5b %>% as.data.frame %>% tibble::rownames_to_column()


statisticsb$n <- setups$n
statisticsb$d0 <- setups$d
statisticsb$r   <- setups$r

tab2 <- left_join(statisticsb, statistics3b) %>%
    left_join(statistics2b)%>%
    left_join(statistics4b)%>%
    left_join(statistics5b)%>%
    
    dplyr::select("n",  "r", "d0",
                  "d", "d_ML", "d_ARMA", "d_50", "d_60", "d_70") %>%
    filter(n %in% c(100, 200, 300), # in the end: 100, 200, 300
           r %in% c(0, 30, 10, 1), 
           d0 %in% c(0.75, 1, 1.75))

# tbd
tab2 <- tab2[with(tab2, order(n, r)),]
tab2$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab2$r[!c(1:36)%in%seq(1, 36, by = 3)] <- NA
tabjoint <- cbind(tab1, tab2[, -(1:3)])
tabjoint$r[tabjoint$r==0]<-NA
saveRDS(tabjoint, file = "./MC/MC_5/tab1_biasrmse_d.RDS")


table <- tabjoint%>% 
    xtable::xtable(digits = c(0, 0, 1, 2, 2, rep(2, 12)), include.rownames=FALSE)
print(table, include.rownames = FALSE)



tab3 <- left_join(statistics, statistics3) %>%
    left_join(statistics2)%>%
    left_join(statistics4)%>%
    left_join(statistics5)%>%
    dplyr::select("n",  "r", "d0", "nu0",
                  "nu", "nu_ML", "nu_ARMA", "nu_i1","nu_i1_ML",
                  "ar1", "ar1_ML", "ar1_ARMA", "ar1_i1", "ar1_i1_ML",
                  "ar2", "ar2_ML", "ar2_ARMA", "ar2_i1", "ar2_i1_ML") %>%
    filter(n %in% c(100, 200, 300), # in the end: 100, 200, 300
           r %in% c(0, 30, 10, 1), 
           d0 %in% c(0.75, 1, 1.75))


#tbd 
tab3 <- tab3[with(tab3, order(n, r)),]
tab3$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab3$r[!c(1:36)%in%seq(1, 36, by = 3)] <- NA
tab3$r[tab3$r==0]<-NA
saveRDS(tab3, file = "./MC/MC_5/tab2_rmse_oth.RDS")



tab3[tab3>100 & !is.na(tab3>100)] <- round(tab3[tab3>100& !is.na(tab3>100)], digits = 0)

table <- tab3%>% 
    xtable::xtable(digits = c(0, 0, 1, 2, 2, c(2, 2, 2, 2, 2), rep(2, 10)), include.rownames=FALSE)
print(table, include.rownames = FALSE)


tab3 <- left_join(statisticsb, statistics3b) %>%
    left_join(statistics2b)%>%
    left_join(statistics4b)%>%
    left_join(statistics5b)%>%
    dplyr::select("n",  "r", "d0", "nu0",
                  "nu", "nu_ML", "nu_ARMA", "nu_i1","nu_i1_ML",
                  "ar1", "ar1_ML", "ar1_ARMA", "ar1_i1", "ar1_i1_ML",
                  "ar2", "ar2_ML", "ar2_ARMA", "ar2_i1", "ar2_i1_ML") %>%
    filter(n %in% c(100, 200, 300), # in the end: 100, 200, 300
           r %in% c(0, 30, 10, 1), 
           d0 %in% c(0.75, 1, 1.75))

#tbd 
tab3 <- tab3[with(tab3, order(n, r)),]
tab3$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab3$r[!c(1:36)%in%seq(1, 36, by = 3)] <- NA
tab3$r[tab3$r==0]<-NA
saveRDS(tab3, file = "./MC/MC_5/tab3_bias_oth.RDS")


tab3[tab3>100 & !is.na(tab3>100)] <- round(tab3[tab3>100& !is.na(tab3>100)], digits = 0)
tab3[tab3< -100 & !is.na(tab3< -100)] <- round(tab3[tab3< -100& !is.na(tab3< -100)], digits = 0)

table <- tab3%>% 
    xtable::xtable(digits = c(0, 0, 1, 2, 2, c(2, 2, 2, 2, 2), rep(2, 10)), include.rownames=FALSE)
print(table, include.rownames = FALSE)



### Check goodness of fit for the different models
# CSS
load(file = "./MC/MC_5/CSS/results.RData")
X.hat.CSS <- X.res
C.hat.CSS <- C.res
# ML
load(file = "./MC/MC_5/ML/results_ML.RData")
X.hat.ML <- X.res
C.hat.ML <- C.res
# CSS i1
load(file = "./MC/MC_5/CSS_INTEGER/results_i1.RData")
X.hat.i1 <- X.res
C.hat.i1 <- C.res
# CSS i1 ML
load(file = "./MC/MC_5/ML_INTEGER/results_i1_ML.RData")
X.hat.i1.ML <- X.res
C.hat.i1.ML <- C.res
# ARMA ML
load(file = "./MC/MC_5/ARMA/results_ARMA.RData")
X.hat.ARMA <- X.res
C.hat.ARMA <- C.res




calcRSQ <- function(x.true, x.est){
    Rsq <- mean(sapply(1:ncol(x.true), function(j) summary(lm(x.true[, j] ~ I(1:NROW(x.true)) + x.est[,j]))$r.squared))
}

Rsq.CSS <- sapply(1:nrow(setups), function(j) calcRSQ(X.true[[j]], X.hat.CSS[[j]]))
Rsq.ML <- sapply(1:nrow(setups), function(j) calcRSQ(X.true[[j]], X.hat.ML[[j]]))
Rsq.i1 <- sapply(1:nrow(setups), function(j) calcRSQ(X.true[[j]], X.hat.i1[[j]]))
Rsq.i1.ML <- sapply(1:nrow(setups), function(j) calcRSQ(X.true[[j]], X.hat.i1.ML[[j]]))
Rsq.ARMA <- sapply(1:nrow(setups), function(j) calcRSQ(X.true[[j]], X.hat.ARMA[[j]]))
Rsq.CSS.C <- sapply(1:nrow(setups), function(j) calcRSQ(C.true[[j]], C.hat.CSS[[j]]))
Rsq.ML.C <- sapply(1:nrow(setups), function(j) calcRSQ(C.true[[j]], C.hat.ML[[j]]))
Rsq.i1.C <- sapply(1:nrow(setups), function(j) calcRSQ(C.true[[j]], C.hat.i1[[j]]))
Rsq.i1.ML.C <- sapply(1:nrow(setups), function(j) calcRSQ(C.true[[j]], C.hat.i1.ML[[j]]))
Rsq.ARMA.C <- sapply(1:nrow(setups), function(j) calcRSQ(C.true[[j]], C.hat.ARMA[[j]]))

# calculate nu

#nuc
tab4 <- cbind(setups[, c("n", "ratio", "d")], nuc, Rsq.CSS, Rsq.ML, Rsq.ARMA, Rsq.i1, Rsq.i1.ML, 
              Rsq.CSS.C, Rsq.ML.C, Rsq.ARMA.C, Rsq.i1.C, Rsq.i1.ML.C)
#tbd
tab4 <- tab4[with(tab4, order(n, ratio)),]
tab4$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab4$ratio[!c(1:36)%in%seq(1, 36, by = 3)] <- NA
tab4$ratio[tab4$ratio==0]<-NA
saveRDS(tab4, file = "./MC/MC_5/tab4_rsq.RDS")

#tab4[tab4 > 100] <- round(tab4[tab4>100], digits = 0)
table <- tab4%>% 
    xtable::xtable(digits = c(0, 0, 1, 2, 2, rep(2, 10)), include.rownames=FALSE)
print(table, include.rownames = FALSE)





