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

# load data
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
est.final <- readRDS(file = paste("./app2/mRVOL_corr0_p", p, ".RDS", sep=""))


# =============================================================================
# Parameter analysis
# =============================================================================
theta <- est.final$par
(d <- theta[1])
(Q <- diag(exp(theta[2:3])))
(ar   <- est.final$par[-(1:3)]) 
toComp(-ar)


H <- est.final$Hessian
J <- numDeriv::jacobian(function(x) (exp(x)), est.final$par[2:3])
cov0 <- solve(H)
Cov <- bdiag(1, J[,], diag(p))%*%cov0%*%t(bdiag(1, J[, ], diag(p)))
theta.trans <- c(est.final$par[1], exp(est.final$par[2:3]), est.final$par[-(1:3)])
se.trans <- sqrt(abs(diag(Cov)))
EST <- cbind(theta.trans, se.trans) 
colnames(EST) <- c("par", "se") 
rownames(EST) <- c("d", "Q11", "Q22", paste("ar_", 1:p, sep=""))
EST
EST[,1] / EST[,2]
EST[1,1] + c(-qt(.975, length(y))*EST[1,2], qt(.975, length(y))*EST[1,2])
(EST[4,1]/EST[2,1])


# add further relevant terms
k <- length(est.final$par)
(BIC          <- k * log(length(y)) + 2*est.final$value)
(AIC          <- 2 * k              + 2*est.final$value)

SSR <- (fUC_comp((y), d, Q, ar, corr=TRUE)$v[-1])^2 %>% sum(.)
add <- c(Q[2,2]/Q[1,1], Q[2,1]/Q[1,1], cov2cor(Q)[2,1], est.final$value, SSR, AIC, BIC)
EST <- rbind(EST, cbind(add, NA))
rownames(EST) <- c("d", "Q11", "Q22", paste("ar_", 1:p, sep=""),
                   "nu_1", "nu_2", "rho", "ll", "CSS", "AIC", "BIC")
colnames(EST)<- c("par", "se")
saveRDS(EST, file = "./app2/fUC_ML_Results.RDS")



# =============================================================================
# Trend-cycle decomposition
# =============================================================================

det.trend <- 0
plot(y, type="l")
lines(det.trend, col="2")


TC      <- fUC_smooth(y-det.trend, d, Q, ar, corr=TRUE)
plot(y, type="l")
lines(det.trend + TC$x, col="2")






# nino periods: 1950 onwards: https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php
# before: 





data.name <- "RVOLm"

getSymbols("^GSPC", from = '1927-12-30',
           to = "2024-10-31",warnings = FALSE,
           auto.assign = TRUE)

# Wd, etc
setwd("/Users/tobias/Dokumente/Projekte/Filtering unknown persistence/R/code")
source("./help functions/fUC_arma_approx.R")

# Calculate monthly volatility:
DATE = (index(GSPC))
data <- GSPC %>%
    as.data.frame %>%
    mutate(DATE = DATE) %>%
    mutate(DATE.m = as.yearmon(DATE))


# calculate daily returns

data <- data %>%
    #group_by(DATE.m) %>%
    mutate(return.open = c(NA, diff(log(GSPC.Open))),
           return.close = c(NA, diff(log(GSPC.Close))),
           rsq.open = c(NA, diff(log(GSPC.Open))^2),
           rsq.close = c(NA, diff(log(GSPC.Close))^2),
           crash = return.close < -0.05
    )%>%
    filter(DATE >= as.Date("1960-01-01")) %>% 
    #MVOL.open.adj = ((GSPC.Open - mean(GSPC.Open))/mean(GSPC.Open))^2,
    #     MVOL.close.adj = ((GSPC.Close - mean(GSPC.Close))/mean(GSPC.Close))^2,
    #     MVOL.open = (GSPC.Open - mean(GSPC.Open))^2,
    #     MVOL.close = (GSPC.Close - mean(GSPC.Close))^2,
    #     MVOL.OvC = (GSPC.Open - GSPC.Close)^2,
    #     MVOL.HvL = (GSPC.High - GSPC.Low)^2) %>%
    group_by(DATE.m) %>%
    summarize(RVOL.open_m = mean(rsq.open),
              RVOL.close_m = mean(rsq.close),
              return_m = sum(return.close),
              crash_col = ifelse(any(crash), "grey70", "white"),
              crash_col2 = ifelse(any(return_m < -0.10), "grey70", "white"),
              crash_alph = ifelse(any(crash), 0.15, 0),
              crash_alph2 = ifelse(any(return_m < -0.10), 0.15, 0),
              crash_col_c = ifelse(any(crash) | any(return_m < -0.10), "grey70", "white"),
              crash_alp_c = ifelse(any(crash) | any(return_m < -0.10), 0.15, 0)
              
              # MVOL.open_m.adj = mean(MVOL.open.adj),
              #  MVOL.close_m.adj = mean(MVOL.close.adj),
              # MVOL.hvl_m = mean(MVOL.HvL),
              #  MVOL.ovc_m = mean(MVOL.OvC)
    ) %>%
    mutate(lRVOL.open_m = log(RVOL.open_m) * 100,
           lRVOL.close_m = log(RVOL.close_m) * 100)



data.plot <- data.frame(
    y=y,
    trend = TC$x + det.trend, 
    x = TC$x, 
    cycle = TC$c,
    crash = data$crash_col,
    crash_alp = data$crash_alph,
    crash2 = data$crash_col2,
    crash_alp2 = data$crash_alph2,
    crash_c = data$crash_col_c,
    crash_alp_c = data$crash_alp_c,
    
    crash_min = seq(from = as.Date("1960-01-01"), to = as.Date("2024-10-01"),  by = "month"),
    crash_max = seq(from = as.Date("1960-01-31"), to = as.Date("2024-10-31"),  by = "month"),
    time = seq(from = as.Date("1960-01-01"), to = as.Date("2024-10-01"),  by = "month"))


saveRDS(data.plot, file = "./app2/fUC_ML_TC.RDS")

cbbPalette <- wesanderson::wes_palette("Darjeeling1", n=5)

library(ggplot2)
gg1 <- ggplot(data.plot, aes(x = time, y = y)) + 
    # la nina
    geom_rect(xmin=data.plot$crash_min, 
              xmax=data.plot$crash_max, 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=data.plot$crash_alp_c*0.15, col=data.plot$crash_c) +
    labs(x = "time", y="")+
    geom_line(color="black") + 
    ggtitle(bquote("Long-run realized volatility")) + 
    geom_line(mapping = aes(x = time, y = trend), color = cbbPalette[1], linetype = "twodash") +
    #    geom_line(mapping = aes(x = time, y = tau), color = cbbPalette[1], linetype = "twodash") +
    theme_classic()
gg1
gg12 <- ggplot(data.plot, aes(x = time, y = y)) + 
    # la nina
    geom_rect(xmin=data.plot$crash_min, 
              xmax=data.plot$crash_max, 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=data.plot$crash_alp*0.15, col=data.plot$crash) +
    labs(x = "time", y="")+
    geom_line(color="black") + 
    ggtitle(bquote("Long-run realized volatility")) + 
    geom_line(mapping = aes(x = time, y = trend), color = cbbPalette[1], linetype = "twodash") +
    #    geom_line(mapping = aes(x = time, y = tau), color = cbbPalette[1], linetype = "twodash") +
    theme_classic()


gg2 <- ggplot(data.plot, aes(x = time, y = cycle)) +
    # la nina
    geom_rect(xmin=data.plot$crash_min, 
              xmax=data.plot$crash_max, 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=data.plot$crash_alp_c*0.15, col=data.plot$crash_c) +
    labs(x = "time", y="")+
    geom_line(color="black") + 
    ggtitle(bquote("Short-run realized volatility")) + 
    #    geom_line(mapping = aes(x = time, y = tau), color = cbbPalette[1], linetype = "twodash") +
    theme_classic()
gg2



