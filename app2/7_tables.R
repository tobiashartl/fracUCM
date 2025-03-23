#Packages
gc()
rm(list = ls())
library(fUCpack)
library(dplyr)
library(parallel)
library(CFFpack)
library(xtable)


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





# load estimation results
Est.id <- readRDS(file = "./app2/fUC_ML_Results.RDS")
Est.i1 <- readRDS(file = "./app2/Benchmark_i1_ML_Results.RDS")
Est.i2 <- readRDS(file = "./app2/Benchmark_i2_ML_Results.RDS")

Est.id["d", 1] + c(qnorm(0.975)*Est.id["d", 2], -qnorm(0.975)*Est.id["d", 2])
Est.id["d", 1] + c(qnorm(0.995)*Est.id["d", 2], -qnorm(0.995)*Est.id["d", 2])

Est.id <- rbind(Est.id[1:4,], matrix(NA, nrow = 4,ncol=2), Est.id[-(1:4),])
Est.i1 <- rbind(Est.i1[1:4,], matrix(NA, nrow = 4,ncol=2), Est.i1[-(1:4),])

tab.all <- data.frame(est.id = Est.id[,1],
                      se.id = Est.id[,2],
                      est.i1 = Est.i1[,1],
                      se.i1 = Est.i1[,2],
                      est.i2 = Est.i2[,1],
                      se.i2 = Est.i2[,2]
)

tab.all["ll", ] <- - tab.all["ll", ]

colnames(tab.all) = c("Estimate", "Std. Error", "Estimate", "Std. Error", "Estimate", "Std. Error")
mdat <- matrix(c(4, rep(4,3),rep(3,4), rep(3, (3)), 5, 5, 6, 6),
               nrow = 15, ncol=7, byrow=F)
tab <- xtable::xtable(tab.all, display = c("s",rep("G", ncol(tab.all))), digits=mdat)

print.xtable(tab, type = "latex")






# figures
TC.id <- readRDS(file = "./app2/fUC_ML_TC.RDS")
TC.i1 <- readRDS(file = "./app2/Benchmark_i1_ML_TC.RDS")
TC.i2 <- readRDS(file = "./app2/Benchmark_i2_ML_TC.RDS")

hp <- mFilter::hpfilter(y, freq = 14400, type = "lambda")

# Trend
data.plot <- data.frame(
    y = TC.id$y + y0,
    trend.id  = TC.id$trend+ y0,
    cycle.id  = TC.id$cycle,
    trend.i1  = TC.i1$trend+ y0,
    cycle.i1  = TC.i1$cycle,
    trend.i2  = TC.i2$trend+ y0,
    cycle.i2  = TC.i2$cycle,
    trend.hp  = hp$trend+ y0,
    cycle.hp  = y - hp$trend, 
    time      = TC.id$time
)


library(quantmod)
getSymbols("^GSPC", from = '1927-12-30',
           to = "2024-10-31",warnings = FALSE,
           auto.assign = TRUE)

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
    y = TC.id$y + y0,
    trend.id  = TC.id$trend+ y0,
    cycle.id  = TC.id$cycle,
    trend.i1  = TC.i1$trend+ y0,
    cycle.i1  = TC.i1$cycle,
    trend.i2  = TC.i2$trend+ y0,
    cycle.i2  = TC.i2$cycle,
    trend.hp  = hp$trend+ y0,
    cycle.hp  = y - hp$trend, 
    crash = data$crash_col,
    crash_alp = data$crash_alph,
    crash2 = data$crash_col2,
    crash_alp2 = data$crash_alph2,
    crash_c = data$crash_col_c,
    crash_alp_c = data$crash_alp_c,
    
    crash_min = seq(from = as.Date("1960-01-01"), to = as.Date("2024-10-01"),  by = "month"),
    crash_max = seq(from = as.Date("1960-01-31"), to = as.Date("2024-10-31"),  by = "month"),
    time = seq(from = as.Date("1960-01-01"), to = as.Date("2024-10-01"),  by = "month"))


saveRDS(data.plot, file = "./app2/fUC_ML_TC_cbind.RDS")

cbbPalette <- wesanderson::wes_palette("Darjeeling1", n=5)

library(ggplot2)














gg1 <- ggplot(data.plot, aes(x = time, y = y)) + 
    # la nina
    geom_rect(xmin=data.plot$crash_min, 
              xmax=data.plot$crash_max, 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=data.plot$crash_alp_c, col=data.plot$crash_c) +
    labs(x = "time", y="")+
    geom_line(color="black", lwd = 0.75) + 
    ggtitle(bquote("Long-run realized volatility")) + 
    
    geom_line(mapping = aes(x = time, y = trend.id), color = cbbPalette[1], linetype = "solid", lwd  = 0.75) +
    geom_line(mapping = aes(x = time, y = trend.i1), color = cbbPalette[2], linetype = "longdash", lwd  = 1.25) +
    geom_line(mapping = aes(x = time, y = trend.i2), color = cbbPalette[3], linetype = "F1", lwd  = 1.25) +
    geom_line(mapping = aes(x = time, y = trend.hp), color = "#af8dc3", linetype = "dashed", lwd  = 1.25) +
    # geom_line(mapping = aes(x = time, y = cycle.i1), color = cbbPalette[2], linetype = "longdash", lwd  = 1) +
    # geom_line(mapping = aes(x = time, y = cycle.i2), color = cbbPalette[3], linetype = "F1", lwd  = 1) +
    # geom_line(mapping = aes(x = time, y = cycle.hp), color = "#af8dc3", linetype = "dashed", lwd  = 1) +
    # geom_line(mapping = aes(x = time, y = cycle.id), color = cbbPalette[1], lwd = 0.75, linetype = "solid") +
    # 
        #    geom_line(mapping = aes(x = time, y = tau), color = cbbPalette[1], linetype = "twodash") +
    theme_classic()
gg1


gg2 <- ggplot(data.plot, aes(x = time, y = y)) + 
    # la nina
    geom_rect(xmin=data.plot$crash_min, 
              xmax=data.plot$crash_max, 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=1, col=data.plot$crash_c) +
    labs(x = "time", y="")+
    ggtitle(bquote("Cyclical realized volatility")) + 
    geom_line(mapping = aes(x = time, y = cycle.i1), color = cbbPalette[2], linetype = "longdash", lwd  = 0.75) +
    geom_line(mapping = aes(x = time, y = cycle.i2), color = cbbPalette[3], linetype = "F1", lwd  = 0.75) +
    geom_line(mapping = aes(x = time, y = cycle.hp), color = "#af8dc3", linetype = "dashed", lwd  = 0.75) +
    geom_line(mapping = aes(x = time, y = cycle.id), color = cbbPalette[1], lwd = 1.25, linetype = "solid") +
    
    ylim(-100, 100) + 
    #    geom_line(mapping = aes(x = time, y = tau), color = cbbPalette[1], linetype = "twodash") +
    theme_classic()
gg2





KF_id <- fUC_comp(y, Est.id[1,1], diag((Est.id[2:3, 1])),
                  Est.id[4:4,1], corr=TRUE)

# checks
acf_v <- acf(KF_id$v, type = "correlation", plot = T, demean=F, lag.max=48)
data.plotting <- data.frame(autocorrelation = acf_v$acf[-1],
                            lag = 1:48)

ic_alpha= function(alpha, acf_res){
    return(qnorm((1 + (1 - alpha))/2)/sqrt(acf_res$n.used))
}
hline1 <- ic_alpha(0.05, acf_v)
hline2 <- ic_alpha(0.01, acf_v)

plot_v <- ggplot(data.plotting, aes(x = lag, y = autocorrelation)) +
    labs(x = "Lag", y="ACF")+
    geom_hline(aes(yintercept=0)) + 
    geom_segment(mapping = aes(xend = lag, yend=0))+
    ggtitle(bquote("I(d): Autocorrelation in "~v[t])) + 
    geom_hline(aes(yintercept = hline1), linetype = 2, color = cbbPalette[1])+
    geom_hline(aes(yintercept = -hline1), linetype = 2, color = cbbPalette[1])+
    geom_hline(aes(yintercept = hline2), linetype = 2, color = cbbPalette[2])+
    geom_hline(aes(yintercept = -hline2), linetype = 2, color = cbbPalette[2])+
    ylim(c(-0.06, 0.2))+
    theme_classic()





KF_i1 <- fUC_comp(y, 1, diag(c(Est.i1[2:3, 1])),
                  Est.i1[4,1], corr=TRUE)
acf_vi1 <- acf(KF_i1$v, type = "correlation", plot = T, demean=F, lag.max=48)
data.plotting.i1 <- data.frame(autocorrelation = acf_vi1$acf[-1],
                               lag = 1:48)

hline1 <- ic_alpha(0.05, acf_vi1)
hline2 <- ic_alpha(0.01, acf_vi1)

plot_v.i1 <- ggplot(data.plotting.i1, aes(x = lag, y = autocorrelation)) +
    labs(x = "Lag", y="")+
    geom_hline(aes(yintercept=0)) + 
    geom_segment(mapping = aes(xend = lag, yend=0))+
    ggtitle(bquote("I(1): Autocorrelation in "~v[t])) + 
    geom_hline(aes(yintercept = hline1), linetype = 2, color = cbbPalette[1])+
    geom_hline(aes(yintercept = -hline1), linetype = 2, color = cbbPalette[1])+
    geom_hline(aes(yintercept = hline2), linetype = 2, color = cbbPalette[2])+
    geom_hline(aes(yintercept = -hline2), linetype = 2, color = cbbPalette[2])+
    ylim(c(-0.06, 0.2))+
    theme_classic()

KF_i2 <- fUC_comp(y, 2, diag(c(Est.i2[2:3, 1])),
                  Est.i2[4:8,1], corr=TRUE)
acf_vi2 <- acf(KF_i2$v, type = "correlation", plot = T, demean=F, lag.max=48)
data.plotting.i2 <- data.frame(autocorrelation = acf_vi2$acf[-1],
                               lag = 1:48)

hline1 <- ic_alpha(0.05, acf_vi2)
hline2 <- ic_alpha(0.01, acf_vi2)

plot_v.i2 <- ggplot(data.plotting.i2, aes(x = lag, y = autocorrelation)) +
    labs(x = "Lag", y="")+
    geom_hline(aes(yintercept=0)) + 
    geom_segment(mapping = aes(xend = lag, yend=0))+
    ggtitle(bquote("I(2): Autocorrelation in "~v[t])) + 
    geom_hline(aes(yintercept = hline1), linetype = 2, color = cbbPalette[1])+
    geom_hline(aes(yintercept = -hline1), linetype = 2, color = cbbPalette[1])+
    geom_hline(aes(yintercept = hline2), linetype = 2, color = cbbPalette[2])+
    geom_hline(aes(yintercept = -hline2), linetype = 2, color = cbbPalette[2])+
    ylim(c(-0.06, 0.2))+
    theme_classic()
library(ggpubr)
ggcbind <- ggarrange(plot_v, plot_v.i1, plot_v.i2, ncol = 3)




#### store the relevant graphs
dev.off()
pdf(file = paste("/Users/tobias/Dokumente/Projekte/Filtering unknown persistence/tex/figures/app2_trend.pdf", sep=""), width = 9, height = 4)
gg1
dev.off()
ggsave("/Users/tobias/Dokumente/Projekte/Filtering unknown persistence/tex/figures/app2_trend.png", gg1, 
       width = 9, height = 4, dpi=300)

dev.off()
pdf(file = paste("/Users/tobias/Dokumente/Projekte/Filtering unknown persistence/tex/figures/app2_cycle.pdf", sep=""), width = 9, height = 4)
gg2
dev.off()
ggsave("/Users/tobias/Dokumente/Projekte/Filtering unknown persistence/tex/figures/app2_cycle.png", gg2, 
       width = 9, height = 4, dpi=300)


