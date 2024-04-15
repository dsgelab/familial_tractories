setwd("/home/ivm/Desktop/t1d")
require("survival")
library(dplyr)

date = '20240325'
data <- read.csv(paste0("cox_remove_20240305.csv"))
birthyear = 1984

cox_all <- coxph(Surv(time_to_event, T1D_EARLY) ~ group + PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 + birth_yr, data=data, weights=weight)
# summary(cox_all)

pred0 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 0))
pred1 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 1))
pred2 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 2))
pred3 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 3))

svg(paste0("surv_geno_all_",date,".svg"))
plot(pred0$time, (1-pred0$surv)*100, col = "#edadca", lwd = 3, type = "l", ylim = c(0, 15), xlab ="Children's age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "#f0758a", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#c40a0a", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "darkred", lwd = 3)
# legend("bottomleft",leg = c("0-10th: n=790","11-50th: n=3,383","51-90th: n=4,609","91-100th: n=2,005"),
# legend("topleft",leg = c("0-10th PRS","11-50th PRS","51-90th PRS","91-100th PRS"),
legend("topleft",leg = c("0-50th mid-parent PRS","50-90th mid-parent PRS","90-99th mid-parent PRS","99-100th mid-parent PRS"),
       col = c("#edadca","#f0758a","#c40a0a","darkred"), lwd = 2)
dev.off()

data_b <- data %>% filter(sex == 0)
cox_b <- coxph(Surv(time_to_event, T1D_EARLY) ~ group + PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 + birth_yr, data=data_b, weights=weight)

pred0 = survfit(cox_b, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 0))
pred1 = survfit(cox_b, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 1))
pred2 = survfit(cox_b, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 2))
pred3 = survfit(cox_b, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 3))

svg(paste0("surv_geno_boy_",date,".svg"))
plot(pred0$time, (1-pred0$surv)*100, col = "#edadca", lwd = 3, type = "l", xlim = c(0, 20), ylim = c(0, 15), xlab ="Sons' age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "#f0758a", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#c40a0a", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "darkred", lwd = 3)
legend("topleft",leg = c("0-50th mid-parent PRS","50-90th mid-parent PRS","90-99th mid-parent PRS","99-100th mid-parent PRS"),
       col = c("#edadca","#f0758a","#c40a0a","darkred"), lwd = 2)
dev.off()


data_g <- data %>% filter(sex == 1)
cox_g <- coxph(Surv(time_to_event, T1D_EARLY) ~ group + PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 + birth_yr, data=data_g, weights=weight)

pred0 = survfit(cox_g, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 0))
pred1 = survfit(cox_g, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 1))
pred2 = survfit(cox_g, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 2))
pred3 = survfit(cox_g, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 3))

svg(paste0("surv_geno_girl_",date,".svg"))
plot(pred0$time, (1-pred0$surv)*100, col = "#edadca", lwd = 3, type = "l", xlim = c(0, 20), ylim = c(0, 15), xlab ="Daughters' age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "#f0758a", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#c40a0a", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "darkred", lwd = 3)
legend("topleft",leg = c("0-50th mid-parent PRS","50-90th mid-parent PRS","90-99th mid-parent PRS","99-100th mid-parent PRS"),
       col = c("#edadca","#f0758a","#c40a0a","darkred"), lwd = 2)
dev.off()







cox_all <- coxph(Surv(time_to_event, T1D_EARLY) ~ group + PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 +BL_AGE +BL_YEAR, data=data, weights=weight)
pred0 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 0))
pred1 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 1))
pred2 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 2))
pred3 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 3))
svg(paste0("surv_geno_all_",date,".svg"))
plot(pred0$time, (1-pred0$surv)*100, col = "#edadca", lwd = 3, type = "l", ylim = c(0, 15), xlab ="Children's age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "#f0758a", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#c40a0a", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "darkred", lwd = 3)
# legend("bottomleft",leg = c("0-10th: n=790","11-50th: n=3,383","51-90th: n=4,609","91-100th: n=2,005"),
# legend("topleft",leg = c("0-10th PRS","11-50th PRS","51-90th PRS","91-100th PRS"),
legend("topleft",leg = c("0-50th mid-parent Full-PGS","50-90th mid-parent Full-PGS","90-99th mid-parent Full-PGS","99-100th mid-parent Full-PGS"),
       col = c("#edadca","#f0758a","#c40a0a","darkred"), lwd = 2)
dev.off()

data_b <- data %>% filter(sex == 0)
cox_b <- coxph(Surv(time_to_event, T1D_EARLY) ~ group + PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 +BL_AGE +BL_YEAR, data=data_b, weights=weight)
pred0 = survfit(cox_b, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 0))
pred1 = survfit(cox_b, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 1))
pred2 = survfit(cox_b, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 2))
pred3 = survfit(cox_b, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 3))
svg(paste0("surv_geno_boy_",date,".svg"))
plot(pred0$time, (1-pred0$surv)*100, col = "#edadca", lwd = 3, type = "l", xlim = c(0, 20), ylim = c(0, 15), xlab ="Sons' age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "#f0758a", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#c40a0a", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "darkred", lwd = 3)
legend("topleft",leg = c("0-50th mid-parent Full-PGS","50-90th mid-parent Full-PGS","90-99th mid-parent Full-PGS","99-100th mid-parent Full-PGS"),
       col = c("#edadca","#f0758a","#c40a0a","darkred"), lwd = 2)
dev.off()

data_g <- data %>% filter(sex == 1)
cox_g <- coxph(Surv(time_to_event, T1D_EARLY) ~ group + PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 +BL_AGE +BL_YEAR, data=data_g, weights=weight)
pred0 = survfit(cox_g, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 0))
pred1 = survfit(cox_g, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 1))
pred2 = survfit(cox_g, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 2))
pred3 = survfit(cox_g, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 3))
svg(paste0("surv_geno_girl_",date,".svg"))
plot(pred0$time, (1-pred0$surv)*100, col = "#edadca", lwd = 3, type = "l", xlim = c(0, 20), ylim = c(0, 15), xlab ="Daughters' age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "#f0758a", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#c40a0a", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "darkred", lwd = 3)
legend("topleft",leg = c("0-50th mid-parent Full-PGS","50-90th mid-parent Full-PGS","90-99th mid-parent Full-PGS","99-100th mid-parent Full-PGS"),
       col = c("#edadca","#f0758a","#c40a0a","darkred"), lwd = 2)
dev.off()






cox_all <- coxph(Surv(time_to_event, T1D_EARLY) ~ pa_improved_T1D_STRICT + PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 +BL_AGE +BL_YEAR, data=data, weights=weight)
summary(cox_all)

pred0 = survfit(cox_all, newdata = data.frame(PC1 = -0.000212, PC2 = -0.000451, PC3 = 0.000289, PC4 = 0.000283, PC5 = -0.000149, PC6 = 0.000866,
PC7 = -0.000040, PC8 = 0.000184, PC9 = 0.000372, PC10 = 0.000239, BL_AGE = 28.53, BL_YEAR = 2005, group = 0))
pred1 = survfit(cox_all, newdata = data.frame(PC1 = -0.000212, PC2 = -0.000451, PC3 = 0.000289, PC4 = 0.000283, PC5 = -0.000149, PC6 = 0.000866,
PC7 = -0.000040, PC8 = 0.000184, PC9 = 0.000372, PC10 = 0.000239, BL_AGE = 28.53, BL_YEAR = 2005, group = 1))
pred2 = survfit(cox_all, newdata = data.frame(PC1 = -0.000212, PC2 = -0.000451, PC3 = 0.000289, PC4 = 0.000283, PC5 = -0.000149, PC6 = 0.000866,
PC7 = -0.000040, PC8 = 0.000184, PC9 = 0.000372, PC10 = 0.000239, BL_AGE = 28.53, BL_YEAR = 2005, group = 2))
pred3 = survfit(cox_all, newdata = data.frame(PC1 = -0.000212, PC2 = -0.000451, PC3 = 0.000289, PC4 = 0.000283, PC5 = -0.000149, PC6 = 0.000866,
PC7 = -0.000040, PC8 = 0.000184, PC9 = 0.000372, PC10 = 0.000239, BL_AGE = 28.53, BL_YEAR = 2005, group = 3))

svg(paste0("surv_geno_all_",date,".svg"))
plot(pred0$time, (1-pred0$surv)*100, col = "#edadca", lwd = 3, type = "l", ylim = c(0, 10), xlab ="Children's age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "#f0758a", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#c40a0a", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "darkred", lwd = 3)
# legend("bottomleft",leg = c("0-10th: n=790","11-50th: n=3,383","51-90th: n=4,609","91-100th: n=2,005"),
# legend("topleft",leg = c("0-10th PRS","11-50th PRS","51-90th PRS","91-100th PRS"),
legend("topleft",leg = c("0-50th PRS","50-90th PRS","90-99th PRS","99-100th PRS"),
       col = c("#edadca","#f0758a","#c40a0a","darkred"), lwd = 2)
dev.off()


pred0 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, g1 = 1, g2 = 0, g3 = 0, g4 = 0))
pred1 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, g1 = 0, g2 = 1, g3 = 0, g4 = 0))
pred2 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, g1 = 0, g2 = 0, g3 = 1, g4 = 0))
pred3 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, g1 = 0, g2 = 0, g3 = 0, g4 = 1))

pred0 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 0))
pred1 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 1))
pred2 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 2))
pred3 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 3))

pred0 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 20, BL_YEAR = 2004, group = 0))
pred1 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 20, BL_YEAR = 2004, group = 1))
pred2 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 20, BL_YEAR = 2004, group = 2))
pred3 = survfit(cox_all, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2004, group = 3))

svg(paste0("surv_geno_boy_",date,".svg"))
plot(pred0$time, (1-pred0$surv)*100, col = "#edadca", lwd = 3, type = "l", ylim = c(0, 20), xlab ="Children's age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "#f0758a", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#c40a0a", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "darkred", lwd = 3)
legend("topleft",leg = c("0-50th PRS","50-90th PRS","90-99th PRS","99-100th PRS"),
       col = c("#edadca","#f0758a","#c40a0a","darkred"), lwd = 2)
dev.off()


data_g <- data %>% filter(sex == 1)
cox_g <- coxph(Surv(time_to_event, T1D_EARLY) ~ group + PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 +BL_AGE +BL_YEAR, data=data_g, weights=weight)
summary(cox_g)

pred0 = survfit(cox_g, newdata = data.frame(PC1 = -0.000123, PC2 = -0.000792, PC3 = 0.000334, PC4 = 0.000324, PC5 = -0.000180, PC6 = 0.000932,
PC7 = 0.000008, PC8 = 0.000246, PC9 = 0.000345, PC10 = 0.000275, BL_AGE = 29.81, BL_YEAR = 2010, group = 0), weights=weight)
pred1 = survfit(cox_g, newdata = data.frame(PC1 = -0.000123, PC2 = -0.000792, PC3 = 0.000334, PC4 = 0.000324, PC5 = -0.000180, PC6 = 0.000932,
PC7 = 0.000008, PC8 = 0.000246, PC9 = 0.000345, PC10 = 0.000275, BL_AGE = 29.81, BL_YEAR = 2010, group = 1), weights=weight)
pred2 = survfit(cox_g, newdata = data.frame(PC1 = -0.000123, PC2 = -0.000792, PC3 = 0.000334, PC4 = 0.000324, PC5 = -0.000180, PC6 = 0.000932,
PC7 = 0.000008, PC8 = 0.000246, PC9 = 0.000345, PC10 = 0.000275, BL_AGE = 29.81, BL_YEAR = 2010, group = 2), weights=weight)
pred3 = survfit(cox_g, newdata = data.frame(PC1 = -0.000123, PC2 = -0.000792, PC3 = 0.000334, PC4 = 0.000324, PC5 = -0.000180, PC6 = 0.000932,
PC7 = 0.000008, PC8 = 0.000246, PC9 = 0.000345, PC10 = 0.000275, BL_AGE = 29.81, BL_YEAR = 2010, group = 3), weights=weight)

pred0 = survfit(cox_g, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 0))
pred1 = survfit(cox_g, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 1))
pred2 = survfit(cox_g, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 2))
pred3 = survfit(cox_g, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, BL_AGE = 21, BL_YEAR = 2005, group = 3))

svg(paste0("surv_geno_girl_",date,".svg"))
plot(pred0$time, (1-pred0$surv)*100, col = "#edadca", lwd = 3, type = "l", xlim = c(0, 20), ylim = c(0, 10), xlab ="Daughters' age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "#f0758a", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#c40a0a", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "darkred", lwd = 3)
#legend("bottomleft",leg = c("0-10th: n=488","11-50th: n=2,103","51-90th: n=2,738","91-100th: n=1,070"),
# legend("topleft",leg = c("0-10th PRS","11-50th PRS","51-90th PRS","91-100th PRS"),
legend("topleft",leg = c("0-50th PRS","50-90th PRS","90-99th PRS","99-100th PRS"),
       col = c("#edadca","#f0758a","#c40a0a","darkred"), lwd = 2)
dev.off()


data_g <- data %>% filter(sex == 0)
cox_g <- coxph(Surv(time_to_event, T1D_EARLY) ~ group + PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 +BL_AGE +BL_YEAR, data=data_g, weights=weight)
summary(cox_g)

pred0 = survfit(cox_g, newdata = data.frame(PC1 = -0.000440, PC2 = 0.000092, PC3 = 0.000205, PC4 = 0.000194, PC5 = -0.000118, PC6 = 0.000736,
PC7 = -0.000126, PC8 = 0.000106, PC9 = 0.000410, PC10 = 0.000194, BL_AGE = 26.135, BL_YEAR = 2003, group = 0), weights=weight)
pred1 = survfit(cox_g, newdata = data.frame(PC1 = -0.000440, PC2 = 0.000092, PC3 = 0.000205, PC4 = 0.000194, PC5 = -0.000118, PC6 = 0.000736,
PC7 = -0.000126, PC8 = 0.000106, PC9 = 0.000410, PC10 = 0.000194, BL_AGE = 26.135, BL_YEAR = 2003, group = 1), weights=weight)
pred2 = survfit(cox_g, newdata = data.frame(PC1 = -0.000440, PC2 = 0.000092, PC3 = 0.000205, PC4 = 0.000194, PC5 = -0.000118, PC6 = 0.000736,
PC7 = -0.000126, PC8 = 0.000106, PC9 = 0.000410, PC10 = 0.000194, BL_AGE = 26.135, BL_YEAR = 2003, group = 2), weights=weight)
pred3 = survfit(cox_g, newdata = data.frame(PC1 = -0.000440, PC2 = 0.000092, PC3 = 0.000205, PC4 = 0.000194, PC5 = -0.000118, PC6 = 0.000736,
PC7 = -0.000126, PC8 = 0.000106, PC9 = 0.000410, PC10 = 0.000194, BL_AGE = 26.135, BL_YEAR = 2003, group = 3), weights=weight)

svg(paste0("surv_geno_boy_",date,".svg"))
plot(pred0$time, (1-pred0$surv)*100, col = "#edadca", lwd = 3, type = "l", xlim = c(0, 20), ylim = c(0, 10), xlab ="Sons' age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "#f0758a", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#c40a0a", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "darkred", lwd = 3)
# legend("bottomleft",leg = c("0-10th: n=302","11-50th: n=1,280","51-90th: n=1,871","91-100th: n=935"),
# legend("topleft",leg = c("0-10th PRS","11-50th PRS","51-90th PRS","91-100th PRS"),
legend("topleft",leg = c("0-50th PRS","50-90th PRS","90-99th PRS","99-100th PRS"),
       col = c("#edadca","#f0758a","#c40a0a","darkred"), lwd = 2)
dev.off()


setwd("/home/ivm/Desktop/t1d")
library(survival)
library(survminer)

date = '20240128'
data <- read.csv(paste0("cox_data_",date,".csv"))
data_fr <- read.csv("df_fr.csv")
fit2 <- coxph(Surv(time_to_event, T1D_EARLY) ~ PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 +BL_AGE +BL_YEAR + strata(group), data = data)
summary(fit2)

# conditional balancing in groups
ggadjustedcurves(fit2, data = data, method = "conditional", variable = "group", reference = data_fr)+coord_cartesian(ylim = c(0.8, 1))
curve <- surv_adjustedcurves(fit2, data = data, method = "conditional", variable = "group", reference = data_fr)


data_date = '20240202' #'20240122'
plot_date = '20240202' #'20240122'
data <- read.csv(paste0("family_",date,".csv"))
birthyear = 1984
y_len = 1.5

plotting <- function(model, x_label, disease){
    pred0 = survfit(model, newdata = data.frame(ch_year = birthyear, group_t1d = 0))
    pred1 = survfit(model, newdata = data.frame(ch_year = birthyear, group_t1d = 1))
    pred2 = survfit(model, newdata = data.frame(ch_year = birthyear, group_t1d = 2))
    pred3 = survfit(model, newdata = data.frame(ch_year = birthyear, group_t1d = 3))

    svg(paste0("surv_reg_all_",plot_date,".svg"))
    plot(pred0$time, (1-pred0$surv)*100, col = "#006e97", lwd = 3, type = "l", xlim = c(0, 20), ylim = c(0, y_len),
        xlab = x_label, ylab = "Cumulative incidence of T1D (%)")
    lines(pred1$time, (1-pred1$surv)*100, col = "cadetblue", lwd = 3)
    lines(pred2$time, (1-pred2$surv)*100, col = "#ff8c00", lwd = 3)
    lines(pred3$time, (1-pred3$surv)*100, col = "#8b0000", lwd = 3)
    legend("topleft",leg = c(paste0("No parental ", disease), paste0("Maternal ", disease),
        paste0("Paternal ", disease), paste0("Paternal and maternal ", disease)),
        col = c("#006e97","cadetblue","#ff8c00","#8b0000"), lwd = 2)
    dev.off()
}