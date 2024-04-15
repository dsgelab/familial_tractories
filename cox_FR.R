setwd("~/familial_analysis")
require("survival")
library(dplyr)

date = '20240202' #'20240122'
data <- read.csv(paste0("family_",date,".csv"))
birthyear = 1984
y_len = 1.5

cox_all <- coxph(Surv(time_to_event, T1D_EARLY) ~ group_t1d + ch_year, data=data)
#summary(cox_all)
#median(data$ch_year)

pred0 = survfit(cox_all, newdata = data.frame(ch_year = birthyear, group_t1d = 0))
pred1 = survfit(cox_all, newdata = data.frame(ch_year = birthyear, group_t1d = 1))
pred2 = survfit(cox_all, newdata = data.frame(ch_year = birthyear, group_t1d = 2))
pred3 = survfit(cox_all, newdata = data.frame(ch_year = birthyear, group_t1d = 3))

svg("surv_reg_all_20240306.svg")
plot(pred0$time, (1-pred0$surv)*100, col = "#006e97", lwd = 3, type = "l", xlim = c(0, 20), ylim = c(0, y_len), xlab ="Children's age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "cadetblue", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#ff8c00", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "#8b0000", lwd = 3)
#legend("bottomleft",leg = c("No T1D: n=1,863,055","Maternal: n=2,586","Paternal: n=3,403","Both T1D: n<15"),
legend("topleft",leg = c("No parental T1D","Maternal T1D","Paternal T1D","Paternal and maternal T1D"),
       col = c("#006e97","cadetblue","#ff8c00","#8b0000"), lwd = 2)
dev.off()

data_g <- data %>% filter(sex == 1)
cox_g <- coxph(Surv(time_to_event, T1D_EARLY) ~ group_t1d + ch_year, data=data_g)
#summary(cox_g)

pred0 = survfit(cox_g, newdata = data.frame(ch_year = birthyear, group_t1d = 0))
pred1 = survfit(cox_g, newdata = data.frame(ch_year = birthyear, group_t1d = 1))
pred2 = survfit(cox_g, newdata = data.frame(ch_year = birthyear, group_t1d = 2))
pred3 = survfit(cox_g, newdata = data.frame(ch_year = birthyear, group_t1d = 3))

svg("surv_reg_girl_20240306.svg")
plot(pred0$time, (1-pred0$surv)*100, col = "#006e97", lwd = 3, type = "l", xlim = c(0, 20), ylim = c(0, y_len), xlab ="Daughters' age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "cadetblue", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#ff8c00", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "#8b0000", lwd = 3)
legend("topleft",leg = c("No parental T1D","Maternal T1D","Paternal T1D","Paternal and maternal T1D"),
       col = c("#006e97","cadetblue","#ff8c00","#8b0000"), lwd = 2)
dev.off()


data_b <- data %>% filter(sex == 0)
cox_b <- coxph(Surv(time_to_event, T1D_EARLY) ~ group_t1d + ch_year, data=data_b)
#summary(cox_b)
#median(data_b$ch_year)

pred0 = survfit(cox_b, newdata = data.frame(ch_year = birthyear, group_t1d = 0))
pred1 = survfit(cox_b, newdata = data.frame(ch_year = birthyear, group_t1d = 1))
pred2 = survfit(cox_b, newdata = data.frame(ch_year = birthyear, group_t1d = 2))
pred3 = survfit(cox_b, newdata = data.frame(ch_year = birthyear, group_t1d = 3))

svg("surv_reg_boy_20240306.svg")
plot(pred0$time, (1-pred0$surv)*100, col = "#006e97", lwd = 3, type = "l", xlim = c(0, 20), ylim = c(0, y_len), xlab ="Sons' age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "cadetblue", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#ff8c00", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "#8b0000", lwd = 3)
legend("topleft",leg = c("No parental T1D","Maternal T1D","Paternal T1D","Paternal and maternal T1D"),
       col = c("#006e97","cadetblue","#ff8c00","#8b0000"), lwd = 2)
dev.off()



data_a <- data %>% filter(group_t1d == 0)
cox_all <- coxph(Surv(time_to_event, T1D_EARLY) ~ group_aid + ch_year, data=data_a)
pred0 = survfit(cox_all, newdata = data.frame(ch_year = birthyear, group_aid = 0))
pred1 = survfit(cox_all, newdata = data.frame(ch_year = birthyear, group_aid = 1))
pred2 = survfit(cox_all, newdata = data.frame(ch_year = birthyear, group_aid = 2))
pred3 = survfit(cox_all, newdata = data.frame(ch_year = birthyear, group_aid = 3))
svg("surv_aid_all_20240306.svg")
plot(pred0$time, (1-pred0$surv)*100, col = "#006e97", lwd = 3, type = "l", xlim = c(0, 20), ylim = c(0, y_len), xlab ="Children's age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "cadetblue", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#ff8c00", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "#8b0000", lwd = 3)
legend("topleft",leg = c("No parental AID","Maternal AID(s)","Paternal AID(s)","Paternal and maternal AID(s)"),
       col = c("#006e97","cadetblue","#ff8c00","#8b0000"), lwd = 2)
dev.off()


data_g <- data_a %>% filter(sex == 1)
cox_g <- coxph(Surv(time_to_event, T1D_EARLY) ~ group_aid + ch_year, data=data_g)
pred0 = survfit(cox_g, newdata = data.frame(ch_year = birthyear, group_aid = 0))
pred1 = survfit(cox_g, newdata = data.frame(ch_year = birthyear, group_aid = 1))
pred2 = survfit(cox_g, newdata = data.frame(ch_year = birthyear, group_aid = 2))
pred3 = survfit(cox_g, newdata = data.frame(ch_year = birthyear, group_aid = 3))
svg("surv_aid_girl_20240306.svg")
plot(pred0$time, (1-pred0$surv)*100, col = "#006e97", lwd = 3, type = "l", xlim = c(0, 20), ylim = c(0, y_len), xlab ="Daughters' age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "cadetblue", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#ff8c00", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "#8b0000", lwd = 3)
legend("topleft",leg = c("No parental AID","Maternal AID(s)","Paternal AID(s)","Paternal and maternal AID(s)"),
       col = c("#006e97","cadetblue","#ff8c00","#8b0000"), lwd = 2)
dev.off()


data_g <- data_a %>% filter(sex == 0)
cox_g <- coxph(Surv(time_to_event, T1D_EARLY) ~ group_aid + ch_year, data=data_g)
pred0 = survfit(cox_g, newdata = data.frame(ch_year = birthyear, group_aid = 0))
pred1 = survfit(cox_g, newdata = data.frame(ch_year = birthyear, group_aid = 1))
pred2 = survfit(cox_g, newdata = data.frame(ch_year = birthyear, group_aid = 2))
pred3 = survfit(cox_g, newdata = data.frame(ch_year = birthyear, group_aid = 3))
svg("surv_aid_boy_20240306.svg")
plot(pred0$time, (1-pred0$surv)*100, col = "#006e97", lwd = 3, type = "l", xlim = c(0, 20), ylim = c(0, y_len), xlab ="Sons' age", ylab = "Cumulative incidence of T1D (%)")
lines(pred1$time, (1-pred1$surv)*100, col = "cadetblue", lwd = 3)
lines(pred2$time, (1-pred2$surv)*100, col = "#ff8c00", lwd = 3)
lines(pred3$time, (1-pred3$surv)*100, col = "#8b0000", lwd = 3)
legend("topleft",leg = c("No parental AID","Maternal AID(s)","Paternal AID(s)","Paternal and maternal AID(s)"),
       col = c("#006e97","cadetblue","#ff8c00","#8b0000"), lwd = 2)
dev.off()