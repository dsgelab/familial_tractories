setwd("/home/ivm/Desktop/t1d")
require("survival")
library(dplyr)

data_date = '20240202'
plot_date = '20240325'
birthyear = 1984 # median birth year of the cohort
y_len = 15

plotting <- function(model, file_name, x_label){
    pred0 = survfit(model, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
    PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 0))
    pred1 = survfit(model, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
    PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 1))
    pred2 = survfit(model, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
    PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 2))
    pred3 = survfit(model, newdata = data.frame(PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0, PC6 = 0,
    PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0, birth_yr = birthyear, group = 3))

    svg(paste0(file_name,plot_date,".svg"))
    plot(pred0$time, (1-pred0$surv)*100, col = "#edadca", lwd = 3, type = "l", xlim = c(0, 20), ylim = c(0, y_len),
        xlab = x_label, ylab = "Cumulative incidence of T1D (%)")
    lines(pred1$time, (1-pred1$surv)*100, col = "#f0758a", lwd = 3)
    lines(pred2$time, (1-pred2$surv)*100, col = "#c40a0a", lwd = 3)
    lines(pred3$time, (1-pred3$surv)*100, col = "darkred", lwd = 3)
    legend("topleft",leg = c("0-50th mid-parent PRS", "50-90th mid-parent PRS", "90-99th mid-parent PRS",
        "99-100th mid-parent PRS"), col = c("#edadca","#f0758a","#c40a0a","darkred"), lwd = 2)
    dev.off()
}

data <- read.csv(paste0("cox_remove_",data_date,".csv"))
cox_all <- coxph(Surv(time_to_event, T1D_EARLY) ~ group + PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 + birth_yr,
    data=data, weights=weight)
# summary(cox_all)
plotting(cox_all, "surv_geno_all_", "Children's age")

data_b <- data %>% filter(sex == 0)
cox_b <- coxph(Surv(time_to_event, T1D_EARLY) ~ group + PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 + birth_yr, data=data_b, weights=weight)
plotting(cox_b, "surv_geno_boy_", "Sons' age")

data_g <- data %>% filter(sex == 1)
cox_g <- coxph(Surv(time_to_event, T1D_EARLY) ~ group + PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 + birth_yr, data=data_g, weights=weight)
plotting(cox_g, "surv_geno_girl_", "Daughters' age")