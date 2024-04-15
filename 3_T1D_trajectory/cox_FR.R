setwd("~/familial_analysis")
require("survival")
library(dplyr)

data_date = '20240202' #'20240122'
plot_date = '20240312'
birthyear = 1984 # median birth year of the cohort
y_len = 15

plotting <- function(model, file_name, x_label, disease){
    pred0 = survfit(model, newdata = data.frame(ch_year = birthyear, group_t1d = 0))
    pred1 = survfit(model, newdata = data.frame(ch_year = birthyear, group_t1d = 1))
    pred2 = survfit(model, newdata = data.frame(ch_year = birthyear, group_t1d = 2))
    pred3 = survfit(model, newdata = data.frame(ch_year = birthyear, group_t1d = 3))

    svg(paste0("surv_",file_name,"_",plot_date,".svg"))
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

data <- read.csv(paste0("family_",data_date,".csv"))
cox_all <- coxph(Surv(time_to_event, T1D_EARLY) ~ group_t1d + ch_year, data=data)
#summary(cox_all)
#median(data$ch_year)
plotting(cox_all, "t1d_all", "Children's age", "T1D")

data_g <- data %>% filter(sex == 1)
cox_g <- coxph(Surv(time_to_event, T1D_EARLY) ~ group_t1d + ch_year, data=data_g)
plotting(cox_g, "t1d_girl", "Daughters' age", "T1D")

data_b <- data %>% filter(sex == 0)
cox_b <- coxph(Surv(time_to_event, T1D_EARLY) ~ group_t1d + ch_year, data=data_b)
plotting(cox_b, "t1d_boy", "Sons' age", "T1D")

data_a <- data %>% filter(group_t1d == 0)
cox_all <- coxph(Surv(time_to_event, T1D_EARLY) ~ group_aid + ch_year, data=data_a)
plotting(cox_all, "aid_all", "Children's age", "AID(s)")

data_g <- data_a %>% filter(sex == 1)
cox_g <- coxph(Surv(time_to_event, T1D_EARLY) ~ group_aid + ch_year, data=data_g)
plotting(cox_g, "aid_girl", "Daughters' age", "AID(s)")

data_b <- data_a %>% filter(sex == 0)
cox_b <- coxph(Surv(time_to_event, T1D_EARLY) ~ group_aid + ch_year, data=data_b)
plotting(cox_b, "aid_boy", "Sons' age", "AID(s)")