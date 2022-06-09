setwd("~/familial_analysis")
library(dplyr)
library(survival)
library(rjson)

OUTCOME = "T1D_STRICT"
m.data <- read.csv(paste0(OUTCOME,"/data_",OUTCOME,".csv"))
m.data$subclass <- as.factor(m.data$subclass)
eps <- fromJSON(file=paste0(OUTCOME,"/eps_",OUTCOME,".json"))


get_stats <- function(who, number , note, dataframe) {
    ep_col_name <- paste0(substring(who, 1, 2),'_ep',as.character(number-1))
    data <- select(dataframe, ep_col_name, "ch_year", "mo_year", "fa_year", "number_of_sib", "outcome", "subclass")
    names(data)[names(data) == ep_col_name] <- "exposure"
    n_cases <- sum(data$exposure)
    if (n_cases < 20) {
        res <- c(eps[number],who,NaN,NaN,NaN,NaN,NaN,note,n_cases)
    } else {
        model <- clogit(outcome ~ exposure + ch_year + mo_year + fa_year + number_of_sib + strata(subclass), data, method = "exact")
        se <- summary(model)$coeff["exposure","se(coef)"]
        pval <- summary(model)$coeff["exposure","Pr(>|z|)"]
        hr <- summary(model)$conf.int["exposure","exp(coef)"]
        hr_025 <- summary(model)$conf.int["exposure","lower .95"]
        hr_975 <- summary(model)$conf.int["exposure","upper .95"]
        res <- c(eps[number],who,se,pval,hr,hr_025,hr_975,note,n_cases)
    }
    return(res)
}


results <- NULL

note = "all"
for (number in seq(1,length(eps))){
    # loop for the diseases from mothers
    results <- rbind(results,get_stats("mother", number , note, m.data))
    # loop for the diseases from fathers
    results <- rbind(results,get_stats("father", number , note, m.data))
}

note = "boy"
for (number in seq(1,length(eps))){
    # loop for the diseases from mothers
    results <- rbind(results,get_stats("mother", number , note, m.data %>% filter(sex==0)))
    # loop for the diseases from fathers
    results <- rbind(results,get_stats("father", number , note, m.data %>% filter(sex==0)))
}

note = "girl"
for (number in seq(1,length(eps))){
    # loop for the diseases from mothers
    results <- rbind(results,get_stats("mother", number , note, m.data %>% filter(sex==1)))
    # loop for the diseases from fathers
    results <- rbind(results,get_stats("father", number , note, m.data %>% filter(sex==1)))
}

RES <- data.frame(results)
colnames(RES) <- c("endpoint","who","se","pval","hr","hr_025","hr_975","note","n_cases")
write.csv(RES,file=paste0(OUTCOME, "/results_",OUTCOME,"_r.csv"), row.names = FALSE)
