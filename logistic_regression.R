library(survival)
library(MatchIt)

setwd("~/familial_analysis")
data <- read.csv("data_py.csv")

results <- NULL

note = "all"
for (number in seq(1,length(eps))){
        # loop for the diseases from mothers
        if (eps[number] %in% ep_remove_mo) {
                results <- rbind(results,c(eps[number],"mother",NaN,NaN,NaN,NaN,NaN,note))
        } else {
                ep_col_name <- paste('mo','_ep',as.character(number-1), sep = '')
                data <- select(m.data, ep_col_name, "ch_year", "mo_year", "fa_year", "number_of_sib", "province", "ch_ep0", "subclass")
                names(data)[names(data) == ep_col_name] <- "exposure"
                model <- clogit(ch_ep0 ~ exposure + ch_year + mo_year + fa_year + number_of_sib + province + strata(subclass), data, method = "exact")
                se <- summary(model)$coeff["exposure","se(coef)"]
                pval <- summary(model)$coeff["exposure","Pr(>|z|)"]
                hr <- summary(model)$conf.int["exposure","exp(coef)"]
                hr_025 <- summary(model)$conf.int["exposure","lower .95"]
                hr_975 <- summary(model)$conf.int["exposure","upper .95"]
                results <- rbind(results,c(eps[number],"mother",se,pval,hr,hr_025,hr_975,note))
        }
        # loop for the diseases from fathers
        if (eps[number] %in% ep_remove_fa) {
                results <- rbind(results,c(eps[number],"father",NaN,NaN,NaN,NaN,NaN,note))
        } else {
                ep_col_name <- paste('fa','_ep',as.character(number-1), sep = '')
                data <- select(m.data, ep_col_name, "ch_year", "mo_year", "fa_year", "number_of_sib", "province", "ch_ep0", "subclass")
                names(data)[names(data) == ep_col_name] <- "exposure"
                model <- clogit(ch_ep0 ~ exposure + ch_year + mo_year + fa_year + number_of_sib + province + strata(subclass), data, method = "exact")
                se <- summary(model)$coeff["exposure","se(coef)"]
                pval <- summary(model)$coeff["exposure","Pr(>|z|)"]
                hr <- summary(model)$conf.int["exposure","exp(coef)"]
                hr_025 <- summary(model)$conf.int["exposure","lower .95"]
                hr_975 <- summary(model)$conf.int["exposure","upper .95"]
                results <- rbind(results,c(eps[number],"father",se,pval,hr,hr_025,hr_975,note))
        }
}

note = "boy"
for (number in seq(1,length(eps))){
        # loop for the diseases from mothers
        if (eps[number] %in% ep_remove_mo) {
                results <- rbind(results,c(eps[number],"mother",NaN,NaN,NaN,NaN,NaN,note))
        } else {
                ep_col_name <- paste('mo','_ep',as.character(number-1), sep = '')
                data <- m.data %>% filter(sex==0) %>% select(ep_col_name, "ch_year", "mo_year", "fa_year", "number_of_sib", "province", "ch_ep0", "subclass")
                names(data)[names(data) == ep_col_name] <- "exposure"
                if (sum(data$exposure) < 20) {
                        results <- rbind(results,c(eps[number],"mother",NaN,NaN,NaN,NaN,NaN,note))
                } else {
                        model <- clogit(ch_ep0 ~ exposure + ch_year + mo_year + fa_year + number_of_sib + province + strata(subclass), data, method = "exact")
                        se <- summary(model)$coeff["exposure","se(coef)"]
                        pval <- summary(model)$coeff["exposure","Pr(>|z|)"]
                        hr <- summary(model)$conf.int["exposure","exp(coef)"]
                        hr_025 <- summary(model)$conf.int["exposure","lower .95"]
                        hr_975 <- summary(model)$conf.int["exposure","upper .95"]
                        results <- rbind(results,c(eps[number],"mother",se,pval,hr,hr_025,hr_975,note))
                }

        }
        # loop for the diseases from fathers
        if (eps[number] %in% ep_remove_fa) {
                results <- rbind(results,c(eps[number],"father",NaN,NaN,NaN,NaN,NaN,note))
        } else {
                ep_col_name <- paste('fa','_ep',as.character(number-1), sep = '')
                data <- m.data %>% filter(sex==0) %>% select(ep_col_name, "ch_year", "mo_year", "fa_year", "number_of_sib", "province", "ch_ep0", "subclass")
                names(data)[names(data) == ep_col_name] <- "exposure"
                if (sum(data$exposure) < 20) {
                        results <- rbind(results,c(eps[number],"mother",NaN,NaN,NaN,NaN,NaN,note))
                } else {
                        model <- clogit(ch_ep0 ~ exposure + ch_year + mo_year + fa_year + number_of_sib + province + strata(subclass), data, method = "exact")
                        se <- summary(model)$coeff["exposure","se(coef)"]
                        pval <- summary(model)$coeff["exposure","Pr(>|z|)"]
                        hr <- summary(model)$conf.int["exposure","exp(coef)"]
                        hr_025 <- summary(model)$conf.int["exposure","lower .95"]
                        hr_975 <- summary(model)$conf.int["exposure","upper .95"]
                        results <- rbind(results,c(eps[number],"father",se,pval,hr,hr_025,hr_975,note))
                }
        }
}


note = "girl"
for (number in seq(1,length(eps))){
        # loop for the diseases from mothers
        if (eps[number] %in% ep_remove_mo) {
                results <- rbind(results,c(eps[number],"mother",NaN,NaN,NaN,NaN,NaN,note))
        } else {
                ep_col_name <- paste('mo','_ep',as.character(number-1), sep = '')
                data <- m.data %>% filter(sex==1) %>% select(ep_col_name, "ch_year", "mo_year", "fa_year", "number_of_sib", "province", "ch_ep0", "subclass")
                names(data)[names(data) == ep_col_name] <- "exposure"
                if (sum(data$exposure) < 20) {
                        results <- rbind(results,c(eps[number],"mother",NaN,NaN,NaN,NaN,NaN,note))
                } else {
                        model <- clogit(ch_ep0 ~ exposure + ch_year + mo_year + fa_year + number_of_sib + province + strata(subclass), data, method = "exact")
                        se <- summary(model)$coeff["exposure","se(coef)"]
                        pval <- summary(model)$coeff["exposure","Pr(>|z|)"]
                        hr <- summary(model)$conf.int["exposure","exp(coef)"]
                        hr_025 <- summary(model)$conf.int["exposure","lower .95"]
                        hr_975 <- summary(model)$conf.int["exposure","upper .95"]
                        results <- rbind(results,c(eps[number],"mother",se,pval,hr,hr_025,hr_975,note))
                }

        }
        # loop for the diseases from fathers
        if (eps[number] %in% ep_remove_fa) {
                results <- rbind(results,c(eps[number],"father",NaN,NaN,NaN,NaN,NaN,note))
        } else {
                ep_col_name <- paste('fa','_ep',as.character(number-1), sep = '')
                data <- m.data %>% filter(sex==1) %>% select(ep_col_name, "ch_year", "mo_year", "fa_year", "number_of_sib", "province", "ch_ep0", "subclass")
                names(data)[names(data) == ep_col_name] <- "exposure"
                if (sum(data$exposure) < 20) {
                        results <- rbind(results,c(eps[number],"mother",NaN,NaN,NaN,NaN,NaN,note))
                } else {
                        model <- clogit(ch_ep0 ~ exposure + ch_year + mo_year + fa_year + number_of_sib + province + strata(subclass), data, method = "exact")
                        se <- summary(model)$coeff["exposure","se(coef)"]
                        pval <- summary(model)$coeff["exposure","Pr(>|z|)"]
                        hr <- summary(model)$conf.int["exposure","exp(coef)"]
                        hr_025 <- summary(model)$conf.int["exposure","lower .95"]
                        hr_975 <- summary(model)$conf.int["exposure","upper .95"]
                        results <- rbind(results,c(eps[number],"father",se,pval,hr,hr_025,hr_975,note))
                }
        }
}