"""METHOD 2: disease risk score matching"""

install.packages('MatchIt',repo='file://data/cran/')
library(survival)
library(MatchIt)

setwd("~/familial_analysis")
df <- read.csv("df.csv")
m.out <- matchit(ch_ep0 ~ sex + ch_year + mo_year + fa_year + sib_number + province, data = df, method = "nearest", ratio = 3)
summary(m.out) # for numerical summaries
plot(m.out) # for graphical summaries
m.data <- match.data(m.out)
summary(m.data$ch_ep0)
write.csv(m.data, 'data_r.csv')

eps <- c('T1D_STRICT', 'M13_RHEUMA', 'M13_RELAPSPOLYCHONDR', 'M13_SJOGREN', 'M13_SYSTSLCE', 'M13_DERMATOPOLY', # E4_DM1
       'M13_WEGENER', 'M13_MICROPOLYANG', 'M13_CHURGSTRAUSS', 'D3_ALLERGPURPURA', 'M13_BEHCET', 'M13_MCTD',
       'M13_HYPERANG', 'SLE_FG',  #'M13_SLE',
       'I9_RHEUFEV', 'G6_MS', 'G6_ADEM', 'G6_DISSOTH', 'G6_NARCOCATA', 'AUTOIMMUNE_HYPERTHYROIDISM',
       'E4_THYROIDITAUTOIM', 'E4_AUTOPOLYFAI', 'E4_HYTHY_AI_STRICT', 'E4_GRAVES_OPHT_STRICT', 'E4_ADDISON',
       'AUTOHEP', 'D3_AIHA_DRUG', 'D3_AIHA_OTHER', 'D3_ITP', 'D3_ANAEMIA_B12_DEF', 'K11_COELIAC', 'K11_IBD',
       'G6_MYASTHENIA', 'G6_OTHDEMYEL', 'G6_MYOMUSCINOTH', 'G6_GUILBAR', 'H7_IRIDOCYC_ANTER',  'CHIRBIL_PRIM',
       'L12_PSORIASIS', 'L12_VITILIGO', 'L12_ALOPECAREATA', 'L12_PEMPHIGOID', 'L12_DERMATHERP',
       'N14_HENOCHSCHONLEIN_NEPHRITIS', 'N14_IGA_NEPHROPATHY', 'T2D', 'GEST_DIABETES')
ep_remove_mo <- c('M13_RELAPSPOLYCHONDR',
                  'M13_MICROPOLYANG',
                  'M13_CHURGSTRAUSS',
                  'M13_BEHCET',
                  'M13_HYPERANG',
                  'G6_ADEM',
                  'G6_NARCOCATA',
                  'E4_AUTOPOLYFAI',
                  'AUTOHEP',
                  'D3_AIHA_DRUG',
                  'G6_MYOMUSCINOTH',
                  'N14_HENOCHSCHONLEIN_NEPHRITIS')
ep_remove_fa <- c('M13_RELAPSPOLYCHONDR',
                  'M13_MICROPOLYANG',
                  'M13_CHURGSTRAUSS',
                  'M13_BEHCET',
                  'M13_HYPERANG',
                  'G6_ADEM',
                  'G6_DISSOTH',
                  'G6_NARCOCATA',
                  'E4_THYROIDITAUTOIM',
                  'E4_AUTOPOLYFAI',
                  'E4_GRAVES_OPHT_STRICT',
                  'AUTOHEP',
                  'D3_AIHA_DRUG',
                  'G6_MYOMUSCINOTH',
                  'N14_HENOCHSCHONLEIN_NEPHRITIS',
                  'GEST_DIABETES')

results <- NULL

# use the whole matched data to do analysis
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

# given parents' disease, remove children who had the same disease
note = "ch_no_pa_ep"
for (number in seq(1,length(eps))){
        # loop for the diseases from mothers
        if (eps[number] %in% ep_remove_mo) {
                results <- rbind(results,c(eps[number],"mother",NaN,NaN,NaN,NaN,NaN,note))
        } else {
                ep_col_name <- paste('mo','_ep',as.character(number-1), sep = '')
                if (ep_col_name != 'mo_ep0') {
                        ch_col_name <- paste('ch','_ep',as.character(number-1), sep = '')
                        data <- select(m.data, ep_col_name, "ch_year", "mo_year", "fa_year", "number_of_sib", "province", "ch_ep0", "subclass", ch_col_name)
                        names(data)[names(data) == ep_col_name] <- "exposure"
                        names(data)[names(data) == ch_col_name] <- "exposure_ch"
                        data_to_remove <- data %>% filter((exposure==1) & (exposure_ch==1)) %>% select(subclass)
                        data <- data %>% filter(!subclass %in% data_to_remove)
                } else {
                        data <- select(m.data, ep_col_name, "ch_year", "mo_year", "fa_year", "number_of_sib", "province", "ch_ep0", "subclass")
                        names(data)[names(data) == ep_col_name] <- "exposure"
                }
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
                if (ep_col_name != 'fa_ep0') {
                        ch_col_name <- paste('ch','_ep',as.character(number-1), sep = '')
                        data <- select(m.data, ep_col_name, "ch_year", "mo_year", "fa_year", "number_of_sib", "province", "ch_ep0", "subclass", ch_col_name)
                        names(data)[names(data) == ep_col_name] <- "exposure"
                        names(data)[names(data) == ch_col_name] <- "exposure_ch"
                        data_to_remove <- data %>% filter((exposure==1) & (exposure_ch==1)) %>% select(subclass)
                        data <- data %>% filter(!subclass %in% data_to_remove)
                } else {
                        data <- select(m.data, ep_col_name, "ch_year", "mo_year", "fa_year", "number_of_sib", "province", "ch_ep0", "subclass")
                        names(data)[names(data) == ep_col_name] <- "exposure"
                }
                model <- clogit(ch_ep0 ~ exposure + ch_year + mo_year + fa_year + number_of_sib + province + strata(subclass), data, method = "exact")
                se <- summary(model)$coeff["exposure","se(coef)"]
                pval <- summary(model)$coeff["exposure","Pr(>|z|)"]
                hr <- summary(model)$conf.int["exposure","exp(coef)"]
                hr_025 <- summary(model)$conf.int["exposure","lower .95"]
                hr_975 <- summary(model)$conf.int["exposure","upper .95"]
                results <- rbind(results,c(eps[number],"father",se,pval,hr,hr_025,hr_975,note))
        }
}


RES <- data.frame(results)
colnames(RES) <- c("endpoint","who","se","pval","hr","hr_025","hr_975","note")
write.csv(RES,file="results_ep0.csv", row.names = FALSE)

