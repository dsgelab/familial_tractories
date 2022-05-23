"""METHOD 2: disease risk score matching"""

install.packages('MatchIt',repo='file://data/cran/')
library(survival)
library(MatchIt)
setwd("~/familial_analysis")

MATCH_NUM = 3

df <- read.csv("df.csv")
m.out <- matchit(ch_ep0 ~ sex + ch_year + mo_year + fa_year + sib_number + province,
                 data = df, method = "nearest", ratio = MATCH_NUM)
summary(m.out) # for numerical summaries
plot(m.out) # for graphical summaries
m.data <- match.data(m.out)
z.out <- zelig(ch_ep0 ~ mo_ep1 + x1 + x2, model = mymodel, data = m.data)
# model <- clogit(ch_ep0 ~ mo_ep1 + age + age2 + strata(['parent_id']), data, weight, strata, method="exact")
