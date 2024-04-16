install.packages("/finngen/green/gap.datasets_0.0.5.tar.gz", 
                 "~/R/x86_64-pc-linux-gnu-library/4.2", 
                 repos=NULL, type="source")
install.packages("/finngen/green/gap_1.3-1.tar.gz", 
                 "~/R/x86_64-pc-linux-gnu-library/4.2", 
                 repos=NULL, type="source")
install.packages("/finngen/green/polspline_1.1.22.tar.gz", 
                 "~/R/x86_64-pc-linux-gnu-library/4.2", 
                 repos=NULL, type="source")
install.packages("/finngen/green/Hmisc_4.7-2.tar.gz", 
                 "~/R/x86_64-pc-linux-gnu-library/4.2", 
                 repos=NULL, type="source")
install.packages("/finngen/green/rms_6.3-0.tar.gz", 
                 "~/R/x86_64-pc-linux-gnu-library/4.2", 
                 repos=NULL, type="source")
install.packages("/finngen/green/arsenal_3.6.3.tar.gz", 
                 "~/R/x86_64-pc-linux-gnu-library/4.2", 
                 repos=NULL, type="source")
install.packages("/finngen/green/haplo.stats_1.9.2.tar.gz", 
                 "~/R/x86_64-pc-linux-gnu-library/4.2", 
                 repos=NULL, type="source")
library(gap)
library(haplo.stats)

data(fsnps)

hap(id=fsnps[,1],data=fsnps[,3:10],nloci=4)
dir()
hap.control <- function(mb=0,pr=0,po=0.001,to=0.001,th=1,maxit=100,n=0,
                        ss=0,rs=0,rp=0,ro=0,rv=0,sd=0,mm=0,mi=0,mc=50,ds=0.1,de=0,q=0,
                        hapfile="hap.out",assignfile="assign.out") {
  list(mb=mb,pr=pr,po=po,to=to,th=th,maxit=maxit,n=n,ss=ss,rs=rs,
       rp=rp,ro=ro,rv=rv,sd=sd,mm=mm,mi=mi,mc=mc,ds=ds,de=de,q=q,
       hapfile=hapfile,assignfile=assignfile)
}
control <- hap.control(ss=1,mi=5,hapfile="h",assignfile="a")
hap(id=fsnps[,1],data=fsnps[,3:10],nloci=4,control=control)

hap(id=hla[,1],data=hla[,3:8],nloci=4,control=control)

data(hla)
y<-hla[,2]
geno<-hla[,3:8]
# complete data
hap.score(y,geno,locus.label=c("DRB","DQA","DQB"))

library(gap)
library(haplo.stats)
data <- read.csv("/home/ivm/Desktop/hla_data_for_hap.csv")
y = data$T1D_STRICT

geno = data[, 2:15]
rm(data)
# incomplete genotype data
genes = c("A","B","C","DPB1","DQA1","DQB1","DRB1")
#genes = c("A","B","C","DPB1","DQA1","DQB1","DRB1","DRB345")
res = hap.score(y,geno,locus.label=genes)
#hap.score(y,geno,locus.label=genes,miss.val=NA,handle.miss=1,mloci=1)
res
#unlink("assign.dat")

geno = data[, 4:7]
genes = c("B","C")
res = hap.score(y,geno,locus.label=genes)
res

geno = data[, c(2,3,6,7)]
genes = c("A","C")
res = hap.score(y,geno,locus.label=genes)
res

geno = data[, 12:15]
genes = c("DQB1","DRB1")
res = hap.score(y,geno,locus.label=genes)
res

geno = data[, c(8,9,14,15)]
genes = c("DPB1","DRB1")
res = hap.score(y,geno,locus.label=genes)
res

geno = data[, c(4,5,6,7,10,11,12,13,14,15)]
genes = c("B","C", "DQA1","DQB1","DRB1")
res = hap.score(y,geno,locus.label=genes)
res

data <- read.csv("/home/ivm/Desktop/hla_data_for_hap1.csv")
geno = data[, 52:57]
genes = c("DQA1","DQB1","DRB1")
y = data$M13_RHEUMA
res = hap.score(y,geno,locus.label=genes,miss.val=NA)#,handle.miss=1,mloci=1)
res


res = hap.score(y,geno,
                trait.type = "binomial", # for binary outcome
                x.adj = data[,c(32:43)], # non-genetic covariates to adjust the score stat
                locus.label=genes,miss.val=NA)#,handle.miss=1,mloci=1)




data <- read.csv("/home/ivm/Desktop/hla_data_for_hap1.csv")
geno = data[, 52:57]
genes = c("DQA1","DQB1","DRB1")
y = data$T1D_STRICT
trait.int = 2 # for "binomial"
x.adj = as.matrix(data[,c(32:43)])


miss <- apply(is.na(geno),1,any)
if(!all(is.na(0))) {
  for(mval in 0){
    miss <- miss | apply(geno==mval, 1, any)
  }
}
x.adj <- as.matrix(x.adj)

miss <- miss | is.na(y) 
miss <- miss| apply(is.na(x.adj),1,any)

y <- as.numeric(y[!miss])
geno <- geno[!miss,]
x.adj <- x.adj[!miss,,drop=F]

n.subj = length(y)

#haplo <- gc.em(data=geno, NA, converge.eps=0.00001, 
#               maxiter=5000, handle.miss=0, miss.val=0)
# a multiallelic version of progressive EM algorithm for haplotype inference
haplo <- hap.em(1:n.subj, NA, data=geno, converge.eps=0.00001, 
                maxiter=5000, miss.val=0)
nreps = as.vector(haplo$nreps)
# Vector of posterior probabilities of pairs of haplotypes for a person, given thier marker phenotypes.
post <- haplo$post
hap1 <- haplo$hap1code
hap2 <- haplo$hap2code
indx <- haplo$indx.subj
uhap<-haplo$uhap
# Skip score statistics for haplotypes with frequencies < 0.005.
which.haplo<-haplo$hap.prob>=0.005
uhap<-uhap[which.haplo]
x <- outer(hap1,uhap,"==") + outer(hap2,uhap,"==")
n.x <- ncol(x)
x.post<-matrix(rep(NA, n.subj * n.x), ncol=n.x)
for(j in 1:n.x){
  x.post[,j] <- tapply(x[,j]*post, indx, sum)
}


reg.out <- glm(y ~ x.adj, family="binomial")
x.adj <- cbind(rep(1,n.subj),x.adj)
mu <- reg.out$fitted.values
a=1
v <- mu*(1-mu)


u.mtx  <- (y-mu)*x.post / a
u.score <- apply(u.mtx,2,sum)

# Var matrix for x.adj covariates
v.11 <- t(x.adj * v) %*% x.adj

# Var matrix for covar(x.adj, x.post)
v.21 <- t(x.post) %*% (x.adj * v)

# Var matrix for haplo scores
res <- ( (y - mu)/a ) ^2
t1 <- rep( (v-res) ,nreps) * post
v.22 <- t(x*t1) %*% x + t(u.mtx) %*% u.mtx

# Var matrix for haplo scores, adjusted for x.adj
v.score <- v.22 - v.21 %*% solve(v.11) %*% t(v.21) 



tmp <- haplo.stats::Ginv(v.score)
df <- tmp$rank
g.inv <- tmp$Ginv
score.global <- u.score%*% g.inv %*%u.score
score.haplo <- u.score / sqrt(diag(v.score))
score.max <-  max(score.haplo^2)

score.global.p <- 1 - pchisq(score.global,df)
score.haplo.p <- 1-pchisq(score.haplo^2,1)




library(gap)
library(haplo.stats)

data <- read.csv("/home/ivm/Desktop/hla_data_for_hap1.csv")
trait.int = 2 # for "binomial"
numbers_df = read.csv("~/num_df.csv")

gene_list = c("A","B","C","DPB1","DQA1","DQB1","DRB1")
col_dict = list("A"=c(44,45), "B"=c(46,47), "C"=c(48,49), "DPB1"=c(50,51), 
                "DQA1"=c(52,53), "DQB1"=c(54,55), "DRB1"=c(56,57))

combinitions = lapply(seq_along(gene_list), function(i) combn(gene_list, i, FUN = list))
combinitions = unlist(combinitions, recursive=FALSE) # n=127

possi_haplo = c()
poten_haplo = c()
global_scores = c()
RES = NULL

# 69-70 75 76 91 99 100 101 105-107 109-111 120-122 124 125 127
for (i in c(71:73,77:93)) {
  tryCatch(
    expr = {
      curr_list = combinitions[[i]]
      col_list=c()
      for (num in c(1:length(curr_list))) {
        col_list = c(col_list, unlist(col_dict[curr_list[num]]))
        
      }
      geno = data[,col_list]
      genes = curr_list
      
      #res = hap.score(y,geno,locus.label=curr_list)
      y = data$T1D_STRICT
      x.adj = as.matrix(data[,c(32:43)])
      
      miss <- apply(is.na(geno),1,any)
      if(!all(is.na(0))) {
        for(mval in 0){
          miss <- miss | apply(geno==mval, 1, any)
        }
      }
      x.adj <- as.matrix(x.adj)
      
      miss <- miss | is.na(y) 
      miss <- miss| apply(is.na(x.adj),1,any)
      
      y <- as.numeric(y[!miss])
      n.subj = length(y)
      geno <- geno[!miss,]
      x.adj <- x.adj[!miss,,drop=F]
      
      haplo <- gc.em(data=geno, NA, converge.eps=0.00001, 
                     maxiter=5000, handle.miss=0, miss.val=0)
      possi_haplo = c(possi_haplo,c(length(haplo$hap.prob)))
      # a multiallelic version of progressive EM algorithm for haplotype inference
      #haplo <- hap.em(1:n.subj, NA, data=geno, converge.eps=0.00001, 
      #                maxiter=5000, miss.val=0)
      nreps = as.vector(haplo$nreps)
      # Vector of posterior probabilities of pairs of haplotypes for a person, given thier marker phenotypes.
      post <- haplo$post
      hap1 <- haplo$hap1code
      hap2 <- haplo$hap2code
      indx <- haplo$indx.subj
      uhap<-haplo$uhap
      # Skip score statistics for haplotypes with frequencies < 0.005.
      which.haplo<-haplo$hap.prob>=0.005
      uhap<-uhap[which.haplo]
      poten_haplo = c(poten_haplo,c(length(uhap)))
      x <- outer(hap1,uhap,"==") + outer(hap2,uhap,"==")
      n.x <- ncol(x)
      x.post<-matrix(rep(NA, n.subj * n.x), ncol=n.x)
      for(j in 1:n.x){
        x.post[,j] <- tapply(x[,j]*post, indx, sum)
      }
      
      
      reg.out <- glm(y ~ x.adj, family="binomial")
      x.adj <- cbind(rep(1,n.subj),x.adj)
      mu <- reg.out$fitted.values
      a=1
      v <- mu*(1-mu)
      
      
      u.mtx  <- (y-mu)*x.post / a
      u.score <- apply(u.mtx,2,sum)
      
      # Var matrix for x.adj covariates
      v.11 <- t(x.adj * v) %*% x.adj
      
      # Var matrix for covar(x.adj, x.post)
      v.21 <- t(x.post) %*% (x.adj * v)
      
      # Var matrix for haplo scores
      res <- ( (y - mu)/a ) ^2
      t1 <- rep( (v-res) ,nreps) * post
      v.22 <- t(x*t1) %*% x + t(u.mtx) %*% u.mtx
      
      # Var matrix for haplo scores, adjusted for x.adj
      v.score <- v.22 - v.21 %*% solve(v.11) %*% t(v.21) 
      
      
      tmp <- haplo.stats::Ginv(v.score)
      df <- tmp$rank
      g.inv <- tmp$Ginv
      score.global <- u.score%*% g.inv %*%u.score
      global_scores = c(global_scores,score.global)
      score.haplo <- u.score / sqrt(diag(v.score))
      score.max <-  max(score.haplo^2)
      
      score.global.p <- 1 - pchisq(score.global,df)
      score.haplo.p <- 1-pchisq(score.haplo^2,1)
      
      haplos=haplo$haplotype[which.haplo,]
      rownames(haplos) = 1:nrow(haplos)
      colnames(haplos) = curr_list
      
      haplotype = data.frame(matrix(ncol=length(gene_list)+3,nrow=length(uhap), 
                                    dimnames=list(NULL,c(gene_list,
                                                         'Freq','Score','pval'))))
      
      for (j in gene_list) {
        if (j %in% curr_list) {
          # col_name = paste0('loc-',as.character(which(curr_list == 'A')))
          # haplo$haplotype[which.haplo,][col_name]
          haplotype[j] = haplos[,j]
        } 
      }
      haplotype['pval'] = score.haplo.p
      haplotype['Freq'] = haplo$hap.prob[which.haplo]
      haplotype['Score'] = score.haplo
      # haplotype = cbind(haplotype, )
      
      RES = rbind(RES, haplotype)
      new_row = c(possi_haplo=c(length(haplo$hap.prob)),
                  poten_haplo=c(length(uhap)),
                  global_scores=score.global,
                  haplosss=paste(combinitions[[i]], collapse='-'))
      numbers_df = rbind(numbers_df, new_row)
      print(paste(i, "is done."))
    }, error = function(e){
      new_row = c(possi_haplo=0,poten_haplo=0,global_scores=0,
                  haplosss=paste(combinitions[[i]], collapse='-'))
      numbers_df = rbind(numbers_df, new_row)
      print(paste(i, paste(combinitions[[i]], collapse='-')))
    }
    
  )
}

# 1-pchisq(score.haplo^2,1, lower.tail = FALSE, log.p = TRUE)/log(10)

write.csv(RES,file="~/Desktop/t1d/complete_hap_res.csv", row.names = FALSE)


numbers_df = data.frame(haplosss=sapply(combinitions[69:100], paste, collapse='-'), 
                        possi_haplo=possi_haplo, 
                        poten_haplo=poten_haplo, 
                        global_scores=global_scores)
write.csv(numbers_df,file="Desktop/t1d/num_df.csv", row.names = FALSE)

new_row = c(haplosss="A-B-DPB1-DQA1", 
            possi_haplo=0,poten_haplo=0,global_scores=0)
numbers_df = rbind(numbers_df, new_row)


a = read.csv("num_df.csv")
a$haplosss = sapply(combinitions[8:50], paste, collapse='-')
a = rbind(a, numbers_df)
write.csv(a,file="num_df.csv", row.names = FALSE)


b = read.csv("complete_hap_res.csv")
b = rbind(b, RES)
write.csv(b,file="complete_hap_res.csv", row.names = FALSE)
