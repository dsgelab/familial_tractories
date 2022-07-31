############################################################
#      Format and perform logistic regression with R       #
############################################################

setwd("~/Desktop")
library(data.table)
library(R.utils)
library(dplyr)
library("stringr")

'%!in%' <- function(x,y)!('%in%'(x,y))


## read in and format HLA haplotype data --------------------
hla_path <- 'finngen/shared/r9_hla_imputation_version2/20220226_012206/files/ritarij/R9_imputed_HLAs_v2.vcf.gz'
data <- fread(hla_path, sep="\t", header=T)
dim(data)   # 187 356,222

# only keep useful variables ----------
d_names <- colnames(d)
d_keep_names <- c("ID", d_names[substr(d_names,1,2)=="FG"])
d <- d %>% select(d_keep_names)
dim(d)   # 33 392,658

# transpose to a data frame with less columns ----------
dat <- t(d)
hla_alleles <- gsub("\\*","_",gsub(":","_",gsub("-","_",dat[1,])))
colnames(dat) <- hla_alleles
dat <- dat[-1,]
dat <- data.frame(dat)

dat[,"ID"] <- rownames(dat)
dim(dat)  # 356,213     33

# from vcf format to genotype ----------
for (hla_sub_i in hla_alleles) {
	print(hla_sub_i)
	dat[,hla_sub_i] <- substr(dat[,hla_sub_i],1,3)
	dat[,paste0(hla_sub_i,"_GP")] <- ifelse(dat[,hla_sub_i]=="0/0",0,ifelse(dat[,hla_sub_i]=="1/1",2,ifelse(dat[,hla_sub_i] %in% c("1/0","0/1"),1,NA)))
}
dim(dat)  # 356,213     67
save(dat, file=paste0("finngen_R9_imputed_HLA.Rdata"))


## summary stats of HLA for the manuscript ------------
dat <- data.frame(get(load(paste0("finngen_R9_imputed_HLA.Rdata"))))
dim(dat)   # 392,649    429

dat_names <- colnames(dat)[substr(colnames(dat),1,3)=="HLA"]
hla_alleles <- dat_names[substr(dat_names,nchar(dat_names)-2,nchar(dat_names)) %in% c("_GP")]

# matrix to store names of original GP, locus, and gene
hla_info <- matrix(NA, nrow=length(hla_alleles), ncol=3)
colnames(hla_info) <- c("original", "locus", "gene")
hla_info[,"original"] <- hla_alleles
hla_info[,"locus"] <- gsub("_GP", "", hla_alleles)


# filtering
hla_info <- hla_info %>% data.frame() %>%
                  mutate(original=as.character(original), locus=as.character(locus), gene=as.character(gene)) %>%
                  filter(substr(locus, nchar(locus)-2,nchar(locus)) %!in% c("_ng")) %>%
                  filter(str_count(locus, "_")==3) %>%
                  filter(locus!="HLA_DRB4_01_03N")

# add gene name
hla_genes <- gsub("HLA_","", hla_info$locus)
hla_genes <- unlist(str_split(hla_genes, "_", n=2))
hla_info[,"gene"] <- hla_genes[substr(hla_genes,1,1) %!in% 0:9]
nrow(hla_info)   # 172


# format locus name
hla_info <- hla_info %>% mutate(locus=sub("HLA_","",locus), locus=sub("_","*",locus), locus=sub("_",":",locus))
11:54
########################################################################################
#                 Association with mLOX, mLOY, and Autosomal mCA                       #
########################################################################################

mod_NA <- matrix(NA, ncol=19, nrow=1)
colnames(mod_NA) <- c("Estimate", "SE", "Z", "P_val", "OR", "OR_025", "OR_975",
                      "Model", "Exposure", "Population", "Outcome",
                      "N_mCACases", "N_mCAControls",
                      "N_mCACases_HLA00",    "N_mCACases_HLA01",    "N_mCACases_HLA11",
                      "N_mCAControls_HLA00", "N_mCAControls_HLA01", "N_mCAControls_HLA11")


## read in and format HLA haplotype data --------------------
# for (mCA in c("LOX", "LOY", "Autosomal_mCA", "Autosomal_mCA_female", "Autosomal_mCA_male")) {
for (mCA in c("LOX")) {
	if (mCA=="LOX"){ pheno <- read.table("Pheno_LOX_NoCorrectSmoke_PhenoCov.r9.20220218.tsv", header=T) }
	if (mCA=="LOY"){ pheno <- read.table("Pheno_LOY_NoCorrectSmoke_PhenoCov.2Mb.r9.20220225.tsv", header=T) }
	if (substr(mCA, 1, 13)=="Autosomal_mCA"){ pheno <- read.table("FinnGenR9Affymetrix.Autosomal_mCA.withPhe.r9.20220508.tsv", header=T) }

	summaryDF <- data.frame()
	dat <- data.frame(get(load(paste0("finngen_R9_imputed_HLA.Rdata")))) %>%
	       inner_join(pheno, by=c("ID"="FINNGENID")) %>%
	       rename(mCA=substr(mCA, 1, 13)) %>% filter(!is.na(mCA))
	dim(dat)   # 168,838    525

	if (mCA=="Autosomal_mCA_female") { dat <- dat %>% filter(SEX=="female") }
	if (mCA=="Autosomal_mCA_male") { dat <- dat %>% filter(SEX=="male") }


	for (hla_sub_i in hla_info$original) {
		print(hla_sub_i)
		if (length(unique(dat[ ,hla_sub_i]))>=2){
			if (mCA=="Autosomal_mCA") {
				mylogit <- glm(mCA ~ dat[ ,hla_sub_i] + SEX + BL_AGE + BL_AGE2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family="binomial", data=dat)
			} else {
				mylogit <- glm(mCA ~ dat[ ,hla_sub_i] + BL_AGE + BL_AGE2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family="binomial", data=dat)
			}

			mod <- as.data.frame(t(as.data.frame(summary(mylogit)$coeff["dat[, hla_sub_i]",])))
			colnames(mod) <- c("Estimate", "SE", "Z", "P_val")
			mod$OR <- exp(mod[1,1])
			ci <- exp(confint(mylogit,"dat[, hla_sub_i]"))
			mod$OR_025 <- ci[1]
			mod$OR_975 <- ci[2]
		} else {
			mod <- mod_NA
		}

		if (mCA=="Autosomal_mCA") {
			mod[,"Model"] <- "mCA=HLA+SEX+BL_AGE+BL_AGE_2+PC1_10"
		} else {
			mod[,"Model"] <- "mCA=HLA+BL_AGE+BL_AGE_2+PC1_10"
		}
		mod[,"Exposure"] <- gsub("_GP","",hla_sub_i)
		mod[,"Population"] <- "All"
		mod[,"Outcome"] <- mCA

		mod[,"N_mCACases"] <- length(which(dat[,"mCA"]==1))
		mod[,"N_mCAControls"] <- length(which(dat[,"mCA"]==0))
		mod[,"N_mCACases_HLA00"] <- length(which(dat[,"mCA"]==1 & dat[,hla_sub_i]==0))
		mod[,"N_mCACases_HLA01"] <- length(which(dat[,"mCA"]==1 & dat[,hla_sub_i]==1))
		mod[,"N_mCACases_HLA11"] <- length(which(dat[,"mCA"]==1 & dat[,hla_sub_i]==2))
		mod[,"N_mCAControls_HLA00"] <- length(which(dat[,"mCA"]==0 & dat[,hla_sub_i]==0))
		mod[,"N_mCAControls_HLA01"] <- length(which(dat[,"mCA"]==0 & dat[,hla_sub_i]==1))
		mod[,"N_mCAControls_HLA11"] <- length(which(dat[,"mCA"]==0 & dat[,hla_sub_i]==2))
		rownames(mod) <- "HLA"
		print(mod)
		summaryDF <- rbind(summaryDF, mod)
	}
	summaryDF %>% mutate(P_val=as.numeric(as.character(P_val)), OR=round(as.numeric(as.character(OR)),2)) %>% filter(P_val<0.05/206) %>% arrange(-desc(P_val)) %>% select(Exposure,P_val,OR)
	write.table(summaryDF, paste0("FineMap_HLA.",mCA,".finngenR9.20220524.tsv"), append=F, quote=F, sep="\t", row.names=F, col.names=T)
}