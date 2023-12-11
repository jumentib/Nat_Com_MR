###### Generic Script used to our MVMR analysis ######

# For each protein having colocalized with the outcomes (AAM or ANM) (H4 > 0.8)

library(MVMR)

# For each MVMR, we take all the significant SNPs for the two exposures (Proteins and BMI): 
# the leading SNP in the protein GWAS 
# as well as all the significant SNPs in the BMI GWAS (pediatric BMI for AAM and adult BMI for ANM)
# Finally we find these sets of SNPs in the GWAS of outcomes (AAM and ANM)

# Exposure 1 :GWAS Protein level 
# contains SNPs leader of protein GWAS + all the significant SNPs of BMI GWAS
# column contains : rsid, beta, se
expo1 <- "read files"
prot <- "protein level tested"

# Exposure 2 : GWAS BMI (pediatric of adult)
# contains SNPs leader of protein GWAS + all the significant SNPs of BMI GWAS
# column contains : rsid, beta, se
expo2 <- "read files"

# Outcome : GWAS AAM or ANM 
# contains SNPs leader of protein GWAS + all the significant SNPs of BMI GWAS
# column contains : rsid, beta, se
outc <- "read files"

# merge data frame
# be careful to the order of rsid
tmp <- cbind(expo1, expo2, outc)

# format data for MVMR
fdat <- format_mvmr(BXGs = tmp[, c("prot_beta", "bmi_beta")],
                    BYG = tmp[, "out_beta"], 
                    seBXGs = tmp[, c("prot_se", "bmi_se")], 
                    seBYG = tmp[, "out_se"], 
                    RSID = "rsid")


res <- as.data.frame(ivw_mvmr(r_input = fdat))

# add Condidence Interval (95%)
res$CI_sup <- res$Estimate + res$`Std. Error` * qnorm(0.975)
res$CI_inf <- res$Estimate - res$`Std. Error` * qnorm(0.975)

# add protein level name
res$prot <- prot

# saveRDS(res, "data/mvmr....")

