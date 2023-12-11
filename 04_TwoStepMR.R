###### Generic Script used to our Two Step MR analysis ######

# For each protein having colocalized with the outcomes (AAM or ANM) (H4 > 0.8), 
# we test whether the BMI is a mediator of this relationship

# For AAM outcome, we use pediatric BMI GWAS
# For ANM outcome, we use adult BMI GWAS
# See our Manuscript, to have the link to download the GWAS

# Mediation function (adapt from bda package)
mediation <- function(a, sa, b, sb) {
  
  # estimate indirect effect
  ab <- a * b
  
  # sobel formule
  tmp1 <- b^2 * sa^2 + a^2 * sb^2
  tmp2 <- sa^2 * sb^2
  
  # calculate z-score and pval
  zsob <- a * b/sqrt(tmp1)
  psob <- pnorm(-abs(zsob)) * 2
  
  # calculate se from beta and z-score
  se <- ab / zsob
  
  # calculate CI 95%
  CI_sup <- ab + se * qnorm(0.975)
  CI_inf <- ab - se * qnorm(0.975)
  
  res <- c(ab, se, CI_inf, CI_sup, zsob, psob)
  names(res) <- c("ab", "se", "CI_inf", "CI_sup", "z", "pv") 
  
  return(res)
}

# Need the following MR results:

# MR : protein --> BMI
# same script 01_TwoSampleMR but you need to modify the outcome
mr_a <- "MR : prot -> BMI"

# MR : BMI --> Outcome (AAM or ANM)
# 
# # simple MR study between to GWAS (simplify code)
# g1 <- "GWAS BMI"
# g1 <- filter(g1, p_val <= 5e-8) 
# g1 <- format_data(g1)
# g1 <- clump_data(g1)
# g2 <- "GWAS outcome"
# g2 <- format_data(g2)
# g <- harmonise(g1, g2)
# mr_b <- mr(g)

# only keep IVW 
mr_a <- filter(mr_a, method == "Inverse variance weighted")
mr_b <- filter(mr_b, method == "Inverse variance weighted")

# arrange mr_b
mr_b <- mr_b[, 5:9]
colnames(mr_b) <- c("B_method", "B_nsnp", "B_beta", "B_se", "B_pval")

# arrange mr_a
mr_a <- mr_a[, c(1:3, 8:12)]
colnames(mr_a)[4:8] <- c("A_method", "A_nsnp", "A_beta", "A_se", "A_pval")  

# paste mr
mr <- cbind(mr_a, mr_b)

# estimate mediation effect from mr_a and mr_b
med <- NULL

for (i in 1:nrow(mr)) {
  tmp <- mr[i, ]
  tmp <- mediation(a = tmp$A_beta, sa = tmp$A_se, b = tmp$B_beta, sb = tmp$B_se)
  
  med <- rbind(med, tmp)
}

med <- as.data.frame(med)

# all together 
RES <- cbind(mr, med)

# saveRDS(tab_res, "data/mediation....")
