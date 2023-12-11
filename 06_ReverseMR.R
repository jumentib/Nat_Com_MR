###### Generic Script used to our Reverse MR analysis ######

# For each protein having colocalized with the outcomes (AAM or ANM) (H4 > 0.8)

library(TwoSampleMR)
library(tidyverse)

# Reverse MR : AAM (or ANM) as exposure and protein level as outcome

# Load full GWAS of exposure (AAM or ANM)
expo <- "read GWAS of exposure"

# Select SNPs with significant effect (p < 5e-8)
expo <- filter(expo, P_value <= 5e-8)

# Format data for TwoSampleMR package (need to be careful with the name of each colunm)
expo <- format_data(expo)

# Clump the dataset of exposure (independant significant effect of SNP)
expo <- clump_data(expo)

# set some name
protI <- "protein name"

# Load full GWAS of outcome (protein level)
outc <- "read GWAS outcome"

# Format outcome data (need to be careful with the name of each colunm)
outc <- format_data(outc, type = "outcome")

# Harmonise data between exposure and outcome GWAS
dat1 <- harmonise_data(exposure_dat = expo, outcome_dat = outc, action = 1)

# Lauch MR of dat1
res1 <- mr(dat1)

# save the result
tmp <- list(mr = res1,
            harmonise = dat1,
            prot = protI)

# saveRDS(tmp, "data/reverse....")