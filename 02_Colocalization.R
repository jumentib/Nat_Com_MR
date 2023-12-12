###### Generic Script used to do all of our Colocalization analysis ######

# Load Library 

library(coloc)
library(tidyverse)
library(GWAS.utils)
library(reshape2)
library(susieR)

# We have lauch this script for each significant PROTEIN in our MR study

# save the result
tab_res <- data.frame()

# TO MODIFY ----------------------
# 
# GET THE RIGHT FILES OF GWAS
# GET THE RIGHT RSID OF INTEREST
# GET THE RIGHT GENE
# GET THE RIGHT CHR
# 
# exposure GWAS (protein levels)
# need to download all the GWAS of significant PROTEIN in our MR study 
# See our Manuscript, to have the link to download the GWAS
files <- "files GWAS PROT"

# rsid of interest
rsidI <- "NAME OF RSID OF INTEREST"
protI <- "NAME OF PROTEIN"
chrI <- "XXX"

# studyI <- "Emilsson" # NAME of GWAS STUDY

# TO NOT MODIFY ----------------------
# 
# everything is automatic 
# 
# load exposure GWAS
expo <- read.table(files)

# find the position of the rsid of interest
posI <- expo$Pos[expo$rsids == rsidI]

# 1000000 Mb windows
pos_min <- posI - 500000
pos_max <- posI + 500000

# cut the GWAS of exposure
expo_cut <- expo[expo$Pos >= pos_min & expo$Pos <= pos_max, ]


# load GWAS of Outcome (AAM or ANM) => GWAS data (corresponding chr)
outc <- readRDS(paste0("GWAS OUTCOME (previously cut by chr)", chrI, ".rds"))


# merge the 2 GWAS by rsid 
tmp <- merge.data.frame(expo_cut, outc, by.x = 3, by.y = 2)

# rearrange data
tmp <- tmp[, c(1:4, 7, 5, 8, 9, 10, 13, 16, 14)]

colnames(tmp) <- c("rsid", "chr", "position", "expo_eff", "expo_se", "expo_pv", "expo_N", "expo_maf",  
                   "outc_eff", "outc_se", "outc_pv")

# varbeta = square(se)
tmp$expo_varbeta <- (tmp$expo_se)^2
tmp$outc_varbeta <- (tmp$outc_se)^2

# print RSID of Interest
tmp[tmp$rsid == rsidI, ]

# print the min pValue of the outcome
tmp[which.min(tmp$outc_pv), ]

# PLOT GWAS  

# data for ggplot (rsid, position, pv_expo, pv_outc)
# dat <- melt(tmp[, c(1, 3, 6, 12)], id.vars = 1:2) 
# 
# # label 
# dat$lab <- ifelse(dat$rsid == rsidI, dat$rsid, NA)
# dat$size <- ifelse(dat$rsid == rsidI, "1", "0.5")
# 
# levels(dat$variable) <- c(protI, "AAM")
# 
# # plot
# ggplot(dat, aes(position, -log10(value), color = lab, size = size)) +
#   geom_vline(xintercept = posI) +
#   geom_point() +
#   facet_wrap(~ variable, nrow = 2) +
#   ggtitle(paste0("RSID of Interest = ", rsidI)) +
#   xlab("Chromosomal Position") +
#   theme_bw() +
#   theme(legend.position = "none") +
#   scale_size_manual(values = c(0.5, 2))

# COLOC ANALYSIS

# 2 dataset for analyse Prot - ANM
D1_expo <- list() 
D2_outc <- list()

# expo
D1_expo$type <- "quant"
D1_expo$snp <- tmp$rsid
D1_expo$position <- tmp$position
D1_expo$beta <- tmp$expo_eff
D1_expo$varbeta <- tmp$expo_varbeta
D1_expo$sdY <- 1 # (prot level are scale)

check_dataset(D1_expo)

# expo
D2_outc$type <- "quant"
D2_outc$snp <- tmp$rsid
D2_outc$position <- tmp$position
D2_outc$beta <- tmp$outc_eff
D2_outc$varbeta <- tmp$outc_varbeta
D2_outc$sdY <- "XXX" # (FROM GWAS STUDY : 1.3 for AAM and 4 for ANM)

# check_dataset(D2_outc)
# problem, some varbeta == 0
# replace by mean
# D2_outc$varbeta[D2_outc$varbeta == 0] <- mean(D2_outc$varbeta)

check_dataset(D2_outc)

# Run coloc analysis 
res <- coloc.abf(dataset1 = D1_expo, dataset2 = D2_outc)


# save the result 
tab_res <- rbind(tab_res, c(studyI, protI, rsidI, as.vector(res$summary)))
colnames(tab_res) <- c("Study", "Gene", "Rsid", "nSNP", "H0", "H1", "H2", "H3", "H4")

tab_res$H0 <- as.numeric(tab_res$H0)
tab_res$H1 <- as.numeric(tab_res$H1)
tab_res$H2 <- as.numeric(tab_res$H2)
tab_res$H3 <- as.numeric(tab_res$H3)
tab_res$H4 <- as.numeric(tab_res$H4)

# saveRDS(tab_res, "data/coloc....")

# Susie PLUG-IN -------------------
# need LD matrix for the set of tested RSID 
# we estimated the LD matrix from 1000 genomes using plink
ld <- "read ld files"

# calculate Z-score for exposure data and outcome data
z_expo <- expo_beta / expo_se
z_outc <- outc_beta / outc_se

# calcul of lambda
# concordance between the LD matrix and the z-scores
lambda_expo <- estimate_s_rss(z_expo, ld, n = N_expo)
lambda_outc <- estimate_s_rss(z_outc, ld, n = N_outc)

# use susis_rss same as runsusie (runsusie add some names)
fit_expo <- susie_rss(z_expo, ld, n = N_expo)
fit_outc <- susie_rss(z_outc, ld, n = N_outc)

# coloc susie
res <- coloc.susie(fit_expo, fit_outc)

# saveRDS(res, "data/coloc_susie....")