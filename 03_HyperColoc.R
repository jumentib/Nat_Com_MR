###### Generic Script used to do all of our HyperColocalization analysis ######

# Load Library 

library(hyprcoloc)
library(tidyverse)
library(reshape2)

# For each protein that was found significant in our two MR analyzes (AAM and ANM) 
# we carried out a hypercolocalization analysis

# Load GWAS
PROT <- readRDS("PROTEIN LEVEL GWAS")

AAM <- readRDS("AAM GWAS")

ANM <- readRDS("ANM GWAS")

# rsid of Interest
rsidI <- "rsXXX"
pos <- "position of rsidI"
info_rsidI <- PROT[PROT$POSI == pos, ]

# VERIFY position AND rsidI in GWAS PROT
info_rsidI

# 500 Mb < rsidI > 500 Mb
pos_min <- pos - 500000
pos_max <- pos + 500000

# cut GWAS
PROT_cut <- PROT[PROT$POSI >= pos_min & PROT$POSI <= pos_max, ]

# merge the 3 GWAS (BY POSITION rather than RSID)
tmp <- merge.data.frame(x = tmp_PROT, y = AAM)
tmp <- merge.data.frame(x = tmp, y = ANM)

# data manage for HyprColoc package
betas <- tmp[, c("beta_prot", "beta_aam", "beta_anm")]
ses <- tmp[, c("se_prot", "se_aam", "se_anm")]
traits <- c("PROT", "AAM", "ANM")
rs <- tmp$rsid 

# names
colnames(betas) <- traits
colnames(ses) <- traits
row.names(betas) <- rs
row.names(ses) <- rs

# as.matrix
betas <- as.matrix(betas)
ses <- as.matrix(ses)

# HyperColocalization analysis  
res <- hyprcoloc(effect.est = betas, 
                 effect.se = ses,
                 trait.names = traits, 
                 snp.id = rs, 
                 snpscores = T)

# saveRDS(tab_res, "data/hypercoloc....")