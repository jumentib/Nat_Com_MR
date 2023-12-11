###### Generic Script used to do all of our MR ######

# Load required packages
library(MRInstruments)
library(TwoSampleMR)
library(MRPRESSO)
library(LDlinkR)

# Set working directory to where your files are located
setwd("/path")
# Read the outcome file --> AAM or ANM
# See our Manuscript, to have the link to download the GWAS
outcome <- read.table("outcome", sep="\t", header=TRUE)

# Define the column names
colnames <- c("SNP", "Outcome", "Exposure", "Study", "Method", "Number.SNPs", "Beta", "SE", "Pval")

# Create an empty data frame with the desired column names
results <- data.frame(matrix(ncol=length(colnames), nrow=0))
colnames(results) <- colnames

# Read the exposure files
exposure_files <- c("file1","file2","...")

# Loop through each file
for (exposure_file in exposure_files) {
  
  exposure_list <- readLines(exposure_file)  # Read the file as a list of lines 
  
  # Group the exposure data by Protein
  exposure_data_list <- split(exposure_list[-1], sapply(exposure_list[-1], function(x) strsplit(x, "\t")[[1]][1]))
  
  # Loop through each protein group
  for (protein in names(exposure_data_list)) {
    
    # Extract the exposure data for the current protein group
    exposure_data_lines <- exposure_data_list[[protein]]
    
    # Create a data frame with the exposure data
    exposure_data <- data.frame(
      Protein = protein,
      SNP = sapply(exposure_data_lines, function(x) strsplit(x, "\t")[[1]][2]),
      Beta = as.numeric(sapply(exposure_data_lines, function(x) strsplit(x, "\t")[[1]][10])),
      SE = as.numeric(sapply(exposure_data_lines, function(x) strsplit(x, "\t")[[1]][11])),
      Pos = sapply(exposure_data_lines, function(x) strsplit(x, "\t")[[1]][5]),
      Allele.1 = sapply(exposure_data_lines, function(x) strsplit(x, "\t")[[1]][6]),
      Allele.2 = sapply(exposure_data_lines, function(x) strsplit(x, "\t")[[1]][7]),
      EAF = as.numeric(sapply(exposure_data_lines, function(x) strsplit(x, "\t")[[1]][8])),
      Pval = as.numeric(sapply(exposure_data_lines, function(x) strsplit(x, "\t")[[1]][12]))
    )
    
    # Rename the exposure data columns to match the format required by TwoSampleMR
    exposure_data$beta <- exposure_data$Beta
    exposure_data$se <- exposure_data$SE
    exposure_data$effect_allele <- exposure_data$Allele.1
    exposure_data$other_allele <- exposure_data$Allele.2
    exposure_data$eaf <- exposure_data$EAF
    exposure_data$id.exposure <- exposure_data$SNP
    exposure_data$pval <- exposure_data$Pval
    exposure_data$pval.exposure <- exposure_data$pval
    
    # remove commas and quotes from Pos column and new the column position
    exposure_data$position <- gsub("[\",]", "", exposure_data$Pos)
    
    # convert position column to numeric
    exposure_data$position <- as.numeric(exposure_data$position)
    
    # Prepare the exposure data for harmonisation
    exp_data <- format_data(exposure_data, snp_col = "SNP", type = "exposure")
    exp_data$id.exposure <- exposure_data$Protein
    
    write.csv(exp_data, paste0(protein,".tsv"), quote = FALSE, row.names = FALSE)
    
    if (nrow(exp_data) >= 2) {
      # Clump the exposure data for harmonisation
      new_exp_dat<-clump_data(exp_data)
    } else {
      # If there is only one row, skip the clumping step
      new_exp_dat <- exp_data
    }
    
    
    # Rename the outcome data columns to match the format required by TwoSampleMR
    outcome$SNP <- outcome$variant_id
    outcome$pval <- outcome$p_value
    outcome$se <- outcome$standard_error
    outcome$eaf <- outcome$effect_allele_frequency
    
    # Prepare the outcome data for harmonisation
    out_data <- format_data(outcome, snp_col = "SNP", type = "outcome")
    out_data$id.outcome <- "Trait"
    
    # Harmonise the exposure and outcome data
    dat <- harmonise_data(exposure_dat = new_exp_dat, outcome_dat = out_data)
    
    # Perform MR analysis
    res <- mr(dat)
    
    # If you want to perform Steiger Test
    
    # get R2 for exposure
    # r_exp <- get_r_from_bsen(b = beta_expo, 
    #                          se = se_expo, 
    #                          n = N_expo)
    # 
    # # get R2 for outcome
    # r_out <- get_r_from_bsen(b = beta_out, 
    #                          se = se_out, 
    #                          n = N_out)
    # 
    # # steiger
    # res <- mr_steiger(p_exp = p_val_expo, 
    #                   p_out = p_val_out, 
    #                   n_exp = N_expo, 
    #                   n_out = N_out, 
    #                   r_exp = r_exp,
    #                   r_out = r_out)
    # 
    # # Then save the result
    
    if (length(res) == 0)  next # if res is empty skip to next exposure file
    
    new_res <- data.frame(
      SNP=dat$SNP, 
      Outcome=res$id.outcome, 
      Exposure=res$id.exposure, 
      Study=exposure_file, 
      Method=res$method, 
      Number.SNPs=res$nsnp, 
      Beta=res$b, 
      SE=res$se, 
      Pval=res$pval, 
      stringsAsFactors=FALSE)
    
    # Append the new results to the existing data frame
    results <- rbind(results, new_res)
    
  }
  
}


# Write the results to a single file with the all exposures
# saveRDS(results, "data/MR....")
