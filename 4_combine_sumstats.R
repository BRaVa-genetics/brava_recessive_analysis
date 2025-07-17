# Script to collate summary statistics for BRaVa recessive meta-analysis
# ensure that all files are there and they have the right format otherwise some results might be skipped!

# Note: this script assumes the presence of 
#  > derived/studies_to_drop.txt, with specific studies to skip (e.g. for small sample sizes or inflation)
#  > derived/pilot_brava_n_eff_using_googlesheets.txt, file with effective sample size per phenotype per cohort (see Main Brava)

# dev by FL (Sept 2024), refactored by GK (March 2025)

# devtools::load_all("/Users/flassen/Projects/03_bravastring/bravastring/")
setwd("FIXTHIS/project_BRaVa_Recessive/brava_recessive_summary_statistics")
source("00_utils.R")
library(bravastring)
library(data.table)

sex_specific <- c("BenCervUterNeo","CervCanc", "BreastCanc", "EFRMB", "FemInf", "MatHem")
cts_traits <- c("ALT", "AlcCons", "AST", "BMI", "HDLC", "Height", "LDLC", "TChol", "TG", "WHRBMI", "CRP")

# read which biobanks to process 
# args_biobank <- commandArgs(trailingOnly = TRUE)
args_biobank <- c('UKB','AOU','GEL','GNH','BioMe','BBJ')
# Check if at least one argument is provided
if (length(args_biobank) == 0) {
  stop("No study provided. Please specify at least one study as an argument.")
}
fout=paste("derived/BRaVa_recessive.sumstats.combined.",paste(args_biobank,collapse='_'),format(Sys.Date(), ".%d%m%y"),".txt.gz", sep='')
cat("Will save output to", fout, "\n")

# Dynamically define constants for the analysis
study_tables <- lapply(args_biobank, function(x) {
  get_biobank_files_dt(x, use_reskat = FALSE)
})
cat("Studies found:", length(study_tables), "\n")

# combine all the biobanks
# files_dt <- rbind(ukb, bbj, gnh, biome, aou, gel)
files_dt <- do.call(rbind, study_tables)
files_dt[ , ID := paste(biobank,ancestry,trait,sep="_")]
files_dt$ancestry <- toupper(files_dt$ancestry)
files_dt$description <- get_brava_description(files_dt$trait)

# only keep recessive + additive
files_dt <- files_dt[files_dt$encoding %in% c("additive","recessive")]

# deal with bug where regex grabs zero when there are no 
# cases controls in the title of the sumstat file
# TODO: check if needed!
# files_dt$cases[files_dt$cases=="0" & is.na(files_dt$controls)] <- NA 

# subset to pilot traits
pilot_traits <- fread(pilot_traits_path(), header=FALSE)$V1
files_dt <- files_dt[files_dt$trait %in% pilot_traits,]
cat("Files found with Burden results:", length(files_dt$ID), "\n")

# Drop BioMe traits (Did not converge / something went wrong during fitting)
to_drop <- fread("studies_to_drop.txt")
to_drop[ , ID_TO_DROP := paste(biobank,ancestry,trait,sep="_")]
# biome[ , ID := paste(biobank,ancestry,trait,sep="_")]
n_to_drop <- sum(files_dt$ID %in% to_drop$ID_TO_DROP)
cat("Results to remove:", n_to_drop, "\n")
files_dt <- files_dt[!files_dt$ID %in% to_drop$ID_TO_DROP,]

# deal with a few incorrectly stored file names
files_dt$sex[files_dt$sex=="both"] <- "ALL"

# we only want to use females for sex-specific traits, the remaining should be "ALL".
# i.e. we are not running female only analyses when sex=ALL.
files_dt <- files_dt[(files_dt$trait %in% sex_specific & files_dt$sex == "F") |  
                       (!(files_dt$trait %in% sex_specific) & files_dt$sex == "ALL")]

# remove dominance encoding for binary traits
binary_traits <- pilot_traits[!(pilot_traits %in% cts_traits)]
files_dt <- files_dt[!(files_dt$trait %in% binary_traits & files_dt$encoding == "dominance"),]
files_dt$is_binary <- !(files_dt$trait %in% cts_traits)

# extract cases/controls from google sheet document when needed
dt_cases <- fread("derived/pilot_brava_n_eff_using_googlesheets.txt")
dt_cases$biobank <- toupper(dt_cases$biobank)
dt_cases$ancestry <- toupper(dt_cases$ancestry)
dt_cases <- dt_cases[cases>0 & binary==TRUE,]

# map to cases when need be
dt_cases[ , ID := paste(biobank,ancestry,trait,sep="_")]
id_to_cases <- dict(dt_cases$ID, dt_cases$cases)
id_to_controls <- dict(dt_cases$ID, dt_cases$controls)

# remove any case/controls traits wiht less than 100 cases
dt_cases <- dt_cases[dt_cases$cases<100]
files_dt <- files_dt[!files_dt$ID %in% dt_cases$ID,]
cat("Number of files after dropping those with low case numbers:", length(files_dt$ID), "\n")

# check for duplicates
files_dt[ ,DUP_ID := paste0(biobank,"_",trait, "_", ancestry, "_", sex, "_",encoding,"_",annotation)]
(dups <- files_dt[duplicated(files_dt$DUP_ID)])
stopifnot(nrow(dups)==0)

# we use different functions to read
# NOTE, that GEL files uses plain ole fread
read_fun_list <- list(
  "SAIGE" = fread_saige,
  "REGENIE" = fread_regenie
)

# Process files in parallel if possible
dt_list <- lapply(1:nrow(files_dt), function(i) {

  if (i %% 20 == 0) { message(sprintf("Processing %d of %d", i, nrow(files_dt))) }
  
  # Extract metadata once
  file_info <- files_dt[i, ]

  # get the software to read 
  software <- get_brava_analysis_software(file_info$full_path)
  read_fun <- if (is.null(read_fun_list[[software]])) fread else read_fun_list[[software]]

  # Read the data
  d <- read_fun(file_info$full_path)
  
  # Convert numeric columns
  d$N_CASE <- as.numeric(d$N_CASE)
  d$N_CTRL <- as.numeric(d$N_CTRL)
  d$P <- as.numeric(d$P)
  
  # Add metadata columns
  metadata_cols <- c("biobank", "ancestry", "sex", "trait", "encoding", "annotation")
  for (col in metadata_cols) {
    d[[toupper(col)]] <- file_info[[col]]
  }
  
  # Remove NA P-values
  d <- d[!is.na(d$P), ]
  
  # Handle missing case/control counts for binary traits
  if (file_info$is_binary && any(is.na(d$N_CASE))) {
    cases <- as.numeric(id_to_cases[file_info$ID])
    controls <- as.numeric(id_to_controls[file_info$ID])
    
    if (is.na(cases) || cases <= 0) {
      stop(sprintf("Invalid case count for ID %s", file_info$ID))
    }
    
    d[, N_CASE := cases]
    d[, N_CTRL := controls]
  }
  
  # Filter by sample size
  if (all(d$BINARY)) {
    d <- d[d$N_CASE >= 100, ]
  } else {
    d <- d[d$N >= 100, ]
  }
  
  # Calculate effective sample size if not present
  if (!"N_EFF" %in% colnames(d)) {
    d[, N_EFF := ifelse(file_info$is_binary, 
                        4 / ((1 / N_CASE) + (1 / N_CTRL)), 
                        N)]
  }
  
  # Subset to AC >= 10 except for AoU,GEL
  if (!file_info$biobank %in% c("AOU","GEL")) {
    d <- d[d$AC >= 10, ]
  }
  
  # Standardize encoding values
  encoding_map <- c("recessive" = "R", "additive" = "A", "dominance" = "D")
  d$ENCODING <- encoding_map[d$ENCODING]
  
  # Keep only necessary columns
  cols_to_keep <- c("BIOBANK", "ANCESTRY", "SEX", "TRAIT", "ENCODING",
                    "ANNOTATION", "ID", "BETA", "SE", "P", "N_CASE",'N_CTRL', "N_EFF")
  d <- d[, ..cols_to_keep]
  #d$path <- file_info$full_path
  
  return(d)
})

# Remove NULL results (failed processing)
dt_list <- dt_list[!sapply(dt_list, is.null)]

# Combine all data tables (if needed)
combined_dt <- rbindlist(dt_list, fill = TRUE)
# combined_dt$DUP_ID <- NULL

fwrite(combined_dt , fout, sep="\t", na="NA", quote=FALSE)

# end-of-script
