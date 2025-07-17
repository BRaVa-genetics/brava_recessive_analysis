# Script to perform meta-analysis of BRaVa-recessive
# INPUT: file with all results across cohorts and annotations, output from 04_combine_sumstats.R
# OUTPUT: meta-analysis and heterogeneity results

# dev by FL (Sept 2024), refactored by GK (March 2025)

library(data.table)
library(dplyr)
library(bravastring)
source("scripts/00_utils.R")

# setwd("~/Projects/02_brava_recessive/brava_recessive_summary_statistics/")
# devtools::load_all("/Users/flassen/Projects/03_bravastring/bravastring")
# source("~/Projects/02_brava_recessive/brava_recessive_summary_statistics/scripts/00_utils.R")

# out_path="~/Dropbox/BRaVa/meta/"
out_path='derived/'
out_hetr=paste(out_path,"BRaVa_recessive.",format(Sys.Date(), "%d%m%y"),".heterogeneity.txt.gz", sep='') # withoutnonsyn.
out_meta=paste(out_path,"BRaVa_recessive.",format(Sys.Date(), "%d%m%y"),".meta_cauchy.txt.gz", sep='')
# read sumstat file
#d_full <- fread("derived/BRaVa_recessive.sumstats.combined.UKB_AOU_GEL_GNH_BBJ_BioMe.110425.txt.gz")
d_full[ ,ANALYSIS_ID := paste0(BIOBANK,"_",TRAIT, "_", ANCESTRY, "_", SEX, "_",ENCODING,"_",ANNOTATION)]
d_full[ ,BIOBANK_ANCESTRY := paste0(BIOBANK, "_", ANCESTRY)]
stopifnot(sum(duplicated(d_full))==0)

# some checks
stopifnot(any(!is.na(d_full$N_EFF)))
stopifnot(any(!is.na(d_full$P)))
cts_traits <- c("ALT", "AlcCons", "AST", "BMI", "HDLC", "Height", "LDLC", "TChol", "TG", "WHRBMI", "CRP")
# drop (any) withdrawn phenotypes
d_full <- d_full[!d_full$TRAIT %in% c("AlcCons","WHRBMI"),]

traits <- d_full %>% pull(TRAIT) %>% unique()
biobanks <- d_full %>% pull(BIOBANK) %>% unique()

cat("Total number of results:", dim(d_full)[1], '\n')
cat("Number of traits in analysis:",length(traits),'\n')
cat("Biobanks under consideration:", biobanks, '\n')
cat("Annotation masks under consideration:",d_full %>% pull(ANNOTATION) %>% unique(),'\n')

# prepare grouped file for meta-analysis
grouped_dt <- d_full %>%
  group_by(SEX, TRAIT, ID, ANNOTATION, ENCODING)
# note: we could now delete `d_full` to save memory

# ensure that these columns are numerics!
grouped_dt$BETA <- as.numeric(grouped_dt$BETA)
grouped_dt$SE <- as.numeric(grouped_dt$SE)
grouped_dt$P <- as.numeric(grouped_dt$P)
grouped_dt$N_EFF <- as.numeric(grouped_dt$N_EFF)

# which ones are missing
stopifnot(nrow(grouped_dt[is.na(grouped_dt$P),])==0)
stopifnot(nrow(grouped_dt[is.na(grouped_dt$N_EFF),])==0)
#stopifnot(nrow(stopifnot(nrow(grouped_dt[is.na(grouped_dt$BETA),])==0)
#stopifnot(nrow(grouped_dt[is.na(grouped_dt$SE),])==0)

# run stouffer #
################
meta_stouffer <- run_stouffer_dplyr(
  grouped_dt = grouped_dt,
  n_eff_name = "N_EFF",
  weighted_Z_name = "Z",
  input_pvalues = "P",
  input_beta = "BETA",
  two_tail = TRUE,
  output_meta_pvalue = "P_meta",
  input_study = "BIOBANK",
  output_num_cohorts = "number_of_pvals"
)
cat("Done with Stouffer meta-analysis.\n")
# note: `run_stouffer_dplyr` returns columns with number of input pvals and min-pvalue (and min-p-id) which we should keep
# note: we need to change the name of the output to P (we needed "P_meta" above to correctly count the number of input p-values)
setnames(meta_stouffer, "P_meta", "P")

# truncate ultra-low P-values that are zero:
meta_stouffer_for_cauchy <- meta_stouffer
min_p <-  min(meta_stouffer_for_cauchy$P[meta_stouffer_for_cauchy$P!=0.000000e+00])
setDT(meta_stouffer_for_cauchy)
meta_stouffer_for_cauchy[P==0.000000e+00, P := min_p]

####################################
# Cauchy combination across groups #
####################################
meta_stouffer_cauchy <- meta_stouffer_for_cauchy %>%
  filter(ANNOTATION != "synonymous") %>%
  group_by(SEX, TRAIT, ID, ENCODING) %>% 
  run_cauchy("CCT", "P", "P")
cat("Done with Cauchy combination.\n")

# prepare for combination by adding missing columns
meta_stouffer_cauchy$ANNOTATION <- "cauchy"
meta_stouffer_cauchy$min_p_id <- NA
meta_stouffer_cauchy$min_pvalue <- NA
# Cauchy doesnt give Z and Stouffer doesnt give CCT:
meta_stouffer_cauchy$Z <- NA
meta_stouffer$CCT <- NA
# meta_stouffer$number_of_pvals <- 1
# meta_stouffer$num_cohorts <- NULL
# meta_stouffer$min_pvalue <- NULL

# merge and re-order columns:
meta_stouffer_cauchy <- 
  rbind(meta_stouffer,
        meta_stouffer_cauchy, fill=TRUE) %>%
  arrange(ENCODING, SEX, TRAIT, ID, ANNOTATION)
# add some columns:
meta_stouffer_cauchy$BIOBANK <- "META"
meta_stouffer_cauchy$ANCESTRY<- "META"
# write output
fwrite(meta_stouffer_cauchy, out_meta, sep="\t", na="NA", quote=FALSE)

#################################
# get meta P value heterogenity #
#################################
meta_hetereogeneity <- data.table(run_heterogeneity_test(
  calc_weights_meta(grouped_dt, FALSE, se_name="SE", n_eff_name="N_EFF"),
  input_beta="BETA", output_meta_beta="BETA_HET"
),  key=c("SEX", "TRAIT", "ID", "ANNOTATION", "ENCODING"))
cat("Done with heterogeneity testing.\n")
# write output
fwrite(meta_hetereogeneity, out_hetr, sep="\t", na="NA", quote=FALSE)
