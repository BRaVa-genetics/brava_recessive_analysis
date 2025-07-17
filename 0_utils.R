# R-utils
library(dplyr)
library(data.table)
library(bravastring)
# requires 'bravastring' pacakge (see frhl/github)
# devtools::load_all("/Users/flassen/Projects/28_bravastring/bravastring")

get_bbj_dir <- function(use_reskat_path=FALSE){
  if (use_reskat_path){
    relative_path <- "data/sumstats_reskat/biobankjapan"
  } else {
    relative_path <- "data/sumstats/biobankjapan"
  }
  return(normalizePath(relative_path, winslash = "/", mustWork = TRUE))
}

get_ukb_dir <- function(use_reskat_path=FALSE){
  if (use_reskat_path){
    relative_path <- "data/sumstats_reskat/uk-biobank"
  } else {
    relative_path <- "data/sumstats/uk-biobank"
  }
  return(normalizePath(relative_path, winslash = "/", mustWork = TRUE))
}

get_gnh_dir <- function(use_reskat_path=FALSE){
  if (use_reskat_path){
    relative_path <- "data/sumstats_reskat/genes-and-health"
  } else {
    relative_path <- "data/sumstats/genes-and-health"
  }
  return(normalizePath(relative_path, winslash = "/", mustWork = TRUE))
}

get_gel_dir <- function(use_reskat_path=FALSE){
  if (use_reskat_path){
    relative_path <- "data/sumstats_reskat/genomics-england"
  } else {
    relative_path <- "data/sumstats/genomics-england"
  }
  return(normalizePath(relative_path, winslash = "/", mustWork = TRUE))
}

get_biome_dir <- function(use_reskat_path=FALSE){
  if (use_reskat_path){
    relative_path <- "data/sumstats_reskat/biome_biobank"
  } else {
    relative_path <- "data/sumstats/biome-biobank"
  }
  return(normalizePath(relative_path, winslash = "/", mustWork = TRUE))
}

get_aou_dir <- function(use_reskat_path=FALSE){
  if (use_reskat_path){
    relative_path <- "data/sumstats_reskat/all-of-us"
  } else {
    relative_path <- "data/sumstats/all-of-us"
  }
  return(normalizePath(relative_path, winslash = "/", mustWork = TRUE))
}

# get a data.table of association paths
# as well as meta data. Note, that this
# requires the package 'bravastring' from
# #devtools::install_github("frhl/bravastring")
get_biobank_files_dt <- function(biobank, use_reskat=FALSE){
  require(bravastring)
  # biobank <- toupper(biobank)
  stopifnot(biobank %in% get_biobanks())
  if (biobank == "UKB") dir <- get_ukb_dir(use_reskat)
  if (biobank == "GNH") dir <- get_gnh_dir(use_reskat)
  if (biobank == "BBJ") dir <- get_bbj_dir(use_reskat)
  if (biobank == "GEL") dir <- get_gel_dir(use_reskat)
  if (biobank == "BioMe") dir <- get_biome_dir(use_reskat)
  if (biobank == "AOU") dir <- get_aou_dir(use_reskat)
  files <- list.files(dir, full.names=TRUE)

  # setup data.table with out information
  d <- data.table(
	  biobank=str_extract_biobank(basename(files)),
	  ancestry=str_extract_ancestry(basename(files)),
	  sex=str_extract_sex(basename(files)),
	  trait=str_extract_brava_trait(basename(files)),
	  encoding=str_extract_encoding(basename(files)),
	  annotation=str_extract_annotation(basename(files)),
	  cases=str_extract_cases(basename(files)),
	  controls=str_extract_controls(basename(files)),
	  size=file.size(files),
	  full_path=files
  )

  # return d
  return(d)
}

map_transcript_to_gene_id <- function(df, gene_names, ID_tag) {
  # header <- as.character(data[nrow(data), ]); data <- data[-nrow(data), ]
  # df$transcript_id <- df$Region
  df <- df %>% 
    select( everything(), transcript_id=all_of(ID_tag) ) %>%
    left_join( gene_names, by='transcript_id' ) %>%
    select( everything(), ID='gene_id', -hgnc_symbol, -transcript_id )
  return(df)
}

filter_out_genes_based_on_burden <- function( df, path, biobank, annotation, ac_cutoff=10 ) {
  # stick to genes with MORE than ac_cutoff/2 carriers, based on the corresponding burden analysis
  # this way we avoid genes not processed by the burden analysis (which could have been skipped due to even fewer carriers)
  df_ref <- fread(path)
  if ('CHROM' %in% colnames(df_ref)) {
    # this is from regenie
    good_genes <- df_ref %>% select( ID,N, A1FREQ ) %>% filter( A1FREQ * N * 2>=ac_cutoff ) %>% .$ID
  } else {
    # this is from saige
    good_genes <- df_ref %>% select( MarkerID, AC_Allele2 ) %>% filter( AC_Allele2>=ac_cutoff ) %>% .$MarkerID
  }
  # cat("Bad genes:",length(bad_genes), "\n")
  df <- df %>% mutate( 
    Pvalue = case_when( 
      GROUP == annotation & !(ID %in% good_genes) ~ NaN, 
      TRUE ~ Pvalue ) )
  return(df)
}


get_brava_analysis_software <- function(path){
  header <- system(paste("(zcat < ", path, "| head -n1) 2>&1 | grep -v 'Broken pipe'"), intern = TRUE)
  if (grepl(header, pattern="GENPOS")){
    return("REGENIE")
  } else if (grepl(header, patter="MarkerID")){
    return("SAIGE")
  } else {
    return(NA)
  }
}


fread_saige <- function(f) {
  
  d <- fread(f)
  
  is_binary <- !("N" %in% colnames(d))
  
  if (is_binary) {
    d[, (setdiff(c("AF_case", "AF_ctrl"), names(d))) := NA]
    d <- d[, .(
      CHROM = CHR,
      ID = MarkerID,
      AF = AF_Allele2,
      AC = AC_Allele2,
      AF_CASE = AF_case,
      AF_CTRL = AF_ctrl,
      N = N_case + N_ctrl,
      N_CASE = N_case,
      N_CTRL = N_ctrl,
      BETA,
      SE,
      CHISQ = (BETA/SE)^2,
      P = p.value,
      BINARY = TRUE,
      TOOL = "SAIGE",
      CONVERGE = ifelse(p.value.NA<0.01 | p.value<0.01, Is.SPA, TRUE)
    )]
  } else {
    d <- d[, .(
      CHROM = CHR,
      ID = MarkerID,
      AF = AF_Allele2,
      AC = AC_Allele2,
      AF_CASE = NA,
      AF_CTRL = NA,
      N = N,
      N_CASE = NA,
      N_CTRL = NA,
      BETA,
      SE,
      CHISQ = (BETA/SE)^2,
      P = p.value,
      BINARY = FALSE,
      TOOL = "SAIGE",
      CONVERGE = TRUE
    )]
  }
  
  return(d)
}

fread_saige_reskat <- function(f) {
  
  d <- fread(f)
  is_binary <- ("MAC_case" %in% colnames(d))
  
  if (is_binary) {
    d <- d[, .(
      ID = Region,
      AF = NA,
      AC = NA,
      AF_CASE = NA,
      AF_CTRL = NA,
      N = NA,
      N_CASE = NA,
      N_CTRL = NA,
      Pvalue,
      Pvalue_Burden,
      Pvalue_SKAT,
      BETA_Burden,
      SE_Burden,
      MAC,
      MAC_case,
      MAC_control,
      Number_rare,
      Number_ultra_rare,
      BINARY = TRUE,
      TOOL = "SAIGE",
      CONVERGE = NA,
      GROUP = Group
    )]
  } else {
    d <- d[, .(
      ID = Region,
      AF = NA,
      AC = NA,
      AF_CASE = NA,
      AF_CTRL = NA,
      N = NA,
      N_CASE = NA,
      N_CTRL = NA,
      Pvalue,
      Pvalue_Burden,
      Pvalue_SKAT,
      BETA_Burden,
      SE_Burden,
      MAC,
      MAC_case = NA,
      MAC_control = NA,
      Number_rare,
      Number_ultra_rare,
      BINARY = FALSE,
      TOOL = "SAIGE",
      CONVERGE = TRUE,
      GROUP = Group
    )]
  }
  
  return(d)
}


fread_regenie <- function(f) {
  
  d <- fread(f)
  
  # this is needed beacause bbj does not have case column
  # but gnh do have case columns
  has_a1freq_col  <- "A1FREQ_CASES" %in% colnames(d)
  has_cases_col <- "N_CASES" %in% colnames(d)
  
  d <- d[, .(
    CHROM = NA,
    ID,
    AF = A1FREQ,
    AC = round(A1FREQ * N * 2),
    AF_CASE = ifelse(has_a1freq_col, A1FREQ_CASES, NA),
    AF_CTRL = ifelse(has_a1freq_col, A1FREQ_CONTROLS, NA),
    N,
    N_CASE = ifelse(has_cases_col, N_CASES, NA),
    N_CTRL = ifelse(has_cases_col, N_CONTROLS, NA),
    BETA,
    SE,
    CHISQ,
    P = 10^(-LOG10P),
    BINARY = has_a1freq_col,
    TOOL = "REGENIE",
    CONVERGE = TRUE
  )]
  
  return(d)
}


run_stouffer_dplyr <- function(
    grouped_dt, n_eff_name, weighted_Z_name,
    input_pvalues, input_study, output_meta_pvalue,
    output_num_cohorts = "num_cohorts", output_min_pvalue = "min_pvalue", output_min_p_id = "min_p_id",
    two_tail = FALSE, input_beta = NULL
) {
  if (two_tail) {
    grouped_dt <- grouped_dt %>%
      mutate("{input_pvalues}" := .data[[input_pvalues]]/2)
  } else {
    input_beta <- "beta_dummy"
    grouped_dt <- grouped_dt %>% mutate("{input_beta}" := 1)
  }
  
  stopifnot(is.numeric(grouped_dt[[input_pvalues]]))
  stopifnot(is.numeric(grouped_dt[[n_eff_name]]))
  
  stopifnot(!any(is.na(grouped_dt[[input_pvalues]])))
  stopifnot(!any(is.na(grouped_dt[[n_eff_name]])))
  
  result <- grouped_dt %>%
    mutate(
      weighted_Z_numerator = (
        sqrt(.data[[n_eff_name]]) *
          (-qnorm(.data[[input_pvalues]])) *
          sign(.data[[input_beta]])
      )
    )
  
  if (two_tail) {
    result <- result %>%
      summarise(
        "{weighted_Z_name}" := sum(weighted_Z_numerator) /
          sqrt(sum(.data[[n_eff_name]])),
        "{output_meta_pvalue}" := 2 * pnorm(abs(.data[[weighted_Z_name]]), lower.tail=FALSE),
        "{output_num_cohorts}" := sum(!is.na(.data[[input_pvalues]])),
        "{output_min_pvalue}" := min(.data[[input_pvalues]], na.rm=TRUE),
        "{output_min_p_id}" := .data[[input_study]][which.min(.data[[input_pvalues]])]  # Study with min p-value
      )
  } else {
    result <- result %>%
      summarise(
        "{weighted_Z_name}" := sum(weighted_Z_numerator) /
          sqrt(sum(.data[[n_eff_name]])),
        "{output_meta_pvalue}" := pnorm(.data[[weighted_Z_name]], lower.tail=FALSE),
        "{output_num_cohorts}" := sum(!is.na(.data[[input_pvalues]])),
        "{output_min_pvalue}" := min(.data[[input_pvalues]], na.rm=TRUE),
        "{output_min_p_id}" := .data[[input_study]][which.min(.data[[input_pvalues]])]  # Study with min p-value
      )
  }
  return(result)
}



run_inv_var_dplyr <- function(
    grouped_dt, input_beta_name, input_se_name,
    output_beta_meta, output_se_meta,
    output_meta_pvalue
) {
  dt <- grouped_dt %>%
    mutate(weight = 1/(.data[[input_se_name]]**2))
  
  dt <- dt %>% mutate(effs_inv_var = .data[[input_beta_name]] * weight)
  dt <- dt %>%
    summarise(
      "{output_beta_meta}" := sum(effs_inv_var) / sum(weight),
      "{output_se_meta}" := sqrt(1/sum(weight))) %>%
    mutate("{output_meta_pvalue}" := 2 * pnorm(
      abs(.data[[output_beta_meta]] / .data[[output_se_meta]]), lower=FALSE))
  
  # Deal with the edge cases
  # Send the pvalues of the data with infinite standard errors to 1, and the
  # pvalues of the data with standard errors of 0 to 0.
  dt$P[which(is.na(dt$P) & (dt$SE == Inf))] <- 1
  dt$P[which(is.na(dt$P) & (dt$SE == 0))] <- 0
  
  return(dt)
}


run_heterogeneity <- function(
    grouped_dt, n_eff_name, input_beta, output_meta_beta)
{
  grouped_dt <- grouped_dt %>%
    mutate(
      weights = sqrt(.data[[n_eff_name]]),
      beta = .data[[input_beta]]
    )
  summary_dt <- grouped_dt %>%
    summarise(
      sum_weights = sum(weights),
      sum_betas = sum(weights * beta)
    ) %>% mutate("{output_meta_beta}" := sum_betas/sum_weights)
  return(
    merge(grouped_dt, summary_dt) %>%
      mutate(
        deviation = weights * (beta - .data[[output_meta_beta]])^2
      ) %>% group_by(across(as.character(groups(grouped_dt)))) %>%
      summarise(sum_deviation = sum(deviation), n_studies=n()) %>%
      filter(n_studies > 1) %>% mutate(p_het = pchisq(sum_deviation, n_studies-1, lower.tail=FALSE))
  )
}

weights <- function(dt, is_inv_var=TRUE,
                    n_eff_name=NULL, se_name=NULL) {
  if (is_inv_var) {
    return(dt %>% mutate(weights = 1/(.data[[se_name]]^2)))
  } else {
    return(dt %>% mutate(weights = sqrt(.data[[n_eff_name]])))
  }
}

run_heterogeneity_test <- function(
    grouped_dt, input_beta, output_meta_beta) {
  grouped_dt %>%
    mutate(weighted_effects = weights * .data[[input_beta]]) %>%
    summarise(
      df = n()-1,
      sum_weights = sum(weights),
      "{output_meta_beta}" := sum(weighted_effects) / sum_weights,
      chisq_het = sum(
        weights * (.data[[input_beta]] - .data[[output_meta_beta]])^2),
      Pvalue_het = pchisq(chisq_het, df, lower.tail=FALSE)
    )
}

run_cauchy_ungrouped <- function(
    dt, Cauchy_stat_name, pvalue_columns, output_meta_pvalue)
{
  if (!all(pvalue_columns %in% colnames(dt))) {
    stop("Some p-value columns specified are not present in the data.")
  }
  return(
    dt %>%
      rowwise() %>%  # Operate row-by-row
      mutate(
        "{Cauchy_stat_name}" := cauchy_combination(
          na.omit(c_across(all_of(pvalue_columns)))
        ),
        number_of_pvals := sum(!is.na(c_across(all_of(pvalue_columns))))
      ) %>%
      ungroup() %>%
      mutate(
        "{output_meta_pvalue}" := ifelse(
          .data[[Cauchy_stat_name]] > 1e+15,
          (1 / .data[[Cauchy_stat_name]]) / pi,
          pcauchy(.data[[Cauchy_stat_name]], lower.tail=FALSE)
        )
      ) %>% mutate(
        "{output_meta_pvalue}" := ifelse(
          .data[[output_meta_pvalue]] > (1 - 1e-10),
          (1 - 1/number_of_pvals), .data[[output_meta_pvalue]]
        )
      )
  )
}

run_cauchy <- function(
    grouped_dt, Cauchy_stat_name, input_pvalues, output_meta_pvalue)
{
  return(
    grouped_dt %>%
      summarise(
        "{Cauchy_stat_name}" := cauchy_combination(
          .data[[input_pvalues]],
          weights=
        ),
        number_of_pvals := n()
      ) %>% mutate(
        "{output_meta_pvalue}" := ifelse(
          .data[[Cauchy_stat_name]] > 1e+15,
          (1 / .data[[Cauchy_stat_name]]) / pi,
          pcauchy(.data[[Cauchy_stat_name]], lower.tail=FALSE)
        )
      ) %>% mutate(
        "{output_meta_pvalue}" := ifelse(
          .data[[output_meta_pvalue]] > (1 - 1e-10),
          (1 - 1/number_of_pvals), .data[[output_meta_pvalue]]
        )
      )
  )
}

cauchy_combination <- function(p_values, weights=NULL)
{

  is.zero <- sum(p_values == 0) >= 1
  is.one <- sum(p_values > (1 - 1e-14)) >= 1
  
  if (is.zero) {
    return(0)
  }
  
  if (is.one) {
    p <- min(p_values) * length(p_values)
    if (p > 1) {
      return(-Inf)
    } else {
      return(qcauchy(p, lower.tail=FALSE))
    }
    # return(min(1, (min(p_values)) * (length(p_values)))) # redundant
  }
  
  if (is.null(weights)) {
    weights <- rep(1 / length(p_values), length(p_values))
  } else if (length(weights) != length(p_values)) {
    stop("The length of weights should be the same as that of the p-values!")
  } else if (sum(weights < 0) > 0) {
    stop("All the weights must be positive!")
  } else {
    weights <- weights / sum(weights)
  }
  
  is_small <- (p_values < 1e-16)
  if (sum(is_small) == 0) {
    cct_stat <- sum(weights * tan((0.5 - p_values) * pi))
  } else {
    cct_stat <- sum((weights[is_small] / p_values[is_small]) / pi)
    cct_stat <- cct_stat +
      sum(weights[!is_small] * tan((0.5 - p_values[!is_small]) * pi))
  }
  
  return(cct_stat)
}
