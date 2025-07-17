# BRaVa Recessive Analysis Pipeline
#### Scripts needed for BRaVa's recessive manuscript.

This repository contains a series of scripts needed for BRaVa's recessive project, described in *Lassen F. H., Kalantzis G. et al (2025)*. The pipeline identifies compound heterozygous and homozygous gene knockouts across multiple biobanks, performs federated association testing, and then meta-analyzes the summary statistics. Please note that these script are meant to be the backbone of the analysis and some biobank-specific scripts might be missing.

The biobanks included are UK Biobank (UKB), All of Us Research Program (AOU), 100,000 Genomes Project (GEL or 100kGP), Genes & Health (GNH or G&H), BioMe Biobank (BioMe), and  Biobank Japan (BBJ).  We consider the following variant consequences (burdens):
1. pLoF: Predicted loss-of-function variants
2. damaging_missense
3. nonsynonymous: all nonsynonymous variants
4. synonymous (control)

### Pipeline Overview
0. Statistical phasing with SHAPEIT5. For this step check our [snakemake pipeline](https://github.com/BRaVa-genetics/snakemake_pipeline_for_phasing).
1. Compound heterozygote calling 
2. Gene knockout analysis
3. Association testing with SAIGE
4. QC of statistics combination
5. Meta-analysis (Stouffer and Cauchy combination)

### More details for each step

#### 1a. 011_run_call_chets.sh
> This is a wrapper for another module, [call_chets](https://github.com/BRaVa-genetics/biallelic), to identify compound heterozygous variants and create VCF files for different consequence categories, per chromosome.
* Input: Phased genotypes (BCF), variant annotations
* Output: VCF files with additive and recessive encodings for each chromosome
* Dependencies: phase data from another pipeline; annotation files (see Main BRaVa)
* Key Features:
    * Filters variants by posterior probability (PP > 90%) and allele frequency (AF < 5%)
Processes multiple consequence categories: pLoF, damaging_missense, synonymous, nonsynonymous.
    * Generates both additive (`add`) and recessive (`rec`) genetic encodings

#### 1b. 012_concat_vcf_and_make_bgen.sh
> Concatenates chromosome-specific VCF files and converts to BGEN format needed for more efficient association testing.
* Input: Per-chromosome VCF files from previous step
* Output: Genome-wide VCF and BGEN files
* Dependencies: bcftools, plink2, bgenix

#### 2a. 021_extract_gene_counts.py
> Tallies gene knockout counts by genotype and consequence across all chromosomes
* Input: Compound heterozygote calling results
* Output: Summary table of gene knockouts (homozygous, compound heterozygous, cis-heterozygous)

#### 2b. 022_gene_KO_analysis.R
> Comprehensive analysis of gene knockouts across biobanks and comparison with literature
* Input: output from previous step (across biobanks) as well as tables from other studies.
* Output: Table of gene knockouts (Sup. Table 4). It also prints summaries of counts used in the manuscript (e.g. Sup. Table 3).
* Key Features:
    * Combines biallelic counts from all biobanks (BBJ, BioMe, GNH, GEL, UKB, AOU)
    * Maps gene IDs to symbols using HGNC annotations
    * Compares findings with previous studies (Oddsson et al., Sun et al.)
    * Identifies novel gene knockouts unique to BRaVa

#### 3. 03_saige_step2_recessive.r
> Performs gene-based association testing using SAIGE. It assumes step-1 files are already existent. Since we have one score per gene, this is essentially the same as single-variant testing. Note that REGENIE can also be used here.
* Input: BGEN files, phenotype data, genetic relationship matrix (GRM), step-1 files
* Output: Association test results with p-values and effect sizes
* Features:
    * Uses SAIGE's mixed-model approach with Firth correction for rare variants

#### 4. 04_combine_sumstats.R
> Collates and standardizes summary statistics across all biobanks. Note that more preprocessing might be needed depending on output from individual biobanks.
* Input: Individual biobank association results
* Output: Combined summary statistics file
* Dependency: `bravastrings`, file with effective sample sizes
* Key Features:
    * Handles multiple biobanks and ancestry groups
    * Applies quality control filters (minimum case counts, allele counts)
    * Standardizes file formats across different analysis software (SAIGE, REGENIE)

#### 5. 05_meta_analysis.R
> Purpose: Performs meta-analysis with Stouffer's method across biobanks and annotations, followed by the Cauchy combination test. Also performs tests for effect heterogeneity using Cochran's Q statistic.
* Input: Combined summary statistics from step 04
* Output: Meta-analysis results and heterogeneity statistics
Methods:
Stouffer's method for combining p-values across biobanks
Cauchy combination test for combining across annotation categories
Heterogeneity testing 

#### Software Dependencies
* Tools: bcftools (v1.18+); plink2; bgenix; SAIGE; [call_chets](https://github.com/BRaVa-genetics/biallelic)
* R packages: dplyr, data.table, bravastring
* Python packages: pandas, numpy

### Notes
* Check out this [Google Doc](https://docs.google.com/document/d/1fdoGeZGyE7tsXurEyr4LcCJK6J_RwlR1q-Ib_XWdBYE/edit?tab=t.0) with more details on the analysis.