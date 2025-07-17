library(SAIGE)
library(optparse)

option_list = list(
  make_option(c("-o", "--out"), type="character", default="saige_gene", help="prefix for any output"),
  make_option("--vcf", type="character", default="chr11", help="filepath for the VCF file"),
  make_option("--bgen", type="character", default="chr11", help="filepath for the BGEN file"),
  make_option("--model", type="character", help="prefix for rda and varRatio files"),
  make_option("--grm", type="character", help="full path to the GRM"),
  make_option("--sample", type="character", help="full path to a sample file"),
  make_option(c("-c", "--chr"), type="character", default="chr22", help="Chromosome to test for")
  # make_option("--wdir", type="character", default="./", help="working directory"),
  # make_option(c("-a", "--annot"), type="character", help="full path to file with annotation")
)

opt_parser = OptionParser(option_list=option_list)
opts = parse_args(opt_parser)

SPAGMMATtest(bgenFile = opts$bgen,
             bgenFileIndex = paste(opts$bgen, '.bgi', sep=''),
             subSampleFile = opts$sample,
             dosage_ZEROD_MAC_cutoff = 0.5,
             min_MAC = 2,
             min_MAF = 1e-4,
             GMMATmodelFile = paste(opts$model, ".rda", sep=''),
             varianceRatioFile = paste(opts$model, ".varianceRatio.txt", sep=''),
             sparseGRMFile = opts$grm,
             sparseGRMSampleIDFile = paste(opts$grm, ".sampleIDs.txt", sep=''),
             is_Firth_beta = TRUE,
             pCutoffforFirth = 0.1,
             is_output_moreDetails = TRUE,
             is_output_markerList_in_groupTest = TRUE,
             is_single_in_groupTest = FALSE,
             LOCO = FALSE,
             is_fastTest = FALSE,
             SAIGEOutputFile = opts$out
             )
             
print(f'DONE: {opts$out}')
