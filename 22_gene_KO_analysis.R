library(dplyr)
library(data.table)
library(tidyr)
setwd('/Users/gk18/Documents/WORK/project_BRaVa_Recessive/gene_KOs')

gene_names <- read.table("map.hgnc.gene_id_transcript.txt", header = TRUE, sep = "\t", na.strings = c(" ", "NA"), fill = TRUE)
# gene_names <- gene_names %>% select( gene_id = 'Ensembl.gene.ID', Symbol='Approved.symbol')
gene_names <- gene_names %>%
  mutate(transcript_id = sub("\\..*", "", MANE_transcript_ID))

map_geneid_to_symbol <- function(df, gene_names, ID_tag) {
  df <- df %>% 
    left_join( gene_names, by=ID_tag ) %>%
    select( everything(), -MANE_transcript_ID  ) #, Symbol='hgnc_symbol', -transcript_id )
  return(df)
}
map_symbol_geneid <- function(df, gene_names, new_name, old_name) {
  df <- df %>% 
    rename( !!new_name := .data[[old_name]] ) %>%
    left_join( gene_names, by=new_name ) %>%
    select( everything(), -transcript_id, -MANE_transcript_ID )
  return(df)
}

################################################################
# Code to generate the input table from cohort-specific tables #

file_tags <- c('BBJ','BioMe','GNH','GEL','UKBB_all', 'AOU.pLoF_only')
data_brava <- lapply( file_tags, 
                      function(x){ read.table(paste(c('gene_KOs.PP90af05_s50',x,'txt'),collapse='.'), 
                                              sep='\t', header=TRUE)}
)
names(data_brava) <- c('BBJ','BioMe','GNH','GEL','UKB','AOU')

# homogenize all the tables - remember that all consequence cols should be named "consq"
data_brava[['BBJ']]$ancestry = 'EAS'
data_brava[['BBJ']]$biobank = 'BBJ'
data_brava[['GNH']]$ancestry = 'SAS'
data_brava[['GNH']]$biobank = 'GNH'
data_brava[['BioMe']]$biobank = 'BioMe'
data_brava[['UKB']]$biobank = 'UKB'
data_brava[['UKB']]$hgnc_symbol <- NULL
data_brava[['GEL']]$biobank = 'GEL'
data_brava[['GEL']]$ancestry <- toupper(data_brava[['GEL']]$ancestry)
data_brava[['AOU']]$ancestry <- toupper(data_brava[['AOU']]$ancestry)
data_brava[['AOU']]$biobank = 'AOU'

# drop CHR column from UKB:
# data_brava[['UKBB_450k']] <- data_brava[['UKBB_450k']][, !colnames(data_brava[['UKBB_450k']]) %in% "chr"]
data_brava <- lapply(data_brava, function(df) {
  if ("chr" %in% colnames(df)) df <- df[, !colnames(df) %in% "chr"]
  if ("hgnc_symbol" %in% colnames(df)) df <- df[, !colnames(df) %in% "hgnc_symbol"]
  
  df <- df %>%
    mutate(consq = recode(consq,
                          "damaging_missense" = "damaging_missense_or_protein_altering",
                          "pLoF_damaging_missense" = "pLoF_damaging"
    ))
  
  return(df)
})

for (x in names(data_brava)){
  # print(colnames(data_brava[[x]]))
  # print(unique(data_brava[[x]]$ancestry))
  print(unique(data_brava[[x]]$consq))
}

# combine all tables in one
dt_brava_all <- rbindlist(data_brava, use.names = TRUE) #fill = TRUE, idcol = "COHORT"
dt_brava_all <- dt_brava_all %>%
  select(everything(), gene_id='gene', csqs='consq' )

dt_brava_all %>%
  # filter( biobank=="UKB") %>%
  mutate(across(c('hom','chet'), ~ replace_na(., 0))) %>%
  mutate(kos = hom + chet) %>%
  filter(csqs == "pLoF") %>%
  group_by(gene_id, ancestry) %>%
  summarise(kos = sum(kos), .groups = "drop") %>%
  pivot_wider(names_from = ancestry, values_from = kos, values_fill = 0) %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  summarise(across(everything(), ~ sum(.x > 0)))

# Main work ---------------------------------------------------------------

# select knocked-out genes based on pLoF, replase NAs with zeros, and group by ancestry
brava_plofs = dt_brava_all %>%
  mutate(across(c('hom','chet'), ~ replace_na(., 0))) %>% 
  mutate(kos = hom + chet) %>%
  filter(csqs == "pLoF") %>%
  group_by(gene_id, ancestry) %>%
  summarise(kos = sum(kos), .groups = "drop") %>%
  pivot_wider(names_from = ancestry, values_from = kos, values_fill = 0) %>%
  mutate(total = rowSums(across(where(is.numeric))))

# get a summary per ancestry
brava_plofs %>% summarise(across(everything(), ~ sum(.x > 0)))

# bring in Symbol names
brava_plofs <- map_geneid_to_symbol(brava_plofs, gene_names, 'gene_id')
# brava_plofs <- map_symbol_geneid(brava_plofs, gene_names,'Symbol', 'gene_id')
fwrite(brava_plofs, 'BRaVa.gene_KOs.pLoF.txt', sep="\t", na="NA", quote=FALSE)


#################################
# comparison with other studies #
data_oddsson = read.table('Oddson_etal_geneKOs.csv', header=TRUE, sep=',')
data_sunetal = read.table('Sun_etal_RGCME_geneKOs.csv', header=TRUE, sep=',')
# WARNING: Sun et al did not work with MANE Selece trascripts, so better to stick to gene Symbols
data_oddsson <- data_oddsson %>% select( everything(), gene_id = gene)
data_sunetal <- data_sunetal %>% select( everything(), gene_id = GeneName)
brava_plofs <- brava_plofs %>% select( everything(), ensmbl=gene_id, gene_id = Symbol)

# data_oddsson <- map_symbol_geneid(data_oddsson,gene_names,'Symbol','gene')
# data_sunetal <- map_symbol_geneid(data_sunetal,gene_names,'Symbol','GeneName')
# data_sunetal <- map_symbol_geneid(data_sunetal,gene_names,'transcript_id','TranscriptId')

nKOs_sunetal = dim(data_sunetal %>% filter(sum_KO>0))[1]
nKOs_brava = dim(brava_plofs %>% filter(total>0))[1]
cat(paste(nKOs_brava, "vs", nKOs_sunetal))

KOs_brava = brava_plofs %>% filter(total>0) %>% pull(gene_id)
# KOs_brava = brava_plofs %>% filter(total>0) %>% pull(Symbol)
KOs_rgvme = data_sunetal$gene_id
n_missed = length(setdiff(KOs_rgvme, KOs_brava))
n_new    = length(setdiff(KOs_brava, KOs_rgvme))
n_shared = length(intersect(KOs_brava,KOs_rgvme))
cat(paste(length(KOs_brava), n_shared, n_new, n_missed))

# make a list of all genes observed across studies
genes_all <- brava_plofs %>% pull(gene_id)
for (x in c( "Sulem", "Narasimhan","Saleheen", "gnomAD", "Oddsson" )){
  genes_new <- data_oddsson %>% filter( data_oddsson[[x]] > 0 ) %>% pull(gene_id)
  cat(x, length(genes_new), length(setdiff(genes_new,genes_all)), '\n')
  genes_all <- union(genes_new,genes_all)
}

genes_all <- union(genes_all,data_sunetal$gene_id)
cat('Total number of genes across all studies:', length(genes_all))

missing_genes <- setdiff(genes_all, brava_plofs$gene_id)
# missing_genes <- setdiff(genes_all, brava_plofs$Symbol)
if (length(missing_genes) > 0) {
  missing_df <- data.frame(
    gene_id = missing_genes,
    EUR = 0, EAS = 0, SAS = 0, AFR = 0, AMR = 0, total = 0  )
  # Bind to original table
  brava_plofs <- bind_rows(brava_plofs, missing_df)
}

# now add info from other studies
for (x in c( "Sulem", "Narasimhan","Saleheen", "gnomAD", "Oddsson" )){
  genes_new <- data_oddsson %>% filter( data_oddsson[[x]] > 0 ) %>% pull(gene_id)
  brava_plofs <- brava_plofs %>%
    mutate( !!x:= gene_id %in% genes_new )
}

brava_plofs <- brava_plofs %>%
  mutate( BRaVA = total>0) %>%
  # mutate(Regeneron = Symbol %in% data_sunetal$GeneName)
  mutate(Regeneron = gene_id %in% data_sunetal$gene_id)

brava_plofs %>%
  summarise(across( everything(), ~ sum(.x > 0)))

# almost done, now prettify
brava_plofs <- brava_plofs %>%
  mutate(across(c(,"BRaVA","Sulem","Narasimhan","Saleheen","gnomAD","Oddsson","Regeneron"), ~ as.integer(.))) %>%
  select( everything(), Gene='gene_id', BRaVA_total='total', ensmbl_id='ensmbl') %>%
  mutate( BRaVa_novel = as.integer( (BRaVA==1) & (Sulem+Narasimhan+Saleheen+gnomAD+Oddsson+Regeneron==0) ) ) %>%
  filter( BRaVA+Sulem+Narasimhan+Saleheen+gnomAD+Oddsson+Regeneron>0 ) %>%
  arrange(Gene)
  # mutate( shared = as.integer( (BRaVA==1) & (Sulem+Narasimhan+Saleheen+gnomAD+Oddsson+Regeneron==1) ) ) 

brava_plofs <- brava_plofs[c("Gene","BRaVA_total","AFR","AMR","EUR","EAS","SAS","BRaVA",
                             "Sulem","Narasimhan","Saleheen","gnomAD","Oddsson","Regeneron","BRaVa_novel","ensmbl_id","transcript_id")]
                             
fwrite(brava_plofs, 'BRaVa.gene_KOs.pLoF.with_other_studies.txt', sep="\t", na="NA", quote=FALSE)

# Get summaries:

brava_plofs %>%
  filter( BRaVa_novel>0 ) %>%
  summarise(across( everything(), ~ sum(.x > 0)))

brava_plofs %>%
  filter(BRaVA_total>0) %>%
  mutate( Other_studies = as.integer(Sulem+Narasimhan+Saleheen+gnomAD+Oddsson+Regeneron>0 )) %>%
  summarise(
    total_genes = n(),
    found_elsewhere = sum(Other_studies),
    new_genes = total_genes - found_elsewhere,
    fraction_found = found_elsewhere / total_genes,
    fraction_new = new_genes / total_genes
  )

brava_plofs %>%
  mutate( Other_studies = Sulem+Narasimhan+Saleheen+gnomAD+Oddsson+Regeneron ) %>%
  summarise(
    total_genes = n(),
    found_elsewhere = sum(Other_studies>0),
    new_genes = sum(BRaVA_total>0 & Other_studies==0),
    fraction_found = found_elsewhere / total_genes,
    fraction_new = new_genes / total_genes
  )

brava_plofs %>%
  filter(BRaVA_total==0) %>%
  summarise(across(everything(), ~ sum(.x > 0)))

brava_plofs %>%
  mutate( Other_studies = as.integer(Sulem+Narasimhan+Saleheen+gnomAD+Oddsson+Regeneron>0 )) %>%
  filter(Other_studies>0) %>%
  summarise(across(everything(), ~ sum(.x > 0)))

# other stuff -------------------------------------------------------------

# count the increase in testable genes due to CH
annot="pLoF_damaging_missense"
annot="nonsynonymous"
annot='pLoF'

for (annot in unique(dt_brava_all$csqs)){
    
  dt_testable_genes <- dt_brava_all %>%
    mutate( cohort = paste0(biobank,":",ancestry) ) %>%
    mutate( kos = hom) %>%
    # mutate( kos = hom+chet) %>%
    # filter(csqs == annot & kos>4 & biobank != "AOU") %>%
    filter(csqs == annot & kos>0 & biobank != "AOU") %>%
    # filter( csqs == annot & (hom>0 | chet>0) ) %>%
    # filter( csqs == annot & hom>0 ) %>%
    group_by(gene_id, cohort)
  
  cat(paste("Testable genes for",annot,": ",length(unique(dt_testable_genes$gene_id)),'\n'))
}

# genes per cohort; use `kos = hom` and `kos = chet` 
dt_brava_all %>%
  filter( csqs=="pLoF") %>%
  mutate( cohort = paste0(biobank,":",ancestry), kos = chet ) %>%
  # mutate( ) %>%
  pivot_wider(names_from = cohort, values_from = kos, values_fill = 0) %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  summarise(across(everything(), ~ sum(.x > 0))) %>%
  pull(`GNH:SAS`, `UKB:EUR`)

# get the corresponding number of singletons
n_singletons = dim(brava_plofs %>% filter( BRaVA_total==1))[1]
cat("Singletons: ", n_singletons, 100*n_singletons/nKOs_brava)

n_more = dim(brava_plofs %>% filter( BRaVA_total>1))[1]
cat("Two or more: ", n_more, 100*n_more/nKOs_brava)

# get number of genes in one or more biobanks
dt_brava_all %>%
  mutate(across(c('hom','chet'), ~ replace_na(., 0))) %>% 
  mutate(kos = hom + chet) %>%
  filter(csqs == "pLoF" & kos>1) %>%
  group_by(gene_id, biobank) %>%
  summarise(kos = sum(kos), .groups = "drop") %>%
  pivot_wider(names_from = biobank, values_from = kos, values_fill = 0) %>%
  mutate(nonzero_cohorts = rowSums(select(., -gene_id) > 0)) %>%
  filter( nonzero_cohorts==1 )
  # mutate(total = rowSums(across(where(is.numeric))))

nrow( dt_brava_all %>%
  mutate(across(c('hom','chet'), ~ replace_na(., 0))) %>% 
  filter(csqs == "pLoF" & chet>0 & biobank=="UKB" & ancestry=="EUR") )

dt <- dt_brava_all %>%
  mutate( cohort = paste0(biobank,":",ancestry) ) %>%
  filter(csqs == 'pLoF' & chet>0 & biobank != "aou") %>%
  group_by(gene_id, cohort) %>%
  pivot_wider(names_from = cohort, values_from = chet, values_fill = 0)

dt = dt[colnames(dt)[9:21]]
dt %>% summarise(across(everything() , ~ sum(.x > 0)))
