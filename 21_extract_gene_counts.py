#!/usr/bin/env python3
# a script to get a quick tally of the number of KOs per gene, per genotype type
# given the output of `called_chets`

import pandas as pd
import numpy as np

work_dir='/FIX_THIS/recessive_encoding'
tag='PP90af05_s50'

carrier_summaries = {}
for consq in ['pLoF','damaging_missense_or_protein_altering','pLoF_damaging', 'nonsynonymous', 'synonymous']:
    carrier_summaries[consq] = []
    for c in range(1,23):
        df = pd.read_csv(f'{work_dir}/sandbox/chr{c}.gen_{tag}.{consq}.txt',
                         sep='\t', header=None, names=['id','chr','gene','gt','ac','allele'])
        carrier_summaries[consq].append(
            df.pivot_table(index='gene', columns='gt', aggfunc='size').replace({np.nan:0})
        )

    print("Read in all chromosomes for", consq, "; concatenating...")
    carrier_summaries[consq] = pd.concat(carrier_summaries[consq], axis=0, join='inner')

pd.concat( carrier_summaries ).astype(int).to_csv(f'gene_KOs.{tag}.txt', sep='\t', na_rep=0, float_format='%.0f')

# end-of-script
