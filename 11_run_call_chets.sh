#!/bin/bash

# this is a wrapper that combines several tools towards a pipeline for CompHet calling.
# see here https://github.com/frhl/call_chets (needs to be installed).
# input: phased genotypes, annotations in BRaVa/SAIGE format
# output: VCFs based on different consequences, for either additive or recessive effects

# NOTES: 
#   > update all /FIXTHIS/ flags to your local paths
#   > main argument should be chromosome index
#   > remember to choose the correct tag between damaging_missense and damaging_missense_or_protein_altering

# dev by 11/09/2023; refactored July 17/07/25

if [[ -z $1 ]]; then
    CHR=$LSB_JOBINDEX
else
    CHR=$1
fi
echo "Will create VCFs with biallelic encodings for chr-$CHR."

module load common-apps/bcftools/1.18
cpp_dir="/FIXTHIS/call_chets/"
work_dir="/FIXTHIS/phasing"
annot="/FIXTHIS/annotation.chr$CHR.txt"

out_dir="$work_dir/recessive_encoding/vcf_ready"
tmp_dir="$work_dir/recessive_encoding/sandbox"
BCF="$work_dir/phased_genotypes_rare/GNH_39k.chr$CHR.phased_rare.no100trios.bcf"
genotypes="$tmp_dir/chr$CHR.nonref.PP90af05.txt.gz"

# make a list of all samples from the BCF, if not already present
if [ ! -f $tmp_dir/samples.txt ]; then
    bcftools query $BCF --list samples > $tmp_dir/samples.txt
fi

# Step 1: Generate genotypes if not already available (quite slow due to large input size)
if [ ! -f $genotypes ]; then 
    echo "Generating a file with all genotypes and filtering for PP>90%"
    bcftools view $BCF --max-af 0.05 -Ou | bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT %PP\n]' | \
    awk '($4=="." || $4>0.90)' | gzip > $genotypes
fi

# Step 2: Prepare list of variants from annotation file
for consq in pLoF damaging_missense synonymous; do
    python prepare_genemap.py -a $annot -c $consq -o $tmp_dir/gene_map.chr$CHR.$consq.txt
    python prepare_genemap.py -a $annot -c $consq -o $tmp_dir/gene_map.chr$CHR.$consq.txt
done
# repeat for combinations of consequences
python prepare_genemap -a $annot -c pLoF damaging_missense -o $tmp_dir/gene_map.chr$CHR.pLoF_damaging_missense.txt
python prepare_genemap -a $annot -c pLoF damaging_missense other_missense -o $tmp_dir/gene_map.chr$CHR.nonsynonymous.txt

# Step 3: call comphets and encode for additive or recessive effects and generate VCF files
for consq in pLoF pLoF_damaging_missense nonsynonymous synonymous; do
    $cpp_dir/call_chets -g $genotypes -m $tmp_dir/gene_map.chr$CHR.$consq.txt --show-variants > $tmp_dir/chr$CHR.gen_all.$consq.txt
    $cpp_dir/encode_vcf -i $tmp_dir/chr$CHR.gen_all.$consq.txt -s $tmp_dir/samples.txt -m additive  | bgzip > $out_dir/chr$CHR.add.$consq.vcf.gz
    $cpp_dir/encode_vcf -i $tmp_dir/chr$CHR.gen_all.$consq.txt -s $tmp_dir/samples.txt -m recessive | bgzip > $out_dir/chr$CHR.rec.$consq.vcf.gz

    # print some stats (optional)
    for bial in cis chet hom; do
        grep $bial $tmp_dir/chr$CHR.gen_all.$consq.txt > $tmp_dir/chr$CHR.$bial.$consq.txt
        tmp=$(wc -l $tmp_dir/chr$CHR.$bial.$consq.txt | cut -d' ' -f1)
        echo "$bial-$consq events found: $tmp"
    done
done

echo "Script completed for chromosome $CHR."
# end-of-script
