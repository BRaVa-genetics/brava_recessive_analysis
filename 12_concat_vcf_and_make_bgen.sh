#!bin/bash

# Script that concatenates together all chromosome files to one, per annotation and encoding mode
# and also creates a bgen file for each set
# input: generated VCFs by the previous step
# output: one concatenated VCF, and BGEN + .sample file
# dependencies: bcftools, plink2, bgenix

# dev Sept 2023; refactored 17/7/2025

module load common-apps/bcftools/1.18

if [[ -z $1 ]]; then
    echo "Error: no file-tag was given! Maybe use PP90af05 ?"
    exit 1
else
    TAG=$1
    echo "Will processes files with $TAG."
fi
work_dir='/FIXTHIS/recessive_encoding/'

for consq in pLoF pLoF_damaging_missense nonsynonymous synonymous; do
    for mode in add rec; do
        vcf_prefix=$work_dir/vcf_ready/chrALL.${mode}_${TAG}.$consq
        bgen_prefix=$work_dir/bgen_genotypes/chrALL.${mode}_${TAG}.$consq

        # Step 1: concatenate several chromosome-based files to one
        echo $work_dir/vcf_ready/chr{1..22}.${mode}_${TAG}.$consq.vcf.gz | tr ' ' '\n' > $work_dir/files_to_concat.txt 
        bcftools concat -f files_to_concat.txt -Oz -o $vcf_prefix.vcf.gz

        # Step 2: Generate a BGEN file for this set of genotypes
        plink2 --export bgen-1.2 bits=8 \
            --vcf $vcf_prefix.vcf.gz dosage=DS \
            --out $bgen_prefix
        bgenix -index -g $bgen_prefix.bgen

        # Step 3: fix the samples file
        head -2 $bgen_prefix.sample > $bgen_prefix.tmp
        awk '(NR>2){print 1,$2,$3,$4}' $bgen_prefix.sample >> $bgen_prefix.tmp
        mv $bgen_prefix.tmp $bgen_prefix.sample
        rm $work_dir/files_to_concat.txt # redundant, but clean
    done
done

tar -czvf $work_dir/$TAG.rec_and_add.chrALL.tar.gz vcf_ready/chrALL.add_$TAG.* vcf_ready/chrALL.rec_$TAG.*
