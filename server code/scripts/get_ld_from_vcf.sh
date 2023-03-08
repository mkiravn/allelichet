#!/bin/bash

# set paths

bcftools_path="bcftools"
plink_path="plink"
vep_path="vep"
vep_cache_dir="path/to/vep_cache_dir"

pops="/project2/jnovembre/mianniravn/1kg_phase3/haps/integrated_call_samples_v3.20130502.ALL.panel"

# get chromosome and positions from command line arguments
locus_chrom="$1"
locus_start="$2"
locus_end="$3"
vcfpath="/project2/jnovembre/data/external_public/1kg_phase3/haps/ALL.chr$locus_chrom.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
# create output directory
output_dir="chr${locus_chrom}_${locus_start}-${locus_end}"
mkdir -p "$output_dir"

# loop through populations
for pop in GBR YRI AFR EUR
do
    # extract sample IDs for population
    grep $pop $pops | awk '{print $1}' > "$output_dir/$pop.samples.txt"

    # filter vcf by population and locus using bcftools
    bcftools view -Oz -f "PASS" -S "$output_dir/$pop.samples.txt" --force-samples -r "$locus_chrom:$locus_start-$locus_end" $vcfpath > "$output_dir/$pop.vcf.gz"

    bcftools index "$output_dir/$pop.vcf.gz"  # index the output VCF file

    # calculate LD using plink
    plink --vcf "$output_dir/$pop.vcf.gz" --r square0 --out "$output_dir/$pop"

    # make BED file
    plink --vcf "$output_dir/$pop.vcf.gz" --make-bed --out "$output_dir/$pop"

    # get frequencies
    plink --vcf "$output_dir/$pop.vcf.gz" --freq --out "$output_dir/$pop"
done
