#!/bin/bash

# Create a temporary directory for ANNOVAR output files
mkdir annovar_output

# Loop through all VCF files in the current directory
for vcf in *.vcf.gz; do
  bcftools norm -m-both -o step1.vcf $vcf
  # Run ANNOVAR on the VCF file and output the results to the temporary directory
  /project2/jnovembre/mianniravn/software/annovar/table_annovar.pl step1.vcf /project2/jnovembre/mianniravn/software/annovar/humandb/ -buildver hg19 -out annovar_output/${vcf%.vcf} -remove -protocol refGene -operation g -nastring . -vcfinput
done




