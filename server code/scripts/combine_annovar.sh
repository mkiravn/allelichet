#!/bin/bash

# Define the output file
output_file="combined_annotated.txt"

# Loop over the populations and append their data to the output file
for pop in AFR YRI EUR GBR; do
  # Get the current input file
  input_file="annovar_output/${pop}.vcf.gz.hg19_multianno.txt"
  
  # Get the header and first 10 fields of the input file
  header=$(head -n 1 "$input_file" | cut -f 1-10)
  
  # Append the header to the output file if it hasn't been added yet
  if [[ $pop == "AFR" ]]; then
    echo "$header" > "$output_file"
  fi
  
  # Append the data to the output file, removing duplicates
  tail -n +2 "$input_file" | cut -f 1-10 | sort -u >> "$output_file"
done
