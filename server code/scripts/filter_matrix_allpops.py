import numpy as np
import pandas as pd

# Define a function to process each file
def process_file(file_prefix):
    # Load the frequencies file
    freqs = np.genfromtxt(f'{file_prefix}.frq', skip_header=1, usecols=(1,4))

    # Load the LD file
    ld = np.genfromtxt(f'{file_prefix}.ld', skip_header=False, filling_values=np.nan)

    # Load in the SNP positions
    bim = pd.read_csv(f'{file_prefix}.bim', sep='\t', names=['chrom','rsID','pos','ref','alt'])

    # Add frequencies
    bim['freq'] = freqs[:,1]

    # Filter to where frequencies are greater than zero
    nonzero_indices = np.nonzero(freqs[:,1])[0]
    bim_filtered = bim.iloc[nonzero_indices,:]
    
    # Write out filtered bim file
    bim_filtered.to_csv(f'{file_prefix}_filtered.bim', sep='\t', header=True, index=False)

    # Filter the LD matrix to keep only the SNPs with MAF > 0
    filtered_ld = ld[nonzero_indices][:,nonzero_indices]

    # Save the filtered LD matrix to a file
    np.savetxt(f'{file_prefix}_filtered.ld', filtered_ld, fmt='%g')


# Process each file
for file_prefix in ['YRI', 'GBR', 'EUR', 'AFR']:
    process_file(file_prefix)
