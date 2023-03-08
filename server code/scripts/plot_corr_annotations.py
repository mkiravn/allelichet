import pandas as pd

# Read in the combined annotated file
annotated = pd.read_csv('combined_annotated.txt', sep='\t')

# Read in the bim file for each population and concatenate them into one DataFrame
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import cm
import matplotlib as mpl
import matplotlib.colors as colors

populations = ["YRI", "EUR", "GBR", "AFR"]

for pop in populations:
    # Load data
    ld_matrix = np.loadtxt(f"{pop}_filtered.ld")
    np.fill_diagonal(ld_matrix, 0)
    #ld_matrix[np.triu_indices(ld_matrix.shape[0], k=1)] = np.nan
    bim = pd.read_csv(f'{pop}_filtered.bim', sep='\t')
    ld_matrix = pd.DataFrame(ld_matrix, index=bim['pos'], columns=bim['pos'])

    # Merge the annotated file with the bim file based on the position column
    df = pd.merge(annotated, bim, left_on='Start', right_on='pos')
    # create a new column to indicate the region type
    df['region'] = (df['Func.refGene'] != df['Func.refGene'].shift()).cumsum()

    # group by the region column to get separate groups for each contiguous region
    grouped = df.groupby('region')
    r_mat = np.zeros((grouped.ngroups, grouped.ngroups))
    axis_labels=[]
    # Loop over each group
    for name, group in grouped:
        axis_labels.append(group['Func.refGene'].mode()[0])
        for name_y, group_y in grouped:
            if name>=name_y:
                pos = group['pos'].values
                pos_y = group_y['pos'].values
                pairwise_ld = ld_matrix.loc[pos, pos_y]
                try:
                    r_mat[name-1,name_y-1]=np.nanmean(pairwise_ld)
                except:
                    print(f"Error with regions {name} and {name_y}; passing...")
                    pass
                
    r_mat[r_mat==0]=np.nan
    cmap = mpl.cm.get_cmap('seismic').copy()
    cmap.set_under('w')
    cmap.set_over('w')
    # Show the plot
    plt.figure()
    plt.imshow(r_mat, cmap=cmap, vmin=-1, vmax=1)
    labels=axis_labels
    plt.xticks(range(len(labels)), labels, fontsize=6, rotation=90)
    plt.yticks(range(len(labels)), labels, fontsize=6)
    cbar = plt.colorbar()
    cbar.set_label('r') 
    plt.tight_layout()
    plt.title(f'{pop}')
    plt.savefig(pop+'_annotationheatmap.pdf')


    
