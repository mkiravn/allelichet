import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import cm
import matplotlib.colors as colors
import pandas as pd
import itertools

# Read in the combined annotated file
annotated = pd.read_csv('combined_annotated.txt', sep='\t')

# Define a function to group SNPs by distance
def group_distance(distance):
    if distance < 100:
        return '<0.1kb'
    elif distance < 1000:
        return '0.1-1kb'
    else:
        return '>1kb'

def plot_heatmap(populations,annotation):
    # plots a heatmap where we have pairs of functional annotations, and 
    # the colours are the mean value of r for SNPs within those categories
    # also separates by distance
    for pop in populations:
        print(f"Making plot for {annotation} annotation in population {pop}.")
        result = []
        pairs = []
        rs_counts = {}
        # Load data
        ld_matrix = np.loadtxt(f"{pop}_filtered.ld")
        np.fill_diagonal(ld_matrix, 0)

        bim = pd.read_csv(f'{pop}_filtered.bim', sep='\t')
        ld_matrix = pd.DataFrame(ld_matrix, index=bim['pos'], columns=bim['pos'])

        # Merge the annotated file with the bim file based on the position column
        df = pd.merge(annotated, bim, left_on='Start', right_on='pos')
        # create a new column to indicate the region type

        df = df[['Chr', 'pos', annotation]].drop_duplicates()
        df = df[df[annotation] != '.']

        # Compute the mean r for each pair of SNPs with the same annotation and distance
        grouped = df.groupby([annotation])
        pairs = list(itertools.combinations_with_replacement(grouped.groups.keys(), 2))
        for pair in pairs: # extract pairs of groups from the annotations
            name=pair[0]
            name2=pair[1]
            group = grouped.get_group(pair[0])
            group2=grouped.get_group(pair[1])
            pos = group['pos']
            pos2 = group2['pos']
            rs = ld_matrix.loc[pos, pos2]
            # Count the number of non-null values in each cell of the rs matrix
            counts = rs.count()

            # Add the counts to the dictionary
            rs_counts[(name, name2)] = counts

            distance_matrix = pd.DataFrame(0, index=pos, columns=pos2)
            for i in pos:
                for j in pos2:
                    distance_matrix.loc[i, j] = group_distance(abs(i - j))

            distances = distance_matrix.loc[pos, pos2]

            distances_mask = (distances == '<0.1kb') | (distances == '0.1-1kb') | (distances == '>1kb')
            distances = distances[distances_mask]
            rs = rs[distances_mask]
            if sum(rs!=np.nan)<1:
                # handle empty array
                pass
            else:
                for d, r in zip(np.unique(distances), np.nanmean(rs,axis=0)):
                    result.append(
                        {'Annotation_1': name, 'Annotation_2': name2, 'distance': d, 'mean_r': r})

        # Convert the result to a data frame
        result_df = pd.DataFrame(result)
        pivoted_df = result_df.pivot(index=["Annotation_1", "Annotation_2"], columns=["distance"], values="mean_r")
        # specify the desired order of the columns
        column_order = ['<0.1kb', '0.1-1kb', '>1kb']
        # reindex the pivoted dataframe with the desired order of the columns
        pivoted_df = pivoted_df[column_order]
        # make figure
        plt.figure(figsize=(5, 7))
        sns.set(font_scale=0.7)
        heatmap=sns.heatmap(pivoted_df, annot=False, fmt='.3f',
                    cmap='seismic', vmin=-1, vmax=1, linewidths=1, square=True )

        # Adjust the margins to make space for x-axis labels
        plt.subplots_adjust(bottom=0.4)
        plt.tight_layout()
        plt.xlabel('Distance')
        # Set the color bar font size
        plt.ylabel('Functional Annotation')
        # Split the labels into two lines if necessary
        for label in plt.yaxis.get_ticklabels():
            label.set_fontsize(10) # Set the font size
            label.set_wrap(True)   # Allow wrapping
            label.set_linespacing(0.5) # Adjust line spacing
            label.set_text(label.get_text().wrap(10)) # Specify width and split the label text
  
        plt.title(f'{pop} - {annotation}')
        plt.savefig(f"{pop}_{annotation}.pdf")

# define the annotations and populations
annotations = ["Func.refGene", "ExonicFunc.refGene"]
populations = ["YRI","AFR","GBR","EUR"]
# loop over each population and annotation
for anno in annotations:
    plot_heatmap(populations,anno)


