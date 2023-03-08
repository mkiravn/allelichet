import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

populations = ["YRI", "EUR", "GBR", "AFR"]

for pop in populations:
    # Load data
    bim = pd.read_csv(f"{pop}_filtered.bim", sep="\t")
    ld = np.loadtxt(f"{pop}_filtered.ld")
    np.fill_diagonal(ld, 0)

    # Calculate nonzero matrix
    num_chunks = 10
    chunk_size = len(ld) // num_chunks
    print(f"{chunk_size} SNPs per chunk for population {pop}.")
    pos = bim['pos'].values
    pos_labels = []
    nonzero_matrix = np.zeros((num_chunks, num_chunks, 3))
    for i in range(num_chunks):
        start_idx = (i) * chunk_size
        end_idx = start_idx + chunk_size
        pos_end = pos[end_idx - 1]
        pos_labels.append(pos_end)
        for j in range(i + 1):
            start_idx_col = (j) * chunk_size
            end_idx_col = start_idx_col + chunk_size
            chunk = ld[start_idx:end_idx, start_idx_col:end_idx_col]
            num_nonzero = len(chunk[np.nonzero(chunk)])
            try:
                nonzero_matrix[i, j, 0] = len(chunk[chunk < 0]) / num_nonzero
                nonzero_matrix[i, j, 1] = (len(chunk[chunk > 0]) - len(chunk[chunk == 1]))/ num_nonzero
                nonzero_matrix[i, j, 2] = len(chunk[chunk == 1]) / num_nonzero
            except:
                print("Problem: all zero entries or error.")
                for i in range(2):
                    nonzero_matrix[i, j, i] = np.nan
    nonzero_matrix[nonzero_matrix == 0] = np.nan

    # Create plot
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(12, 5))
    for i, label in enumerate(["negative","positive", "perfect"]):
        ax = axs[i]
        mat = nonzero_matrix[:, :, i]
        im = ax.imshow(mat, cmap='coolwarm', vmin=-1, vmax=1)
        ax.set_title(f"{label}ly correlated")
        axs[i].set_xticks(np.arange(num_chunks) - 0.5, minor=True)
        axs[i].set_yticks(np.arange(num_chunks) - 0.5, minor=True)
        axs[i].set_xticks(np.arange(num_chunks))
        axs[i].set_yticks(np.arange(num_chunks))
        # Set the tick labels
        xlabels = [f"{int(bim['pos'][i * chunk_size]):,}" for i in range(num_chunks)]
        ylabels = [f"{int(bim['pos'][i * chunk_size]):,}" for i in range(num_chunks)]
        axs[i].set_xticklabels(xlabels, fontsize=8, rotation=90)
        axs[i].set_yticklabels(ylabels, fontsize=8)
        # Set the tick locations to the end of the cell
        axs[i].tick_params(axis="both", which="both", length=0, labelsize=6, pad=-2, direction="out")
        axs[i].tick_params(axis="both", which="minor", length=5, direction="out")
        # Shift the y-axis tick labels to the end of each chunk
        if i == 0:
            ax.tick_params(axis='y', direction='out', pad=5)
        else:
            ax.tick_params(axis='y', direction='out', pad=20)
        ax.set_yticklabels([])

        # Set the tick label font size
        ax.tick_params(axis='both', labelsize=8)

        # Add labels to the heatmap
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if not np.isnan(mat[i, j]):
                    text = ax.text(j, i, round(mat[i, j], 2),
                                   ha="center", va="center", color="black", fontsize=6)
        for edge, spine in ax.spines.items():
            spine.set_visible(False)

        # Set the overall title of the figure
        fig.suptitle(f"Correlation Probabilities for {pop}: {min(bim['pos'])}-{max(bim['pos'])} \n{chunk_size} SNPs per chunk.", fontsize=12)
        fig.savefig(f"{pop}_corr_probs.pdf", bbox_inches='tight')