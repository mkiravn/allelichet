import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.ticker as mticker
label_format = '{:,.1f}'

populations = ["GBR", "YRI", "AFR", "EUR"]

for pop in populations:
    bim = pd.read_csv(f"{pop}_filtered.bim", sep="\t")
    ld = np.loadtxt(f"{pop}_filtered.ld")
    # set upper triangle to NaN
    ld[np.triu_indices(len(bim))] = np.nan

    n_pairs = 10_000
    random_pairs = np.random.choice(range(ld.shape[0]), size=(n_pairs, 2), replace=True)
    random_pairs = np.sort(random_pairs, axis=1)[: ,::-1] # sort so we get the right part of the matrix

    r_joint_freq=np.zeros(shape=(10,10,2)) # in third dimension, first index for positive, second for negative
    n_joint_freq=np.zeros(shape=(10,10,2))
    for snp1, snp2 in random_pairs:
        r = ld[snp1, snp2]
        f1 = np.around(bim.iloc[snp1, 5],1) * 10 # extract frequency
        f1=f1.astype(int)
        f2 = np.around(bim.iloc[snp2, 5],1) * 10 # extract frequency
        f2 = f2.astype(int)
        freqs=[f1,f2]
        freqs=np.sort(freqs)[::-1]
        # fill in the matrices - add r if positive or negative to right dimension of array
        # and add one to count
        if r>0: 
            r_joint_freq[freqs[0],freqs[1],0]+= r
            n_joint_freq[freqs[0], freqs[1], 0] += 1
        if r<0:
            r_joint_freq[freqs[0],freqs[1],1]+= r
            n_joint_freq[freqs[0], freqs[1], 1] += 1

    n_joint_freq[n_joint_freq==0]=np.nan
    r_joint_freq[n_joint_freq==0]=np.nan
    r_joint_freq = r_joint_freq/n_joint_freq # normalise by number of SNPs
    r_joint_freq[n_joint_freq==0]=np.nan
    # plotting 
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))
    for i, label in enumerate(["positive","negative"]):
        ax = axs[i]
        mat = r_joint_freq[:, :, i]
        nmat= n_joint_freq[:, :, i]
        im = ax.imshow(mat, cmap='coolwarm', vmin=-1, vmax=1)
        ax.set_title(f"{label}ly correlated")
        labels = np.arange(0,1,0.1)
        ax.xaxis.set_ticks(np.arange(0,10))
        ax.yaxis.set_ticks(np.arange(0,10))
        ax.set_yticklabels([label_format.format(x) for x in labels])
        ax.set_xticklabels([label_format.format(x) for x in labels])
        # Add labels to the heatmap
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if not np.isnan(nmat[i, j]):
                    text = ax.text(j, i, nmat[i, j].astype(int),
                                   ha="center", va="center", color="black", fontsize=5)
        for edge, spine in ax.spines.items():
            spine.set_visible(False)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.25, 0.02, 0.5]) # add colourbar
    cbar=fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('r')
    fig.suptitle(f"{pop} r plot (pos: {bim['pos'].min()}-{bim['pos'].max()})")
    plt.savefig(f"{pop}_joint_freq_r.pdf")