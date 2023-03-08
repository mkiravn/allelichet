import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors

populations = ["GBR", "YRI", "AFR", "EUR"]

for pop in populations:
    bim = pd.read_csv(f"{pop}_filtered.bim", sep="\t")
    ld = np.loadtxt(f"{pop}_filtered.ld")
    # set upper triangle to NaN
    ld[np.triu_indices(len(bim))] = np.nan

    # plotting the LD structure
    plt.imshow(ld, cmap='seismic', vmin=-1, vmax=1)
    cbar = plt.colorbar()
    cbar.set_label('r')
    x_tick_positions = np.arange(0, len(bim), step=50)
    y_tick_positions = np.arange(0, len(bim), step=50)
    plt.xticks(x_tick_positions, bim['pos'].iloc[x_tick_positions].astype(int), rotation=90, fontsize=6)
    plt.yticks(y_tick_positions, bim['pos'].iloc[y_tick_positions].astype(int), fontsize=6)
    plt.title(f"{pop} LD (pos: {bim['pos'].min()}-{bim['pos'].max()})")
    plt.xlabel("chromosomal position")
    plt.ylabel("chromosomal position")
    plt.savefig(f"{pop}_ld.pdf")
    plt.close()

    plt.figure()
    # plotting mean frequency versus r for a random subset of SNP pairs
    n_pairs = 10_000
    random_pairs = np.random.choice(range(ld.shape[0]), size=(n_pairs, 2), replace=True)
    random_pairs = np.sort(random_pairs, axis=1)[:,::-1] # sort so we get the right part of the matrix
    bounds = [0, 1000, 10000, 100000, 1000000]
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256, extend='both')

    for snp1, snp2 in random_pairs:
        r = ld[snp1, snp2]
        mnfreq = (bim.iloc[snp1, 5] + bim.iloc[snp2, 5]) / 2  # mean frequency
        dist = bim.iloc[snp1]["pos"] - bim.iloc[snp2]["pos"]
        dist = np.abs(dist)
        plt.scatter(mnfreq, r, c=dist, s=2, cmap=mpl.cm.get_cmap("magma"),norm=norm)

    plt.ylim((-1,1))
    cbar = plt.colorbar()
    cbar.set_label('Distance (bp)')
    plt.title(f"{pop} r plot (pos: {bim['pos'].min()}-{bim['pos'].max()})")
    plt.xlabel("mean allele frequency")
    plt.ylabel("r")
    plt.savefig(f"{pop}_r_plot.pdf")
    plt.close()
