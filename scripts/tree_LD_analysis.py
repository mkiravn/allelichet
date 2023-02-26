"""Functions to help in analysing single locus trees simulated via msprime"""

import numpy as np
import msprime as msp
import tskit
import random
import pandas as pd
import matplotlib.pyplot as plt
import itertools

from msp_sim import SimulationResults


def allele_frequencies(ts):
    """Obtain allele frequencies for all mutations in a tree sequence
    Arguments
    ---------
    ts : msprime.TreeSequence
        tree sequence object
    Returns
    -------
    afs : list
        allele frequencies for the 1:m mutations
    """
    # option to specify samples
    sample_sets = [ts.samples()]
    # get number of samples
    n = np.array([len(x) for x in sample_sets])

    # frequency calculator
    def f(x):
        return x / n

    afs = ts.sample_count_stat(sample_sets,
                               f,  # apply the frequency calculator
                               len(sample_sets),
                               windows='sites',
                               polarised=True,
                               mode='site',
                               strict=False,
                               span_normalise=False)
    return (afs)

def get_row(ts,rnd=True,return_index=False,mask=0,signed=True):
    # get positions and frequencies
    freqs = np.array(allele_frequencies(ts).T.tolist()[0])
    # make a mask of sites which pass frequency cutoff
    mask = freqs > mask
    # sort by frequency
    if rnd:
        muts = random.sample(list(np.where(mask == 1)[0]), 2)  # the indices 2 randomly chosen mutations
    else:
        muts = np.argsort(freqs)[::-1][:2]  # the indices of the top 2 mutations
    fs = freqs[muts]  # their frequencies
    unique_mutations=np.unique(list(ts.tables.mutations.site),return_index=True)[1]
    derived = ts.tables.mutations.derived_state[unique_mutations]  # getting the derived allele
    # make a haplotype matrix by extracting haplotypes and comparing to derived allele
    haps = np.array([np.array(list(map(ord, hap))) == derived for hap in list(ts.haplotypes())])
    # get a signed LD matrix
    cor_matrix = np.corrcoef(haps.T, rowvar=True)
    # square it to get r^2 matrix
    if signed != True:
        cor_matrix = cor_matrix ** 2
    # extract r2 of the two mutations
    corr = cor_matrix[muts[0], muts[1]]
    row = np.array([[corr, fs[0], fs[1]]])  # make a row of (corr,pa,pb)
    row
    if return_index:
        return(row,muts)
    else:
        return(row)

def sim_replicates(sample_size=100,
                   num_replicates=1,
                   pop_size=10_000,
                   seq_length=1000,
                   mut_rate=1e-7,
                   quiet=True,
                   mask=0.01):
    # intialise array
    arr=np.empty([0,3])
    # simulate
    ancestry_reps = msp.sim_ancestry(
        samples=sample_size,
        population_size=pop_size,
        sequence_length=seq_length,
        num_replicates=num_replicates)
    # extract info from tree sequences
    for ts in ancestry_reps: # get each ts
        mutated_ts = msp.sim_mutations(ts, rate=mut_rate) # mutate
        # get the frequencies and r2 and iteratively add to arr
        try:
            r=get_row(mutated_ts,mask=mask,signed=True)
            arr=np.vstack((arr,r))
            arr=pd.DataFrame(arr)
            arr.columns =['r2', 'f1', 'f2']
        # handling the failing due to recurrent mutations at a site
        except Exception as e:
            # in development- troubleshooting here
            if quiet==False:
                print(e,mutated_ts.tables.sites)
            pass
    print(f"Failed {100-(arr.shape[0]/num_replicates) * 100}% of attempts.")
    return(arr)

def get_nodes_above(ts,mut_id):
    """traverses the tree upwards from a mutation to retrieve all nodes avoce"""
    above={}
    u=list(ts.tables.mutations.node)[mut_id]
    for tree in ts.trees():
        path = []
        v = u
        while v != tskit.NULL:
            path.append(v)
            v = tree.parent(v)
        above[u]=path
    return(path)

def check_same_branch(ts,mut_a,mut_b):
    """checks whether two mutations are on the same branch"""
    same_branch=ts.tables.mutations.node[mut_b]==ts.tables.mutations.node[mut_a]
    return(same_branch)

def check_above(ts,mut_a,mut_b):
    """checks whether two mutations are on the same side of a tree"""
    is_above=ts.tables.mutations.node[mut_b] in get_nodes_above(ts,mut_a) or ts.tables.mutations.node[mut_a] in get_nodes_above(ts,mut_b)
    return(is_above)

def find_topology(ts,mut_a,mut_b):
    """checks the configuration of two mutations"""
    same_branch=False
    is_above=False
    opposite_branch=False
    other_side=False
    if check_same_branch(ts,mut_a,mut_b): # check whether they fall on the same branch
        same_branch=True
    elif check_above(ts,mut_a,mut_b): # check whether one is above the other
        is_above=True
    elif sum(allele_frequencies(ts)[[mut_a,mut_b],])<1: # above implies that the two mutations are on different sides
        other_side=True
    else: # the case where the two mutations are on opposite sides of the root on the furthest internal branch
        opposite_branch=True
    topology=[is_above,other_side,same_branch,opposite_branch]
    return(topology)

def find_topology_2(ts,mut_a,mut_b):
    """checks the configuration of two mutations"""
    same_branch=False
    is_above=False
    opposite_branch=False
    other_side=False
    if check_same_branch(ts,mut_a,mut_b): # check whether they fall on the same branch
        same_branch=True
    elif check_above(ts,mut_a,mut_b): # check whether one is above the other
        is_above=True
    else: # the case where the two mutations are on opposite sides of the root on the furthest internal branch
        other_side=True
    topology=[is_above,other_side,same_branch]
    return(topology)


def sim_replicates_topology(sim,
                            reps=1,
                            quiet=True,
                            print_errors=False,
                            mask=0.01,
                            signed=True):
    """
    Params:
        * sim - a Simulation() object of whatever demography you choose
        * reps - number of replicates to run
        * mask - an allele frequency cut-off
        * signed - whether to return signed ld matrix
    Output:
        * arr - an (reps x 3) array with [r^2 or r,f1,f2] where f1 and f2 are the freqs
            of two randomly chosen mutations and r^2 is the r^2 between them
        * tops - a (reps x 4) tuple of topologies, in order counts of
                (is_above,other_side,same_branch,opposite_branch)
        * fracs - (reps x len(pop_sizes)) array containing the fraction
                of the tree made up from internal branches.
    """
    # simulates replicate msprime simulations and extracts topology of mutations and r2
    arr = np.empty([0, 3])
    tops = np.empty([0, 4])
    fracs = np.empty([0, 1])
    stats_res = {"cor": np.empty([0, 3]),
                 "t": np.empty([0, 3]),
                 "c": np.empty([0, 3])}
    # simulate
    # do this outside?
    ress = sim._sim_replicates(samples=100, reps=reps)
    # extract trees
    for replicate_index, mutated_ts in enumerate(ress):
        if mutated_ts.num_sites == 0:
            print("No mutations. Passing....")
        else:
            try:
                r, muts = get_row(mutated_ts, return_index=True, mask=mask, signed=signed)
                arr = np.vstack((arr, r[0]))
                tpl = find_topology(mutated_ts, muts[0], muts[1])
                frac_internal = internal_branch_length(mutated_ts)
                fracs = np.vstack((fracs, frac_internal))
                tops = np.vstack((tops, tpl))
                stats = stats_from_trees(mutated_ts, pop_size=sim.Ne)
                for var in ["cor", "t", "c"]:
                    stats_res[var] = np.vstack([stats_res[var], stats[var]])
            except Exception as e:
                # catch simulations which fail (due to recurrent mutations at a site usually)
                if quiet == False:
                    print(e)
                pass
    arr = pd.DataFrame(arr)
    arr.columns = ['r2', 'f1', 'f2']
    if print_errors:
        print(f"Failed {100 - (arr.shape[0] / reps) * 100}% of attempts.")
    return ((arr, tops, fracs, stats_res))


def topology_pop_size(sim, num_replicates=1, pop_sizes=[10_000],
                      mask=0.01,
                      signed=True):
    """simulates replicates of a given simulation, for different (initial) population sizesParams:
        * sim - a Simulation() object of whatever demography you choose
        * num_replicates - number of replicates to run
        * pop_sizes - list of initial population sizes
    Output: a SimulationResults object with the following attributes
        * rs - a (reps x len(pop_sizes) x 3) array with [r^2,f1,f2] where f1 and f2 are the freqs
            of two randomly chosen mutations and r^2 is the r^2 between them
        * tops - a (reps x len(pop_sizes) x 4) tuple of topologies, in order counts of
                (is_above,other_side,same_branch,opposite_branch)
        * fracs - a (reps x len(pop_sizes) x 1) array of the fraction of internal branch length
        * sim and pop_sizes
    """
    # initialise arrays
    rs = np.zeros(shape=(num_replicates, len(pop_sizes), 3))
    tops = np.zeros(shape=(num_replicates, len(pop_sizes), 4))
    fracs = np.zeros(shape=(num_replicates, len(pop_sizes)))

    for i in range(len(pop_sizes)): # for each of the pop sizes
        ps = pop_sizes[i]
        sim.Ne = ps
        reps = sim_replicates_topology(sim, reps=num_replicates,mask=mask,signed=signed,print_errors=True)
        # get the values of r2
        r2 = reps[0]
        missing=num_replicates-r2.shape[0]
        filler=np.full([missing, 3], np.nan)
        r2 = np.concatenate((r2, filler), axis=0)
        rs[:, i, :] = r2
        # summarise the tree topologies
        filler=np.full([missing, 4], np.nan)
        tps = np.concatenate((reps[1], filler), axis=0)
        tops[:, i, :] = tps
        # fraction internal branch length
        filler=np.full([missing, 1], np.nan)
        frcs = np.concatenate((reps[2], filler), axis=0)
        fracs[:, i] = frcs.T

    results = SimulationResults()
    results.pop_sizes = pop_sizes
    results.sim = sim
    results.r2 = rs
    results.tops = tops
    results.fracs = fracs
    results.reps = num_replicates

    return (results)


def topology_growth_rate(sim,
                         num_replicates=1,
                         growth_rates=[1],
                         plot_dems=False,
                         mask=0.01,
                         signed=True):
    """simulates replicates of a given simulation, for different (initial) population sizesParams:
        * sim - a Simulation() object of whatever demography you choose
        * num_replicates - number of replicates to run
        * growth_rates - list of growth rates
    Output: a SimulationResults object with the following attributes
        * rs - a (reps x len(pop_sizes) x 3) array with [r^2,f1,f2] where f1 and f2 are the freqs
            of two randomly chosen mutations and r^2 is the r^2 between them
        * tops - a (reps x len(pop_sizes) x 4) tuple of topologies, in order counts of
                (is_above,other_side,same_branch,opposite_branch)
        * fracs - a (reps x len(pop_sizes) x 1) array of the fraction of internal branch length
        * also records the simulation object and the growth rates
    """
    # initialise arrays
    rs = np.zeros(shape=(num_replicates, len(growth_rates), 3))
    tops = np.zeros(shape=(num_replicates, len(growth_rates), 4))
    fracs = np.zeros(shape=(num_replicates, len(growth_rates)))
    stats= {}


    for i in range(len(growth_rates)): # for each of the pop sizes

        gr = growth_rates[i]
        sim.r = gr
        sim._make_demography(r=gr)
        if plot_dems:
            sim._demography_plot()
        reps = sim_replicates_topology(sim, reps=num_replicates,
                                       mask=mask,
                                       signed=signed)
        # get the values of r2
        r2 = reps[0]
        missing=num_replicates-r2.shape[0]
        filler=np.full([missing, 3], np.nan)
        r2 = np.concatenate((r2, filler), axis=0)
        rs[:, i, :] = r2
        # summarise the tree topologies
        filler=np.full([missing, 4], np.nan)
        tps = np.concatenate((reps[1], filler), axis=0)
        tops[:, i, :] = tps
        # fraction internal branch length
        filler=np.full([missing, 1], np.nan)
        frcs = np.concatenate((reps[2], filler), axis=0)
        fracs[:, i] = frcs.T
        stats[i]=reps[3]
    # add results to SimulationResults object
    results=SimulationResults()
    results.growth_rates=growth_rates
    results.sim=sim
    results.r2=rs
    results.tops=tops
    results.fracs=fracs
    results.reps = num_replicates
    results.stats = stats

    return (results)

def topology_instant_growth(sim,
                            num_replicates=1,
                            times=[100],
                            s=1,
                            plot_dems=False,
                            mask=0.01,
                            signed=True):
    """simulates replicates of a given simulation, for different (initial) population sizesParams:
        * sim - a Simulation() object of whatever demography you choose
        * num_replicates - number of replicates to run
        * times - list of times of onset of population size change
    Output: a SimulationResults object with the following attributes
        * rs - a (reps x len(pop_sizes) x 3) array with [r^2,f1,f2] where f1 and f2 are the freqs
            of two randomly chosen mutations and r^2 is the r^2 between them
        * tops - a (reps x len(pop_sizes) x 4) tuple of topologies, in order counts of
                (is_above,other_side,same_branch,opposite_branch)
        * fracs - a (reps x len(pop_sizes) x 1) array of the fraction of internal branch length
        * also records the simulation object and the growth rates
    """
    # initialise arrays
    rs = np.zeros(shape=(num_replicates, len(times), 3))
    tops = np.zeros(shape=(num_replicates, len(times), 4))
    fracs = np.zeros(shape=(num_replicates, len(times)))
    sim.s=s

    for i in range(len(times)): # for each of the pop sizes

        T = times[i]
        sim._make_demography(s=s,T=T)
        sim.T = T
        if plot_dems:
            sim._demography_plot()
        reps = sim_replicates_topology(sim, reps=num_replicates,mask=mask,signed=signed)
        # get the values of r2
        r2 = reps[0]
        missing=num_replicates-r2.shape[0]
        filler=np.full([missing, 3], np.nan)
        r2 = np.concatenate((r2, filler), axis=0)
        rs[:, i, :] = r2
        # summarise the tree topologies
        filler=np.full([missing, 4], np.nan)
        tps = np.concatenate((reps[1], filler), axis=0)
        tops[:, i, :] = tps
        # fraction internal branch length
        filler=np.full([missing, 1], np.nan)
        frcs = np.concatenate((reps[2], filler), axis=0)
        fracs[:, i] = frcs.T
    # add results to SimulationResults object
    results=SimulationResults()
    results.times=times
    results.sim=sim
    results.r2=rs
    results.tops=tops
    results.fracs=fracs
    results.reps = num_replicates

    return (results)

def branchwise_afs(ts):
    afs = ts.allele_frequency_spectrum(mode="branch", span_normalise=False)
    return(afs)

def internal_branch_length(ts,fraction=True):
    afs = branchwise_afs(ts)
    extnal=afs[1]
    total=sum(afs)
    intnal= total - extnal # total (sum of all) minus branches which yield singletons
    if fraction:
        return (intnal/total)
    else:
        return(intnal)


def topology_bottleneck(sim,
                         num_replicates=1,
                         btl_times=[100],
                         btl_size=0.5,
                         plot_dems=False,
                         mask=0.01,
                         signed=True):
    """simulates replicates of a given simulation, for different (initial) population sizesParams:
        * sim - a Simulation() object of whatever demography you choose
        * num_replicates - number of replicates to run
        * growth_rates - list of growth rates
    Output: a SimulationResults object with the following attributes
        * rs - a (reps x len(pop_sizes) x 3) array with [r^2,f1,f2] where f1 and f2 are the freqs
            of two randomly chosen mutations and r^2 is the r^2 between them
        * tops - a (reps x len(pop_sizes) x 4) tuple of topologies, in order counts of
                (is_above,other_side,same_branch,opposite_branch)
        * fracs - a (reps x len(pop_sizes) x 1) array of the fraction of internal branch length
        * also records the simulation object and the growth rates
    """
    # initialise arrays
    rs = np.zeros(shape=(num_replicates, len(btl_times), 3))
    tops = np.zeros(shape=(num_replicates, len(btl_times), 4))
    fracs = np.zeros(shape=(num_replicates, len(btl_times)))
    stats={}


    for i in range(len(btl_times)): # for each of the pop sizes

        btl = btl_times[i]
        sim._make_demography(T=btl,s=btl_size)
        if plot_dems:
            sim._demography_plot()
        reps = sim_replicates_topology(sim, reps=num_replicates,
                                       mask=mask,
                                       signed=signed)
        # get the values of r2
        r2 = reps[0]
        missing=num_replicates-r2.shape[0]
        filler=np.full([missing, 3], np.nan)
        r2 = np.concatenate((r2, filler), axis=0)
        rs[:, i, :] = r2
        # summarise the tree topologies
        filler=np.full([missing, 4], np.nan)
        tps = np.concatenate((reps[1], filler), axis=0)
        tops[:, i, :] = tps
        # fraction internal branch length
        filler=np.full([missing, 1], np.nan)
        frcs = np.concatenate((reps[2], filler), axis=0)
        fracs[:, i] = frcs.T
        stats[i]=reps[3]
    # add results to SimulationResults object
    results=SimulationResults()
    results.sim=sim
    results.r2=rs
    results.tops=tops
    results.fracs=fracs
    results.reps = num_replicates
    results.stats=stats

    return (results)



def topology_bottleneck_size(sim,
                         num_replicates=1,
                         btl_time=100,
                         btl_sizes=[0.5],
                         plot_dems=False,
                         mask=0.01,
                         signed=True):
    """simulates replicates of a given simulation, for different (initial) population sizesParams:
        * sim - a Simulation() object of whatever demography you choose
        * num_replicates - number of replicates to run
        * growth_rates - list of growth rates
    Output: a SimulationResults object with the following attributes
        * rs - a (reps x len(pop_sizes) x 3) array with [r^2,f1,f2] where f1 and f2 are the freqs
            of two randomly chosen mutations and r^2 is the r^2 between them
        * tops - a (reps x len(pop_sizes) x 4) tuple of topologies, in order counts of
                (is_above,other_side,same_branch,opposite_branch)
        * fracs - a (reps x len(pop_sizes) x 1) array of the fraction of internal branch length
        * also records the simulation object and the growth rates
    """
    # initialise arrays
    rs = np.zeros(shape=(num_replicates, len(btl_sizes), 3))
    tops = np.zeros(shape=(num_replicates, len(btl_sizes), 4))
    fracs = np.zeros(shape=(num_replicates, len(btl_sizes)))
    stats={}


    for i in range(len(btl_sizes)): # for each of the pop sizes

        btl_size = btl_sizes[i]
        sim._make_demography(T=btl_time,s=btl_size)
        if plot_dems:
            sim._demography_plot()
        reps = sim_replicates_topology(sim, reps=num_replicates,
                                       mask=mask,
                                       signed=signed)
        # get the values of r2
        r2 = reps[0]
        missing=num_replicates-r2.shape[0]
        filler=np.full([missing, 3], np.nan)
        r2 = np.concatenate((r2, filler), axis=0)
        rs[:, i, :] = r2
        # summarise the tree topologies
        filler=np.full([missing, 4], np.nan)
        tps = np.concatenate((reps[1], filler), axis=0)
        tops[:, i, :] = tps
        # fraction internal branch length
        filler=np.full([missing, 1], np.nan)
        frcs = np.concatenate((reps[2], filler), axis=0)
        fracs[:, i] = frcs.T
        stats[i]=reps[3]
        print("done with size",btl_size)
    # add results to SimulationResults object
    results=SimulationResults()
    results.sim=sim
    results.r2=rs
    results.tops=tops
    results.fracs=fracs
    results.reps = num_replicates
    results.stats= stats

    return (results)


def get_mn_tmrca(tree, pop_size):
    """
    This function calculates the mean and variance of branch lengths in a tree, given the population size.

    Parameters:
    tree (tskit.TreeSequence): The tree sequence to analyze.
    pop_size (int): The size of the population.

    Returns:
    list: A list of three values representing the mean, variance, and mean squared over variance of branch lengths.
    """
    nds = list(tree.leaves())
    tmrca = np.zeros(shape=(len(nds), len(nds)))
    tmrca_l = []
    pair_order_list = itertools.combinations(nds, 2)
    for pair in pair_order_list:
        try:
            if pair[1] > pair[0]:
                tmrca[pair[0], pair[1]] = tree.tmrca(pair[0], pair[1]) / (pop_size * 2)
                tmrca_l.append(tree.tmrca(pair[0], pair[1]) / (pop_size * 2))
        except:
            pass
    return ([np.mean(tmrca_l), np.var(tmrca_l), np.mean(tmrca_l) ** 2 / np.var(tmrca_l)])


def get_r_matrix(ts):
    """
    This function calculates the signed LD matrix for a tree sequence.

    Parameters:
    ts (tskit.TreeSequence): The tree sequence to analyze.

    Returns:
    numpy.ndarray: A matrix of Pearson's correlation coefficients representing the signed LD between sites.
    """
    unique_mutations = np.unique(list(ts.tables.mutations.site), return_index=True)[1]
    derived = ts.tables.mutations.derived_state[unique_mutations]  # getting the derived allele
    # make a haplotype matrix by extracting haplotypes and comparing to derived allele
    haps = np.array([np.array(list(map(ord, hap))) == derived for hap in list(ts.haplotypes())])
    # get a signed LD matrix
    cor_matrix = np.corrcoef(haps.T, rowvar=True)
    return (cor_matrix)

def get_afs_from_matrix(ts):
    """
    This function calculates the signed LD matrix for a tree sequence.

    Parameters:
    ts (tskit.TreeSequence): The tree sequence to analyze.

    Returns:
    numpy.ndarray: A matrix of Pearson's correlation coefficients representing the signed LD between sites.
    """
    unique_mutations = np.unique(list(ts.tables.mutations.site), return_index=True)[1]
    derived = ts.tables.mutations.derived_state[unique_mutations]  # getting the derived allele
    # make a haplotype matrix by extracting haplotypes and comparing to derived allele
    haps = np.array([np.array(list(map(ord, hap))) == derived for hap in list(ts.haplotypes())])
    # get a signed LD matrix
    afs = np.sum(haps.T,axis=1)/haps.T.shape[1]
    return (afs)


def get_configuration_probabilities_old(tree):
    # returns the probabilities of different configurations of mutations, for a given tree.
    prob_same_branch=0
    prob_other_side=0
    total = tree.total_branch_length
    # iterate through nodes
    for u in tree.nodes():
        pathlength_abv=0
        pathlength_blw=0
        branchlength=0
        v = u
        # get path to root
        path_to_root = []
        while v != tskit.NULL:
            if v!=u:
                path_to_root.append(v)
            v = tree.parent(v)
        # get list of descendants
        descendants=[]
        for node in tree.nodes():
            if tree.is_descendant(node,u) and node!=u: # check whether the node is a descendant of v
                descendants.append(node)
        # get length of current branch
        branchlength=tree.branch_length(u)
        # get path length below
        for descendant in descendants:
            pathlength_blw+=tree.branch_length(descendant)
        # get path length above
        for ancestor in path_to_root:
            pathlength_abv+=tree.branch_length(ancestor)
        prob_same_branch+=(branchlength/total)**2
        prob_other_side+=(1-((pathlength_abv+pathlength_blw+branchlength)/total))*(branchlength/total)
    return([prob_same_branch,prob_other_side,1-(prob_same_branch+prob_other_side)])


def get_configuration_probabilities(tree):
    """
    Returns the probabilities of different configurations of mutations, for a given tree.

    Parameters:
    tree (tskit tree object): Tree object from which to extract information

    Returns:
    list: A list of 3 probabilities, for each of the following configurations:
        - Probability of mutations being on the same branch
        - Probability of mutations being on different branches, but on the same side of the tree
        - Probability of mutations being on different branches, on different sides of the tree

    """

    # Initialize variables to store the probabilities
    prob_same_branch = 0
    prob_other_side = 0
    prob_same_side = 0
    total = tree.total_branch_length

    # Iterate through all the nodes in the tree
    for u in tree.nodes():
        # Initialize path length variables
        pathlength_abv = 0
        pathlength_blw = 0
        branchlength = 0

        # Get the path to root for node u
        v = u
        path_to_root = []
        while v != tskit.NULL:
            if v != u:
                path_to_root.append(v)
            v = tree.parent(v)

        # Get the descendants of node u
        descendants = []
        for node in tree.nodes():
            if tree.is_descendant(node, u) and node != u:
                descendants.append(node)

        # Get the length of the current branch
        branchlength = tree.branch_length(u)

        # Get the path length below node u
        for descendant in descendants:
            pathlength_blw += tree.branch_length(descendant)

        # Get the path length above node u
        for ancestor in path_to_root:
            pathlength_abv += tree.branch_length(ancestor)

        # Increment the appropriate probabilities based on the path lengths
        prob_same_branch += (branchlength / total) ** 2
        prob_same_side += 2 * (pathlength_blw / total) * (branchlength / total)
        prob_other_side += (1 - ((pathlength_abv + pathlength_blw + branchlength) / total)) * (branchlength / total)

    # Return the list of probabilities
    return ([prob_same_branch, prob_other_side, prob_same_side])


def stats_from_trees(ts, pop_size):
    """
    Calculate statistics from a set of trees.

    Parameters:
        ts (TreeSequence): The input set of trees.
        pop_size (int): The population size.

    Returns:
        dict: A dictionary containing statistics for correlations, time, and configurations.

    """
    # Get the first tree from the set of trees
    tree = ts.first()

    # Calculate the configuration probabilities
    c_stats = get_configuration_probabilities(tree)

    # Calculate the correlation matrix
    cor_matrix = get_r_matrix(ts)

    # Calculate statistics for perfect linkage disequilibrium (r=1),
    # negative correlations (r<0), and positive correlations (r>0)
    r_stats = [
        np.sum(cor_matrix == 1) / (cor_matrix.shape[0] ** 2),
        1 - np.sum(cor_matrix > 0) / (cor_matrix.shape[0] ** 2),
        np.sum(cor_matrix > 0) / (cor_matrix.shape[0] ** 2)
    ]

    # Calculate the mean time to the most recent common ancestor
    t_stats = get_mn_tmrca(tree, pop_size)

    # Return a dictionary containing statistics for correlations, time, and configurations
    return {"cor": r_stats, "t": t_stats, "c": c_stats}



def topology_pop_size(sim, num_replicates=1, pop_sizes=[10_000],
                      mask=0.01,
                      signed=True):
    """simulates replicates of a given simulation, for different (initial) population sizesParams:
        * sim - a Simulation() object of whatever demography you choose
        * num_replicates - number of replicates to run
        * pop_sizes - list of initial population sizes
    Output: a SimulationResults object with the following attributes
        * rs - a (reps x len(pop_sizes) x 3) array with [r^2,f1,f2] where f1 and f2 are the freqs
            of two randomly chosen mutations and r^2 is the r^2 between them
        * tops - a (reps x len(pop_sizes) x 4) tuple of topologies, in order counts of
                (is_above,other_side,same_branch,opposite_branch)
        * fracs - a (reps x len(pop_sizes) x 1) array of the fraction of internal branch length
        * sim and pop_sizes
        * stats - a dictionary containing:
            + 'c' probabilities of mutation configurations, as calculated from branch lengths
            + 't' mean, variance and mean^2/var coalescent times for all pairs of samples
            + 'cor' the MC estimates of mutation configurations from the r matrix
    """
    # initialise arrays
    rs = np.zeros(shape=(num_replicates, len(pop_sizes), 3))
    tops = np.zeros(shape=(num_replicates, len(pop_sizes), 4))
    fracs = np.zeros(shape=(num_replicates, len(pop_sizes)))

    for i in range(len(pop_sizes)): # for each of the pop sizes
        ps = pop_sizes[i]
        sim.Ne = ps
        reps = sim_replicates_topology(sim, reps=num_replicates,mask=mask,signed=signed,print_errors=True)
        # get the values of r2
        r2 = reps[0]
        missing=num_replicates-r2.shape[0]
        filler=np.full([missing, 3], np.nan)
        r2 = np.concatenate((r2, filler), axis=0)
        rs[:, i, :] = r2
        # summarise the tree topologies
        filler=np.full([missing, 4], np.nan)
        tps = np.concatenate((reps[1], filler), axis=0)
        tops[:, i, :] = tps
        # fraction internal branch length
        filler=np.full([missing, 1], np.nan)
        frcs = np.concatenate((reps[2], filler), axis=0)
        fracs[:, i] = frcs.T

    results = SimulationResults()
    results.pop_sizes = pop_sizes
    results.sim = sim
    results.r2 = rs
    results.tops = tops
    results.fracs = fracs
    results.reps = num_replicates

    return (results)









