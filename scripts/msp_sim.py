"""Functions for simulations - borrowed from A. Biddanda
    https://github.com/aabiddanda/aDNA_LD_public/blob/master/src/aDNA_coal_sim.py"""

import numpy as np
import msprime as msp
import tskit
import matplotlib.pyplot as plt
import stdpopsim

class Simulation(object):

    def __init__(self):
        """Create a simulation object."""
        # Defining effective population sizes
        self.Ne = None

        # Defining samples in this case
        self.samples = 100

        # Define demographic events to sample across
        self.demography = None

        # Define demographic events to sample across
        self.mutation_rate = 1e-7


    def _simulate(self, samples,**kwargs):
        """Conduct a simulation using msprime and the parameters we have."""
        # Generate a tree sequence
        tree_seq = msp.sim_ancestry(
            population_size=self.Ne,
            samples=samples,
            sequence_length=1_000,
            recombination_rate=0,
            demography=self.demography,
            **kwargs
        )
        tree_seq = msp.sim_mutations(tree_seq,rate=self.mutation_rate)
        return tree_seq

    def _sim_replicates(self,samples, reps,**kwargs):
        ancestry_reps = msp.sim_ancestry(
            population_size=self.Ne,
            samples=samples,
            sequence_length=1_000,
            recombination_rate=0,
            demography=self.demography,
            num_replicates=reps,
            **kwargs
        )
        for ts in ancestry_reps:
            mutated_ts = msp.sim_mutations(ts, rate=self.mutation_rate)
            yield mutated_ts


class ExpGrowth(Simulation):
    """Simple model of exponential growth."""

    def __init__(self, r=0.1,Ne=10_000):
        """Initialize the model.
        Params:
            * r - exponential growth rate
        """
        super().__init__()
        self.Ne = Ne

        dem=msp.Demography()

        dem.add_population(initial_size=Ne,
                           growth_rate=r)

        self.demography = dem

        self.mutation_rate = None

    def _make_demography(self,r=0):
        Ne=self.Ne
        dem = msp.Demography()

        dem.add_population(initial_size=Ne,
                           growth_rate=r)

        self.demography = dem


    def _simulate(self, samples,**kwargs):
        """Simulate a panel of individuals under the population growth model."""
        assert self.samples is not None
        # Generate a tree sequence
        tree_seq = msp.sim_ancestry(
            samples=samples,
            sequence_length=1_000,
            recombination_rate=0,
            demography=self.demography,
            **kwargs
        )
        tree_seq = msp.sim_mutations(tree_seq,
                                     rate=self.mutation_rate)
        return tree_seq

    def _sim_replicates(self,samples, reps,**kwargs):
        ancestry_reps = msp.sim_ancestry(
            samples=samples,
            sequence_length=1_000,
            recombination_rate=0,
            demography=self.demography,
            num_replicates=reps,
            **kwargs
        )
        for replicate_index, ts in enumerate(ancestry_reps):
            mutated_ts = msp.sim_mutations(ts, rate=self.mutation_rate)
            yield mutated_ts

    def _demography_plot(self):
        """Demography plotting."""
        debug = self.demography.debug()
        t = np.linspace(0, 200, num=20)
        S = debug.population_size_trajectory(t)
        plt.plot(t, S)
        plt.xlabel("Time ago")
        plt.ylabel("Population size");



class InstantGrowth(Simulation):
    """Simple model of constant growth."""

    def __init__(self, s=0.5,Ne=10_000,T=100):
        """Initialize the model.
        Params:
            * T - time of growth onset
            * s - relative size of population before
            * Ne - size of population at t=0 (in the present)
        """
        super().__init__()
        self.Ne = Ne
        self.Ne_before = Ne * s
        self.T=T

        dem=msp.Demography()

        dem.add_population(initial_size=Ne)
        dem.add_population_parameters_change(time=T,
                                             initial_size=self.Ne_before)

        self.demography = dem

        self.mutation_rate = None

    def _make_demography(self,s=0.5,T=100):
        self.Ne_before = self.Ne * s
        self.T=T
        dem=msp.Demography()

        dem.add_population(initial_size=self.Ne)
        dem.add_population_parameters_change(time=T,
                                             initial_size=self.Ne_before)
        self.demography = dem


    def _simulate(self, samples,**kwargs):
        """Simulate a panel of individuals under the population growth model."""
        assert self.samples is not None
        # Generate a tree sequence
        tree_seq = msp.sim_ancestry(
            samples=samples,
            sequence_length=1_000,
            recombination_rate=0,
            demography=self.demography,
            **kwargs
        )
        tree_seq = msp.sim_mutations(tree_seq,
                                     rate=self.mutation_rate)
        return tree_seq

    def _sim_replicates(self,samples, reps,**kwargs):
        ancestry_reps = msp.sim_ancestry(
            samples=samples,
            sequence_length=1_000,
            recombination_rate=0,
            demography=self.demography,
            num_replicates=reps,
            **kwargs
        )
        for replicate_index, ts in enumerate(ancestry_reps):
            mutated_ts = msp.sim_mutations(ts, rate=self.mutation_rate)
            yield mutated_ts

    def _demography_plot(self):
        """Demography plotting."""
        debug = self.demography.debug()
        t = np.linspace(0, 200, num=20)
        S = debug.population_size_trajectory(t)
        plt.plot(t, S)
        plt.xlabel("Time ago")
        plt.ylabel("Population size");


class Bottleneck(Simulation):
    """Simple model of constant growth."""

    def __init__(self, s=0.5,Ne=10_000,T=100):
        """Initialize the model.
        Params:
            * T - time of growth onset
            * s - relative size of population before
            * Ne - size of population at t=0 (in the present)
        """
        super().__init__()
        self.Ne = Ne
        self.Ne_before = Ne * s
        self.T = T

        dem=msp.Demography()

        dem.add_population(initial_size=Ne,name="P")
        dem.add_simple_bottleneck(time=T,
                                  proportion=s,
                                  population="P")

        self.demography = dem

        self.mutation_rate = None

    def _make_demography(self, T=100,s=0.5):
        Ne = self.Ne
        dem=msp.Demography()

        dem.add_population(initial_size=Ne,name="P")
        dem.add_simple_bottleneck(time=T,
                                  proportion=s,
                                  population="P")

        self.demography = dem


    def _simulate(self, samples,**kwargs):
        """Simulate a panel of individuals under the population growth model."""
        assert self.samples is not None
        # Generate a tree sequence
        tree_seq = msp.sim_ancestry(
            samples=samples,
            sequence_length=1_000,
            recombination_rate=0,
            demography=self.demography,
            **kwargs
        )
        tree_seq = msp.sim_mutations(tree_seq,
                                     rate=self.mutation_rate)
        return tree_seq

    def _sim_replicates(self,samples, reps,**kwargs):
        ancestry_reps = msp.sim_ancestry(
            samples=samples,
            sequence_length=1_000,
            recombination_rate=0,
            demography=self.demography,
            num_replicates=reps,
            **kwargs
        )
        for replicate_index, ts in enumerate(ancestry_reps):
            mutated_ts = msp.sim_mutations(ts, rate=self.mutation_rate)
            yield mutated_ts

    def _demography_plot(self):
        """Demography plotting."""
        debug = self.demography.debug()
        t = np.linspace(0, 200, num=20)
        S = debug.population_size_trajectory(t)
        plt.plot(t, S)
        plt.xlabel("Time ago")
        plt.ylabel("Population size");


class SimulationResults(object):

    def __init__(self):
        """Create a results object."""

        self.topologies=None

        self.r2=None

        self.fracs=None

        self.sim=None

        self.growth_rates=None

        self.pop_sizes=None

        self.num_failed=None

        self.reps=None

        self.times=None

    def _get_failed(self,**kwargs):
        """get number of failed attempts for each simulation situation"""
        missing=self.reps-self.r2.shape[1]
        self.num_failed=missing





def _ooa_2():
    id = "OutOfAfrica_2T12"
    description = "Two population out-of-Africa"
    long_description = """
        The model is derived from the Tennesen et al. analysis of the
        jSFS from European Americans and African Americans.
        It describes the ancestral human population in Africa, the out of Africa event,
        and two distinct periods of subsequent European population growth over the past
        23kya. Model parameters are taken from Fig. S5 in Fu et al.
    """
    populations = [
        stdpopsim.Population(id="AFR", description="African Americans"),
        stdpopsim.Population(id="EUR", description="European Americans")
    ]
    citations = [
        _tennessen_et_al,
        stdpopsim.Citation(
            author="Fu et al.",
            year=2013,
            doi="https://doi.org/10.1038/nature11690",
            reasons={stdpopsim.CiteReason.DEM_MODEL})
    ]

    generation_time = 25

    T_AF = 148e3 / generation_time
    T_OOA = 51e3 / generation_time
    T_EU0 = 23e3 / generation_time
    T_EG = 5115 / generation_time

    # Growth rates
    r_EU0 = 0.00307
    r_EU = 0.0195
    r_AF = 0.0166

    # population sizes
    N_A = 7310
    N_AF1 = 14474
    N_B = 1861
    N_EU0 = 1032
    N_EU1 = N_EU0 / math.exp(-r_EU0 * (T_EU0-T_EG))

    # migration rates
    m_AF_B = 15e-5
    m_AF_EU = 2.5e-5

    # present Ne
    N_EU = N_EU1 / math.exp(-r_EU * T_EG)
    N_AF = N_AF1 / math.exp(-r_AF * T_EG)

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_AF, growth_rate=r_AF,
                metadata=populations[0].asdict()),
            msprime.PopulationConfiguration(
                initial_size=N_EU, growth_rate=r_EU,
                metadata=populations[1].asdict())
        ],
        migration_matrix=[
            [0, m_AF_EU],
            [m_AF_EU, 0],
        ],
        demographic_events=[
            msprime.MigrationRateChange(
                time=T_EG, rate=m_AF_EU, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_EG, rate=m_AF_EU, matrix_index=(1, 0)),
            msprime.PopulationParametersChange(
                time=T_EG, growth_rate=r_EU0, initial_size=N_EU1, population_id=1),
            msprime.PopulationParametersChange(
                time=T_EG, growth_rate=0, initial_size=N_AF1, population_id=0),
            msprime.MigrationRateChange(
                time=T_EU0, rate=m_AF_B, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_EU0, rate=m_AF_B, matrix_index=(1, 0)),
            msprime.PopulationParametersChange(
                time=T_EU0, initial_size=N_B, growth_rate=0, population_id=1),
            msprime.MassMigration(
                time=T_OOA, source=1, destination=0, proportion=1.0),
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0)
        ],
        )


class Tennessen(Simulation):
    """Tennessen et al European Demography - stdPopSim Consortium
    https://github.com/popsim-consortium/stdpopsim/blob/c8b557bbfb38ad4371a818dc30acf0e65f15e182/stdpopsim/catalog/homo_sapiens.py#L304"""

    def __init__(self, scale=1.0):
        """Initialize the model with fixed parameters."""
        super().__init__()
        self.Ne = 1e4 * scale
        generation_time = 25
        T_AF = 148e3 / generation_time
        T_OOA = 51e3 / generation_time
        T_EU0 = 23e3 / generation_time
        T_EG = 5115 / generation_time

        # Growth rates
        r_EU0 = 0.00307
        r_EU = 0.0195
        r_AF = 0.0166

        # population sizes
        N_A = 7310 * scale
        N_AF1 = 14474 * scale
        N_B = 1861 * scale
        N_EU0 = 1032 * scale
        N_EU1 = N_EU0 / np.exp(-r_EU0 * (T_EU0 - T_EG))

        # migration rates
        m_AF_B = 15e-5
        m_AF_EU = 2.5e-5

        # present Ne
        N_EU = N_EU1 / np.exp(-r_EU * T_EG)
        N_AF = N_AF1 / np.exp(-r_AF * T_EG)

        population_configurations = [
            msp.PopulationConfiguration(initial_size=N_AF, growth_rate=r_AF),
            msp.PopulationConfiguration(initial_size=N_EU, growth_rate=r_EU),
        ]
        migration_matrix = [[0, m_AF_EU], [m_AF_EU, 0]]
        demographic_events = [
            msp.MigrationRateChange(time=T_EG, rate=m_AF_EU, matrix_index=(0, 1)),
            msp.MigrationRateChange(time=T_EG, rate=m_AF_EU, matrix_index=(1, 0)),
            msp.PopulationParametersChange(
                time=T_EG, growth_rate=r_EU0, initial_size=N_EU1, population_id=1
            ),
            msp.PopulationParametersChange(
                time=T_EG, growth_rate=0, initial_size=N_AF1, population_id=0
            ),
            msp.MigrationRateChange(time=T_EU0, rate=m_AF_B, matrix_index=(0, 1)),
            msp.MigrationRateChange(time=T_EU0, rate=m_AF_B, matrix_index=(1, 0)),
            msp.PopulationParametersChange(
                time=T_EU0, initial_size=N_B, growth_rate=0, population_id=1
            ),
            msp.MassMigration(time=T_OOA, source=1, destination=0, proportion=1.0),
            msp.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0
            ),
        ]

        self.pop_config = population_configurations
        self.migration_matrix = migration_matrix
        self.demography = demographic_events

    def _add_samples(self, samples=100):
        """Add in the samples."""
        samples = [msp.Sample(0, 0) for i in range(samples)]
        # Set samples as the parameters for the object
        self.samples = samples

    def _simulate(self,samples=100, **kwargs):
        """Simulate a panel of individuals under the Tennessen et al European Model."""
        self._add_samples(samples=samples)
        assert self.samples is not None
        # Generate a tree sequence
        tree_seq = msp.simulate(
            samples=self.samples,
            demographic_events=self.demography,
            population_configurations=self.pop_config,
            migration_matrix=self.migration_matrix,
            **kwargs
        )
        tree_seq = msp.sim_mutations(tree_seq,
                                     rate=self.mutation_rate)
        return tree_seq

    def _sim_replicates(self,samples=100, reps=1,**kwargs):
        self._add_samples(samples=samples)
        ancestry_reps = msp.sim_ancestry(
            samples=samples,
            sequence_length=1_000,
            recombination_rate=0,
            demography=self.demography,
            num_replicates=reps,
            **kwargs
        )
        for replicate_index, ts in enumerate(ancestry_reps):
            mutated_ts = msp.sim_mutations(ts, rate=self.mutation_rate)
            yield mutated_ts

    def _demography_plot(self):
        """Demography plotting."""
        dem=msp.Demography(self.demography)
        debug = dem.debug()
        t = np.linspace(0, 200, num=20)
        S = debug.population_size_trajectory(t)
        plt.plot(t, S)
        plt.xlabel("Time ago")
        plt.ylabel("Population size");


class SimulationResults(object):

    def __init__(self):
        """Create a results object."""

        self.topologies=None

        self.r2=None

        self.fracs=None

        self.sim=None

        self.growth_rates=None

        self.pop_sizes=None

        self.num_failed=None

        self.reps=None

        self.times=None

    def _get_failed(self,**kwargs):
        """get number of failed attempts for each simulation situation"""
        missing=self.reps-self.r2.shape[1]
        self.num_failed=missing














