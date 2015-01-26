from functools import partial

import recombination
import selection
from fitness import trivial_fitness, vina_fitness, multi_vina_fitness
from chromosomes import fragment_scaffold
from chromosomes.sequence import GAPeptideChromosome, \
    GAHeadTailPeptideChromosome, GADisulphidePeptideChromosome


recombination_opts = {"SP": recombination.single_point_crossover,
                      "DP": recombination.double_point_crossover,
                      "DIST-BIS": recombination.distance_bisector_crossover,
                      "MP": recombination.multi_point_crossover,
                      "UNI": recombination.uniform_crossover,
                      "SHUFFLE": recombination.shuffle_crossover}

selection_opts = {"PROPORTIONAL": selection.proportional_selection,
                  "LIN-RANK": selection.linear_rank_selection,
                  "BIN-TNMT": selection.binary_tournament_selection,
                  "RANDOM": selection.random_selection,
                  "BEST-FRAC": selection.best_fraction_selection,
                  "Q-TNMT": selection.q_tournament_selection,
                  "STOCH-UNI": selection.stochastic_universal_sampling}

fitness_opts = {'':trivial_fitness.trivial_peptide_fitness,
                'vina_fitness':vina_fitness.VinaFitness,
                'multi_fitness':multi_vina_fitness.MultiVinaFitness}

chromosome_opts = {'peptide': {'':GAPeptideChromosome,
                               'HT':GAHeadTailPeptideChromosome,
                               'SS':GADisulphidePeptideChromosome},
                   'scaffold': fragment_scaffold.scaffold_factory()}