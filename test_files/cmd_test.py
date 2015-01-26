#!/usr/bin/env python
from random import randint, uniform

from pep_ga.ga import GeneticAlgorithm
from pep_ga.options import selection_opts, fitness_opts, recombination_opts


def test_methods():
    for selection in selection_opts.values():
        for recombination in recombination_opts.values():
            for fitness in fitness_opts.values():
                runner = GeneticAlgorithm(fitness_function=fitness,
                                          recombination_function=recombination,
                                          selection_function=selection,
                                          generations=50,
                                          mutation_rate=0.05,
                                          peptide_length=randint(8, 20),
                                          population_size=100,
                                          survivors_per_round=20)
                print "SELECTION METHOD: %s" % selection
                print "RECOMBINATION METHOD: %s" % recombination
                print "FITNESS FUNCTION: %s" % fitness
                runner.run()
                runner.log_status()

if __name__ == '__main__':
    test_methods()


