#!/usr/bin/env python

MIN_SCORE = 0.01
#All fitnesses need to be > 0 and positive for selection functions to work!


def trivial_peptide_fitness(chromosomes):
    return [c.sequence.count('A') + MIN_SCORE for c in chromosomes]


