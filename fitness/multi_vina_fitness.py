import os
import re
import logging
import multiprocessing

from vina_fitness import vina_abs_docking_score, pdbqt_from_smiles, MIN_SCORE


class MultiVinaFitness(object):
    """Calculates a Vina fitness across receptors according to an equation.
    This class allows storing all of the config and receptor information
    on setting up the instance, to allow the fitness to be calculated by
    passing only the list of chromosomes.

    This class can also be pickled, in constrast with alternative approach
    of returning nested functions.
    """

    def __init__(self, fitness_equation, keep_outfile,
                 parallel=False, *receptor_config_pairs):
        self.parallel = parallel
        self.fitness_equation = fitness_equation
        self.keep_outfile  = keep_outfile
        self.receptor_config_pairs = receptor_config_pairs

    def __call__(self, chromosomes):
        fitnesses = []
        if not self.parallel:
            for chromosome in chromosomes:
                fitnesses.append(multiple_fitness(chromosome,
                                                  self.fitness_equation,
                                                  self.receptor_config_pairs,
                                                  self.keep_outfile,
                                                  cpu_count=0))
        else:
            pool = multiprocessing.Pool(multiprocessing.cpu_count())
            arg_supplier = [(c, self.fitness_equation,
                             self.receptor_config_pairs,
                             self.keep_outfile, 1) for c in chromosomes]
            fitnesses = list(pool.map(_starred_multiple_fitness, arg_supplier))
        return fitnesses


def _starred_multiple_fitness(args):
    return multiple_fitness(*args)


def multiple_fitness(chromosome, fitness_equation, receptor_config_pairs,
                     keep_outfile=True, cpu_count=0):
    """
    Calculates a chromosome's vina fitness for multiple receptors.
    Returns a list of fitnesses.
    """
    lig_f = pdbqt_from_smiles(chromosome.smiles, chromosome.idx_name)
    try:
        fitnesses = [vina_abs_docking_score(lig_f, rec_f,
                                            conf_f, cpu_count, keep_outfile)
                     for (rec_f, conf_f) in receptor_config_pairs]
        fitness = max(MIN_SCORE, calc_multi_vina_fitness(fitnesses,
                                                         fitness_equation))
        logging.info("Calculated chromosome multi-fitness:"
                     " Sequence %s, Fitness %s"
                     % (chromosome.idx_name, fitness))
        return fitness
    finally:
        try:
            os.remove(lig_f)
        except OSError:
            pass
    return fitness


def calc_multi_vina_fitness(receptor_fitnesses, fitness_equation):
    """ Evaluates fitness from an equation string:
    Equation string will refer to fitnesses in the form r_0, r_1, etc
    Fitness equation can only contain limited characters: +-*/()r_0123456789
    """
    allowed_chars = '+-*/()r_0123456789 '
    assert (len([c for c in fitness_equation if c in allowed_chars])
            == len(fitness_equation))
    search_patt = r"r_([0-9]+)"
    replace_patt = r"r[\1]"
    fitness_equation = re.sub(search_patt, replace_patt, fitness_equation)
    r = receptor_fitnesses
    return eval(fitness_equation)