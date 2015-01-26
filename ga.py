#!/usr/bin/env python
import random
import logging
import datetime
import cPickle
from time import time


class GeneticAlgorithm(object):
    """Runs genetic algorithm, contains chromosomes and applies fitness,
    selection, recombination and mutation methods to them."""
    def __init__(self, 
                 population_size,
                 survivors_per_round,
                 generations,
                 mutation_rate,
                 peptide_length,
                 fitness_function,
                 selection_function,
                 recombination_function,
                 chromosome_type,
                 chromosome_elements):
        self.population_size = population_size
        self.survivors_per_round = survivors_per_round
        self.generations = generations
        self.mutation_rate = mutation_rate
        self.chromosome_length = peptide_length
        self.current_generation = 0
        self.population = None
        self.selection_func = selection_function
        self.fitness_function = fitness_function
        self.recombination_function = recombination_function
        self.chromosome_type = chromosome_type
        self.elements = chromosome_elements
        self.run_record = []
        self.do_logging = True
        self._init_time = time()
        self._setup_ga()
        
    def run(self, do_logging=True):
        if do_logging:
            self.do_logging = True
            logging.info("-----RUNNING-GA-----")
        while self.current_generation < self.generations:
            self._run_selection_round()
        if self.do_logging:
            logging.info("-----FINISHED------")
            self.log_status()
        self.final_summary()
            
    def log_status(self):
        ave_fitness = lambda:sum(c.fitness for c in
                                 self.population)/self.population_size
        logging.info("\n\n..........")
        logging.info("Population size: %d" % self.population_size)
        logging.info("Survivors per round: %d" % self.survivors_per_round)
        logging.info("Current round %d/%d" % (self.current_generation,
                                              self.generations))
        time_elapsed = datetime.timedelta(seconds=round(time()-self._init_time))
        logging.info("Time elapsed %s" % str(time_elapsed))
        logging.info("Average Fitness %f" % ave_fitness())
        logging.info("..........\n\n")

    @classmethod
    def load_run(cls, f):
        with open(f) as data:
            return cPickle.load(data)

    def save(self):
        pickle_fname = "round%d.pkl" % self.current_generation
        with open(pickle_fname, 'wb') as f:
            cPickle.dump(self, f)
        
    def _setup_ga(self):
        assert self.survivors_per_round >= 2
        assert self.survivors_per_round < self.population_size
        assert 0.0 < self.mutation_rate < 1.0
        if self.do_logging:
            logging.info("\n\n\n-----BEGINNING-GA-RUN-----")
            logging.info("-----CALCULATING-INITIAL-GA-FITNESSES-----")
        self._random_initial_population()
        
    def _random_initial_population(self):
        random_chromosome = self.chromosome_type.random_chromosome
        self.population = []
        while len(self.population) < self.population_size:
            new_chromosome = random_chromosome(self.chromosome_length,
                                               self.mutation_rate,
                                               elements=self.elements)
            assert new_chromosome, "Failed to construct chromosome"
            if new_chromosome not in self.population:
                self.population.append(new_chromosome)
        self._calc_population_fitness()
        initial_statistics = {'initial_population': [(c.idx_name, c.fitness)
                                                     for c in self.population]}
        self.run_record.append(initial_statistics)
    
    def _calc_population_fitness(self):
        fitnesses = self.fitness_function(self.population)
        for chromosome, fitness in zip(self.population, fitnesses):
            #Fitness function modifies chromosome fitness in place
            chromosome.fitness = fitness
            
    def _run_selection_round(self):
        round_statistics = {'selected': [], 'children': []}
        seqs_as_str = lambda l: ", ".join(c.idx_name for c in l)
        if self.do_logging:
            self.log_status()
            logging.info("-----RUNNING-SELECTION-ROUND-----")
        self.current_generation += 1
        selected = []
        #Select survivors
        while len(selected) < self.survivors_per_round:
            selected_chromosomes = self.selection_func(self.population)
            assert [s in self.population for s in selected_chromosomes]
            selected.extend(selected_chromosomes)
            #Delete chosen chromosomes from selection pool (i.e. population)
            for selected_chromosome in selected_chromosomes:
                self.population.remove(selected_chromosome)
        if self.do_logging:
            logging.info("SELECTED CHROMOSOMES: %s" % seqs_as_str(selected))
        round_statistics['selected'] = [(c.idx_name, c.fitness)
                                        for c in selected]
        #Recombination of survivors to make children
        #Randomly select pairs of parents until there are enough children
        children = []
        while len(children) < self.population_size:
            p1, p2 = _random_choose_two(selected)
            new_sequence = self.recombination_function(p1.sequence, p2.sequence)
            #Random mutation of children
            child = self.chromosome_type(new_sequence, self.mutation_rate,
                                         elements=self.elements)
            child.mutate()
            #Do not add duplicate sequences
            if child.sequence not in [c.sequence for c in children]:
                children.append(child)
        if self.do_logging:
            logging.info("NEXT GENERATION: %s" % seqs_as_str(children))
        num_unique_sequences = len(set(c.sequence for c in children))
        assert num_unique_sequences == self.population_size
        self.population = children
        self._calc_population_fitness()
        round_statistics['children'] = [(c.idx_name, c.fitness)
                                        for c in self.population]
        self.run_record.append(round_statistics)
        self.save()

    def final_summary(self):
        with open('ga-summary.txt', 'w') as out:
            line = "Name: %s, Fitness %s"
            print>>out, "Initial Population:"
            for (name, fitness) in self.run_record[0]['initial_population']:
                print>>out, line % (name, fitness)
            print>>out, '\n'
            for i, ga_round in enumerate(self.run_record[1:]):
                print>>out, "Round %d" % (i+1)
                for step in ('selected', 'children'):
                    print>>out, step
                    for (name, fitness) in ga_round[step]:
                        print>>out, line % (name, fitness)
                print>>out


def _random_choose_two(seq):
    a = random.choice(seq)
    b = random.choice(seq)
    while b is not a:
        b = random.choice(seq)
    return a, b
