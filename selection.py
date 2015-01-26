#!/usr/bin/env python
import random

#Selection Functions
def proportional_selection(population):
    '''
    Proportional (roulette wheel) selection: To each
    individual an area on a roulette wheel is assigned
    depending on its fitness value. Individuals with
    higher value have a higher probability to be selected
    than individuals with lower fitness value.
    '''
    return [stochastic_universal_sampling(population, 1)[0]]
            
def stochastic_universal_sampling(population, number_to_keep=4):
    '''
    Similar to the proportional selection,
    every individual obtains a segment on a
    roulette wheel according to its fitness
    value.
    However, it is turned only one time with
    nballs where n is the number of individuals in the
    population.
    https://en.wikipedia.org/wiki/Stochastic_universal_sampling
    '''
    sum_fitnesses = sum([c.fitness for c in population])
    pointer_spread = float(sum_fitnesses) / number_to_keep
    start = random.uniform(0, pointer_spread)
    pointers = [(start+i*pointer_spread) for i in range(number_to_keep)]
    return _roulette_wheel(population, pointers)
            
def _roulette_wheel(population, pointers):
    keep = []
    i = 0
    cum_fitness = population[i].fitness
    for point in pointers:
        while cum_fitness < point:
            i += 1
            cum_fitness += population[i].fitness
        keep.append(population[i])
    #Picks should be unique!
    assert len(keep) == len(set(keep))
    return keep

def linear_rank_selection(population, num_to_return=1):
    '''
    The individuals obtain a
    rank in correspondence to their fitness value. The
    selection is performed on the basis of this rank.
    '''
    assert num_to_return <= len(population)
    return sort_population_by_fitness(population)[:num_to_return]

def binary_tournament_selection(population):
    '''
    Two randomly
    selected individuals compete for their selection
    where the one with the higher fitness value wins
    without any other stochastic influences
    '''
    c1 = random.choice(population)
    c2 = _select_from_excluding(population, [c1])
    selected = c1 if c1.fitness > c2.fitness else c2
    return [selected]

def _select_from_excluding(sequence, excluded):
    ''' Select a random item from a sequence not in "excluded".'''
    assert len(excluded) < len(sequence)
    selection = random.choice(sequence)
    while selection in excluded:
        selection = random.choice(sequence)
    return selection
    
def random_selection(population):
    return [random.choice(population)]
    
def best_fraction_selection(population, n=0.05):
    '''
    The n best fraction of the
    population are chosen straight forwardly.
    '''
    number_to_return = max(int(round(len(population) * n)), 1)
    return sort_population_by_fitness(population)[:number_to_return]
    
def q_tournament_selection(population, q=30, num_to_return=1):
    '''
    All individuals participate
    in q tournaments, where the individuals with the
    most victories are selected.
    import pdb; pdb.set_trace()
    '''
    assert num_to_return <= len(population)
    winner = None
    for chromosome in population:
        wins = 0
        for _ in range(q):
            opponent = _select_from_excluding(population, [chromosome])
            if chromosome.fitness > opponent.fitness:
                wins += 1
        if not winner:
            winner = (chromosome, wins)
        if wins > winner[1]:
            winner = (chromosome, wins)
    return winner[:num_to_return]
    
def sort_population_by_fitness(population):
    return sorted(population, key = lambda c: c.fitness, reverse=True)