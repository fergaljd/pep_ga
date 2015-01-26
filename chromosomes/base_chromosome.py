import random

def basic_mutation(sequence, mutation_rate, building_blocks):
    """
    An amino acid is randomly mutated into any other amino acid.
    """
    mutated_sequence = []
    for element in sequence:
        if random.random() < mutation_rate:
            new_element = random.choice(building_blocks)
            while new_element == element:
                new_element = random.choice(building_blocks)
            mutated_sequence.append(new_element)
        else:
            mutated_sequence.append(element)
    return mutated_sequence


class GAChromosome(object):
    """
    Base chromosome class for use in the genetic algorithm.
    """

    def __init__(self,
                 sequence,
                 mutation_rate,
                 mutation_function=basic_mutation,
                 elements=[]):
        self._sequence = tuple(sequence)
        self.fitness = None
        self.mutation_rate = mutation_rate
        self._mutation_function = mutation_function
        self.elements = elements

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, new_sequence):
        self._sequence = tuple(new_sequence)

    @classmethod
    def random_chromosome(cls,
                          length,
                          mutation_rate,
                          mutation_function=basic_mutation,
                          elements=[]):
        return cls([random.choice(elements) for _ in range(length)],
                   mutation_rate,
                   mutation_function,
                   elements)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.sequence == other.sequence
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def mutate(self):
        self.sequence = self._mutation_function(self.sequence,
                                                self.mutation_rate,
                                                self.elements)
    @property
    def seq_as_string(self):
        return "-".join(self._sequence)

    @property
    def idx_name(self):
        return '-'.join([str(self.elements.index(s)) for s in self.sequence])