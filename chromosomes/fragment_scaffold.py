import random

from base_chromosome import GAChromosome, basic_mutation


class ScaffoldSizeError(Exception):
    pass


class GAScaffoldChromosome(GAChromosome):
    """
    Chromosome representing a chemical scaffold with variable positions.
    The scaffold is specified as a SMILES string, with placeholders
    (default is '[X]') denoting where fragments are to be substituted
    If no scaffold is specified, fragments will be joined head to tail.
    """
    scaffold_smiles = ''
    placeholder = '[X]'
    def __init__(self,
                 sequence,
                 mutation_rate,
                 mutation_function=basic_mutation,
                 elements = []):
        super(GAScaffoldChromosome, self).__init__(sequence,
                                                   mutation_rate,
                                                   mutation_function,
                                                   elements)
        #If no scaffold is specified, assume that fragments are joined
        #head to tail, in order
        self.scaffold_smiles = scaffold_smiles
        self.placeholder = placeholder
        if not self.scaffold_smiles:
            self.scaffold_smiles = self.placeholder * len(sequence)

    @classmethod
    def random_chromosome(cls,
                          length,
                          mutation_rate,
                          mutation_function=basic_mutation,
                          elements=[],
                          scaffold_smiles='',
                          placeholder='[X]'):
        return cls([random.choice(elements) for _ in range(length)],
                   mutation_rate,
                   mutation_function,
                   elements,
                   scaffold_smiles,
                   placeholder)

    @property
    def sequence(self):
        return super(GAScaffoldChromosome, self).sequence

    @sequence.setter
    def sequence(self, new_sequence):
        if len(new_sequence) != self.scaffold_smiles.count(self.placeholder):
            msg = ("Placeholder count  %d differs from sequence length %d"
                   % (self.scaffold_smiles.count(self.placeholder),
                      len(new_sequence)))
            raise ScaffoldSizeError(msg)
        GAChromosome.sequence.fset(self, new_sequence)

    @property
    def smiles(self):
        assert (len(self.sequence)
                == self.scaffold_smiles.count(self.placeholder))
        cur_smiles = self.scaffold_smiles
        for frag in self.sequence:
            cur_smiles = cur_smiles.replace(self.placeholder, frag, 1)
        return cur_smiles

    def __getstate__(self):
        global GAScaffoldSubclass
        GAScaffoldSubclass = self.__class__
        return (self.__dict__, self.scaffold_smiles, self.placeholder)

    def __setstate__(self, state):
        self.__class__ = scaffold_factory(state[1], state[2])
        self.__dict__.update(state[0])


def scaffold_factory(scaff_smiles, pholder='[X]'):
    class GAScaffoldSubclass(GAScaffoldChromosome):
        scaffold_smiles = scaff_smiles
        placeholder = pholder
    return GAScaffoldSubclass


class GAScaffoldSubclass(GAScaffoldChromosome):
    pass


def parse_fragments(fragment_file):
    """ Get SMILES format fragments from a text file, 1 per line."""
    fragments = []
    with open(fragment_file) as raw_fragments:
        for line in raw_fragments:
            fragment = line.strip()
            if not line.startswith('#') and fragment:
                fragments.append(fragment)
    return fragments