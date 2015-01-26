import random

from PepLibGen.StructGen import StructGen as sg

from base_chromosome import GAChromosome

NON_POLAR_AMINOS = "AVMLIPFW"
POSITIVE_AMINOS = "KRH"
POLAR_AMINOS = "TYQGSCN"
NEGATIVE_AMINOS = "ED"
ALL_L_AMINOS = (NON_POLAR_AMINOS + POSITIVE_AMINOS
                + NEGATIVE_AMINOS + POLAR_AMINOS)
ALL_D_AMINOS = ALL_L_AMINOS.lower().replace('g', '')
ALL_AMINOS = ALL_L_AMINOS + ALL_D_AMINOS

#Mutation functions
def amino_acid_mutation(sequence, mutation_rate, elements):
    """
    An amino acid is randomly mutated into any other amino acid.
    """
    mutated_sequence = []
    for amino_acid in sequence:
        if random.random() < mutation_rate:
            possible_aminos = elements.replace(amino_acid, '')
            mutated_sequence.append(random.choice(possible_aminos))
        else:
            mutated_sequence.append(amino_acid)
    return mutated_sequence


def chirality_mutatation(sequence, mutation_rate, elements):
    mutated_sequence = []
    for aa in sequence:
        if random.random() < mutation_rate:
            if (aa.lower() in elements) and (aa.upper() in elements):
                mutated_resi = aa.lower() if aa.islower() else aa.upper()
                mutated_sequence.append(mutated_resi)
            else:
                mutated_sequence.append(aa)
        else:
            mutated_sequence.append(aa)
    return mutated_sequence


def aa_and_chiral_mutation(sequence, mutation_rate, elements):
    mutated_sequence = amino_acid_mutation(sequence, mutation_rate, elements)
    return chirality_mutatation(mutated_sequence, mutation_rate, elements)

def amino_acid_type_mutation(sequence, mutation_rate):
    """
    An amino
    acid is substituted randomly by another amino acid
    which belongs to another amino acid class (nonpolar,
    polar, basic, or acid).
    """
    mutated_sequence = []
    for amino_acid in sequence:
        if random.random() <= mutation_rate:
            similar_aas= set(get_aa_class(amino_acid))
            different_aas = similar_aas.symmetric_difference(set(ALL_AMINOS))
            mutated_sequence.append(random.choice(list(different_aas)))
        else:
            mutated_sequence.append(amino_acid)
    return mutated_sequence


def get_aa_class(aa):
    for aa_class in [NON_POLAR_AMINOS,
                     POSITIVE_AMINOS,
                     POLAR_AMINOS,
                     NEGATIVE_AMINOS]:
        if aa in aa_class:
            return aa_class


def nucleobase_mutation(sequence, mutation_rate):
    """
    A mutation is
    introduced in the triplet of the genetic code. Stop
    codons are avoided.
    Only applicable to natural amino acids!
    """
    from Bio.Seq import Seq
    from Bio.Seq import CodonTable
    from Bio.Alphabet import generic_dna
    nucleotides = 'ATGC'
    mutate = lambda n: random.choice(nucleotides.replace(n, nucleotides))
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    dna = "".join([standard_table.back_table[aa] for aa in sequence])
    new_dna = ""
    for nucleotide in dna:
        if random.random() <= mutation_rate:
            new_dna += mutate(nucleotide)
        else:
            new_dna += nucleotide
    new_sequence = Seq(new_dna, generic_dna).translate()
    if new_sequence.alphabet.stop_symbol in new_sequence:
        #Contains stop codon, try again
        nucleobase_mutation(sequence, mutation_rate)
    else:
        return tuple(new_sequence)


class GAPeptideChromosome(GAChromosome):
    """
    Represents an arbitrary linear peptide sequence chromosome.
    """
    constraint_type = ''
    def __init__(self,
                 sequence,
                 mutation_rate,
                 mutation_function=aa_and_chiral_mutation,
                 elements=ALL_AMINOS):
        super(GAPeptideChromosome, self).__init__(sequence,
                                                  mutation_rate,
                                                  mutation_function,
                                                  elements)

    @classmethod
    def random_chromosome(cls,
                          length,
                          mutation_rate,
                          mutation_function=amino_acid_mutation,
                          elements=ALL_AMINOS):
        cls_name = GAPeptideChromosome
        return super(cls_name, cls).random_chromosome(length,
                                                      mutation_rate,
                                                      mutation_function,
                                                      elements)

    @property
    def smiles(self):
        """The peptide SMILES string."""
        assert self.constraint_type in ['HT', 'SS', '']
        if self.constraint_type == 'SS':
            sequence, constraint_type = self._get_disulphide_seq_const()
        else:
            sequence, constraint_type = self.sequence, self.constraint_type
        _, _, smiles = sg.constrained_peptide_smiles(sequence, constraint_type)
        return smiles

    def _get_disulphide_seq_const(self):
        constraint_type = 'SSC' + 'X'*len(self.sequence) + 'C'
        sequence = list(self.sequence)
        sequence.insert(0, 'C')
        sequence.append('C')
        return sequence, constraint_type

    @property
    def idx_name(self):
        return ''.join(self.sequence)


class GAHeadTailPeptideChromosome(GAPeptideChromosome):
    """
    Represents an arbitrary head-tail peptide sequence chromosome.
    """
    constraint_type = 'HT'


class GADisulphidePeptideChromosome(GAPeptideChromosome):
    """
    Represents an arbitrary disulphide wrapped peptide sequence chromosome.
    """
    constraint_type = 'SS'
