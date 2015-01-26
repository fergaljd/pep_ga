#!/usr/bin/env python
import random


#Recombination functions
def single_point_crossover(seq1, seq2):
    """
    The amino acid sequence
    of parent 1 (P1) and parent 2 (P2) are cleaved at one
    randomly chosen position. The four parts of P1 and
    P2 are recombined at the cleavage point in a way
    that each child contains one part of P1 and the
    other part of P2.
    """
    assert len(seq1) == len(seq2)
    crossover_pos = random.randrange(0, len(seq1))
    new_seq = seq1[:crossover_pos]
    new_seq.extend(seq2[crossover_pos:])
    return new_seq


def double_point_crossover(seq1, seq2): 
    """P1 and P2 are recombined at two randomly chosen positions."""
    return _cross_peptides(seq1, seq2, 
                           _select_n_crossover_positions(2, len(seq1)))


def distance_bisector_crossover(seq1, seq2): 
    """P1 and P2 are recombined in the middle
    of the peptide."""
    assert len(seq1) == len(seq2)
    crossover_pos = len(seq1) / 2
    return seq1[:crossover_pos] + seq2[crossover_pos:]
    

def multi_point_crossover(seq1, seq2):
    """
    P1 and P2 are recombined
    at r randomly chosen positions where r itself is also
    a random number (r < len(p1)).
    """
    assert len(seq1) == len(seq2)
    r = random.randrange(1, len(seq1))
    return _cross_peptides(seq1, seq2, 
                           _select_n_crossover_positions(r, len(seq1)))


def uniform_crossover(seq1, seq2, recombination_threshold=0.5):
    """
    Related to the multi-point
    crossover, however, each position within the
    sequence obtains a randomly assigned probability for
    recombination. If this probability exceeds a certain
    threshold, a recombination at this position occurs
    """
    assert len(seq1) == len(seq2)
    recombination_probabilities = [random.random() for _ in range(len(seq1))]
    new_sequence = [(seq2[idx] if prob > recombination_threshold else seq1[idx])
                    for (idx, prob) in enumerate(recombination_probabilities)]
    return new_sequence


def shuffle_crossover(seq1, seq2):
    """
    Similar to the double point
    crossover P1 and P2 are recombined at two randomly selected positions.
    However, before recombination the amino acids are shuffled in both
    parents. After recombination the amino acids are
    unshuffled.
    """
    assert len(seq1) == len(seq2)
    #Shuffle the list indices
    shuffled_indices = range(len(seq1))
    random.shuffle(shuffled_indices)
    shuffle = lambda seq: tuple(seq[idx] for idx in shuffled_indices)
    shuffled_seq1, shuffled_seq2 = shuffle(seq1), shuffle(seq2)
    shuffled_child = _cross_peptides(shuffled_seq1, shuffled_seq2,
                                     _select_n_crossover_positions(1,
                                                                   len(seq1)))
    return [shuffled_child[shuffled_indices.index(i)] for i in range(len(seq1))]


def _select_n_crossover_positions(n, peptide_length):
    assert n < peptide_length
    positions = []
    while len(positions) < n:
        new_position = random.randrange(0, peptide_length)
        if new_position not in positions:
            positions.append(new_position)
    return sorted(positions)


def _cross_peptides(seq1, seq2, crossover_positions):
    assert crossover_positions
    current_parent = seq1
    offspring_sequence = []
    previous_crossover = 0
    crossover = 0
    for crossover in crossover_positions:
        offspring_sequence.extend(current_parent[previous_crossover:crossover])
        previous_crossover = crossover
        current_parent = seq1 if current_parent == seq2 else seq2
    offspring_sequence.extend(current_parent[crossover:])
    return offspring_sequence
