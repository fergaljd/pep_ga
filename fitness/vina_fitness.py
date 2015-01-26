#!/usr/bin/env python
import os
import re
from os import path
from subprocess import check_call
import logging
import tempfile
from contextlib import closing
import multiprocessing

MIN_SCORE = 0.01
#Fitness functions
#All fitnesses need to be > 0 and positive for selection functions to work!


class VinaFitness(object):
    """Calculates the Vina fitness.
    This class allows storing all of the config and receptor information
    on setting up the instance, to allow the fitness to be calculated by
    passing only the list of chromosomes.

    This class can also be pickled, in constrast with alternative approach
    of returning nested functions.
    """
    def __init__(self, receptor_f, config_f, keep_outfile, parallel=False):
        self.receptor_f = receptor_f
        self.config_f = config_f
        self.keep_outfile = keep_outfile
        self.parallel = parallel
        assert path.exists(receptor_f)
        assert path.exists(config_f)

    def __call__(self, chromosomes):
        fitnesses = []
        if not self.parallel:
            for chromosome in chromosomes:
                f = vina_fitness(chromosome, self.receptor_f,
                                 self.config_f, self.keep_outfile, cpus=1)
            fitnesses.append(f)
        else:
            pool = multiprocessing.Pool(multiprocessing.cpu_count())
            arg_supplier = [(chromosome, self.receptor_f,
                             self.config_f, self.keep_outfile, 1)
                            for chromosome in chromosomes]
            fitnesses = list(pool.map(_starred_vina_fitness, arg_supplier))
        return fitnesses



def _starred_vina_fitness(args):
    return vina_fitness(*args)


def vina_fitness(chromosome, receptor_f, config_f, keep_outfile=True, cpus=0):
    """Dock a chromosome representation to a target to assess fitness.
    Chromosome must have a valid SMILES string as the SMILES property."""
    ligand_f = pdbqt_from_smiles(chromosome.smiles, chromosome.idx_name)
    try:
        fitness = vina_abs_docking_score(ligand_f, receptor_f,
                                         config_f, keep_outfile, cpus)
        logging.info("Calculated chromosome fitness: Sequence %s, Fitness %s,"
                     " Receptor %s"
                     % (chromosome.seq_as_string, fitness,
                        path.basename(receptor_f)))
    finally:
        _silent_delete(ligand_f)
    return fitness

def pdbqt_from_smiles(smiles, mol_name):
    """Write an Autodock4 PDBQT file from a SMILES string.
    Returns the filename."""
    cmdline_babel = False
    try:
        import pybel
    except ImportError:
        cmdline_babel = True
    mol_pdbname = mol_name + '.pdb'
    if not cmdline_babel:
        mol = pybel.readstring("smi", smiles)
        mol.make3D()
        mol.write(format="pdb", filename=mol_pdbname, overwrite=True)
    else:
        in_smi = '-:%s' % smiles
        check_call(['obabel', in_smi, '-O', mol_pdbname, '--gen3d'],
                              stdout=open(os.devnull), stderr=open(os.devnull))
    mol_pdbqtname = mol_pdbname + 'qt'
    if not path.exists(mol_pdbqtname):
        try:
            check_call(['prepare_ligand4.py',
                        '-l', mol_pdbname, '-o', mol_pdbqtname],
                        stdout=open(os.devnull), stderr=open(os.devnull))
        finally:
            os.remove(mol_pdbname)
    return mol_pdbqtname


def _silent_delete(f):
    try:
        os.remove(f)
    except OSError:
        logging.exception("Error deleting file %s" % f)
        pass


def _do_docking(ligand_f, receptor_f, config_f, cpus=0):
    """ Runs a Vina docking run.
    Returns the output PDBQT file name."""
    raw_fname = lambda f: path.splitext(path.basename(f))[0]
    outfile_name = "%s_out_%s.pdbqt" % (raw_fname(ligand_f),
                                        raw_fname(receptor_f))
    args = ['vina',
            '--receptor', receptor_f,
            '--ligand', ligand_f,
            '--config', config_f,
            '--out', outfile_name]
    if cpus:
        args.extend(['--cpu', str(cpus)])
    with closing(tempfile.TemporaryFile()) as output:
        check_call(args, stdout=output, stderr=open(os.devnull))
        output.seek(0)
        random_seed = _parse_vina_stdout(output.read())
    return outfile_name, random_seed


def _parse_vina_stdout(vina_output):
    pattern = re.compile(r"random seed: (-?\d+)")
    match = pattern.search(vina_output)
    random_seed = int(match.group(1))
    return random_seed

def _parse_vina_outfile(outfile):
    """Returns all docking scores from a Vina output pdbqt file. """
    with open(outfile) as raw_f:
        scores = []
        splitter = re.compile(' +')
        for line in raw_f:
            if 'REMARK VINA RESULT' in line:
                raw_score = splitter.split(line)[3]
                scores.append(float(raw_score))
    return scores


def vina_abs_docking_score(ligand_f, receptor_f, config_f,
                            cpus=0, keep_outfile=False):
    vina_outfile, random_seed = _do_docking(ligand_f, receptor_f,
                                            config_f, cpus)
    try:
        docking_score = (_parse_vina_outfile(vina_outfile)[0])
        logging.info("Docked ligand %s with random seed %d and"
                     " raw docking score %f"
                     % (ligand_f, random_seed, docking_score))
        #Make higher scores better, and all scores positive
        docking_score = abs(docking_score) if docking_score < 0 else MIN_SCORE
        return docking_score
    finally:
        if not keep_outfile:
            _silent_delete(vina_outfile)

