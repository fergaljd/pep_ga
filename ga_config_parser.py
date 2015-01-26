#!/usr/bin/env python
import ConfigParser
import argparse
import traceback
import logging
import pdb

from chromosomes.sequence import ALL_AMINOS
from ga import GeneticAlgorithm
from chromosomes.fragment_scaffold import parse_fragments
from options import fitness_opts, selection_opts, \
    recombination_opts, chromosome_opts

DEFAULTS = {'parallel_fitness': False}


class GAConfigError(Exception):
    pass


def ga_from_config():
    """
    Parse config files to setup Genetic Algorithm runs.
    """
    arg_parser = argparse.ArgumentParser()
    group = arg_parser.add_mutually_exclusive_group()
    group.add_argument('--run', metavar="CONFIG_FILE",
                       help="Run GA from config file.")
    group.add_argument('--restart', metavar="GA_STATE_FILE",
                       help="Restart GA run from previous"
                            " GA state (roundX.pkl).")
    arg_parser.add_argument('--logfile', default="ga-run.log",
                            help="Log file for GA run.")
    args = arg_parser.parse_args()
    if not args:
        arg_parser.error("Must specify a config file or a run file to restart."
                         " (-h for help)")
    log_format = "%(asctime)s: %(message)s"
    logging.basicConfig(filename=args.logfile, level=logging.INFO,
                        format=log_format)
    if args.restart:
        ga_instance = GeneticAlgorithm.load_run(args.restart)
        ga_instance.run()
    else:
        ga_params = _parse_config(args.run)
        ga_instance = GeneticAlgorithm(*ga_params, **ga_kw_params)
        ga_instance.run()


def _parse_config(config_file):
    config = ConfigParser.RawConfigParser(allow_no_value=True,
                                          defaults=DEFAULTS)
    config.read(config_file)
    fitness_name = config.get("Fitness", "fitness_func")
    fitness_name = '' if not fitness_name else fitness_name
    fitness_func = fitness_opts[fitness_name]
    if fitness_name:
        if fitness_name == "vina_fitness":
            parallel = config.getboolean("General", "parallel_fitness")
            config_file = config.get("Fitness", "vina_config")
            receptor = config.get("Fitness", "vina_receptor")
            keep_pdbqts = config.getboolean("Fitness", "keep_pdbqts")
            fitness_func = fitness_opts[fitness_name](receptor,
                                                      config_file,
                                                      keep_pdbqts,
                                                      parallel)
        elif fitness_name == "multi_fitness":
            parallel = config.getboolean("General", "parallel_fitness")
            raw_receptors = config.get("Fitness", "receptors")
            receptors = [r.strip() for r in raw_receptors.split(',')]
            raw_configs = config.get("Fitness", "configs")
            configs = [c.strip() for c in raw_configs.split(',')]
            equation = config.get("Fitness", "fitness_equation")
            keep_pdbqts = config.getboolean("Fitness", "keep_pdbqts")
            fitness_func = fitness_opts[fitness_name](equation,
                                                      keep_pdbqts,
                                                      parallel,
                                                      *zip(receptors, configs))
        else:
            msg = "Unknown fitness name received: %s" % fitness_name
            raise GAConfigError(msg)
    c_type = config.get("Chromosome", "type")
    elements = []
    if c_type == 'peptide':
        c_constraint = config.get("Chromosome", "constraint")
        length = config.getint("Chromosome", "peptide_length")
        chromosome = chromosome_opts[c_type][c_constraint]
        elements = ALL_AMINOS
    if c_type == 'scaffold':
        scaffold = config.get("Chromosome", "scaffold_smiles")
        placeholder = config.get("Chromosome", "placeholder")
        length = scaffold.count(placeholder)
        fragment_file = config.get("Chromosome", "fragment_file")
        fragments = parse_fragments(fragment_file)
        chromosome = chromosome_opts[c_type](scaffold, placeholder)
        elements = fragments
    ga_params = [config.getint("General", "population_size"),
                 config.getint("General", "survivors_per_round"),
                 config.getint('General', "max_generations"),
                 config.getfloat("Chromosome", "mutation_rate"),
                 length,
                 fitness_func,
                 selection_opts[config.get("Selection", "selection_method")],
                 recombination_opts[config.get("Recombination",
                                               "recombination_method")],
                 chromosome,
                 elements]
    #logfile = config.get("General", "logfile")
    return ga_params


if __name__ == '__main__':
    try:
        ga_from_config()
    except Exception:
        traceback.print_exc()
        pdb.post_mortem()