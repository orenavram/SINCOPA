# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 18:53:44 2016

@author: Oren
"""

import sys
sys.path.append('/bioseq/bioSequence_scripts_and_constants/')  # ADD file_writer
import os
import logging
import SIMULATION_CONSTANTS as SIM_CONSTS
import subprocess
from sklearn.model_selection import ParameterGrid
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def simulate_loci(parameters_dict, num_of_simulations, window_size, output_path,
                  queue, simulations_per_job, debug, new_line_delimiter='!@#!@#'):

    os.makedirs(output_path, exist_ok=True)

    parametes_grid = ParameterGrid(parameters_dict)
    logger.info(f'Total number of configurations is: {len(parametes_grid)}')
    logger.info(f'Each configuration will be repeated {len(parametes_grid)} times')

    for configuration_number, configuration in enumerate(parametes_grid):

        logger.info(f'Starting to fetch simulation pipeline of configuration {configuration_number}:\n{configuration}\n')
        configuration_params_as_str = '_'.join(str(value) for value in configuration.values())

        configuration_path = f'{output_path}/{configuration_params_as_str}'
        os.makedirs(configuration_path, exist_ok=True)

        configuration_tmp_path = f'{output_path}/{configuration_params_as_str}/tmp'
        os.makedirs(configuration_tmp_path, exist_ok=True)

        cmds = []
        for simulation_number in range(num_of_simulations):

            current_command = f'module load {" ".join(SIM_CONSTS.REQUIRED_MODULES)}; '
            logger.debug(f'Starting to fetch simulation pipeline for simulation number {simulation_number} of configuration {configuration_number}:\n{configuration}\n')
            simulation_path = f'{configuration_path}/simulations/{simulation_number}'
            os.makedirs(simulation_path, exist_ok=True)

            # output files paths
            simulated_loci_path = f'{simulation_path}/all_loci.txt'
            simulated_non_recombining_loci_path = f'{simulation_path}/non_recombining_loci.txt'
            tree_path = f'{simulation_path}/tree.txt'
            sincopa_output = f'{simulation_path}/sweeps_summary.txt'

            # sequences simulation
            if not (os.path.exists(simulated_loci_path) and os.path.exists(simulated_non_recombining_loci_path)):
                control_file_path = create_control_file(configuration, simulation_path)
                current_command += f'{SIM_CONSTS.SIMULATE_LOCI_SCRIPT_PATH} {control_file_path}; '
            else:
                logger.info(f'Skipping simulation. Output files already exist at:\n{simulated_loci_path}')

            # tree reconstruction
            if not os.path.exists(tree_path):
                current_command += f'python {SIM_CONSTS.RECONSTRUCT_TREE_SCRIPT_PATH} {simulated_non_recombining_loci_path} {tree_path}; '
            else:
                logger.info(f'Skipping tree reconstruction. Output file already exists at:\n{tree_path}')

            # sweeps analysis
            if not os.path.exists(sincopa_output):
                current_command += f'python {SIM_CONSTS.SINCOPA_SCRIPT_PATH} {simulated_loci_path} {tree_path} ' + \
                                   f'{simulation_path}/sweeps_analysis --window_size {window_size}; '
            else:
                logger.info(f'Skipping sweeps analysis. Output file already exists at:\n{sincopa_output}')

            # sequence of commands from simulation to sweeps analysis
            cmds.append(current_command)

        # a command for aggregating results
        cmds.append(f'python {SIM_CONSTS.SWEEPS_RESULTS_AGGREGATION_SCRIPT_PATH} {configuration_path} {num_of_simulations}; ')

        for batch_number in range(0, len(cmds), simulations_per_job):
            current_batch = cmds[batch_number:batch_number+simulations_per_job]
            cmds_file = os.path.join(configuration_tmp_path, f'qsub_{batch_number}.cmds')

            job_name = f'LociSimConf{configuration_number}_'
            if batch_number + simulations_per_job < len(cmds):
                job_name += str(batch_number)
            else:
                # last batch
                job_name += 'Last'

            with open(cmds_file, 'w') as f:
                f.write(f'{new_line_delimiter.join(current_batch)}\t{job_name}')

            job_id_file = os.path.join(configuration_tmp_path, f'job_id_{batch_number}.txt')

            # a simple command when using shebang header (#!) in q_submitter_power.py
            submission_cmd = f'/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py {cmds_file} {configuration_tmp_path} -q "{queue}" > {job_id_file}'
            subprocess.call(submission_cmd, shell=True)


def create_control_file(configuration, simulation_path):
    ctrl_file_txt = f'{configuration["population_size"]} #pop_size\n' \
        f'{configuration["num_of_populations"]} #num_of_pops\n' \
        f'{configuration["num_of_generations"]} #num_of_generations\n' \
        f'{configuration["super_generation"]} #super_generation (one population takes over another one)\n' \
        f'{configuration["num_of_loci"]} #num_of_loci\n' \
        f'{configuration["locus_len"]} #locus_len\n' \
        f'{configuration["num_of_non_recombining_loci"]} #num_of_non_recombining_loci\n' \
        f'{configuration["selective_allele_timing_proportion"]} #selective_allele_timing_proportion\n' \
        f'{configuration["selection_intensity"]} #selection_intensity. set to 1 for no selection\n' \
        f'{configuration["p_a"]} #p_a\n' \
        f'{configuration["p_c"]} #p_c\n' \
        f'{configuration["p_g"]} #p_g\n' \
        f'{configuration["p_transition"]} #p_transition\n' \
        f'{configuration["mutation_rate_per_nucleotide_per_step"]} #mutation_rate_per_nucleotide_per_step\n' \
        f'{configuration["recombination_to_mutaiton_rate"]} #r_to_m_rate\n' \
        f'{configuration["between_pops_recombination_rate_to_recombination_rate_per_bacterium_per_step"]} #between_pops_recombination_rate_to_recombination_rate_per_bacterium_per_step\n' \
        f'{simulation_path}/ #output_directory\n'

    control_file_path = f'{simulation_path}/control_file.txt'
    with open(control_file_path, 'w') as f:
        f.write(ctrl_file_txt)

    return control_file_path


if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('output_path', type=lambda path: path.rstrip('/'),
                        help='A path to a directory in which the simulation results will be written to.')
    parser.add_argument('-n', '--num_of_simulations', default=100, type=int,
                        help='number of sweeps score observation for each configuration')
    parser.add_argument('--window_size', type=int, default=50,
                        help='The size of a window to which a sweeps score will be computed')
    parser.add_argument('-q', '--queue', default='pupkolabr', help='The queue to which the jobs will be submitted to')
    parser.add_argument('-j', '--simulations_per_job', default=1, type=int, help='How many simulations will be fetched in each job')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('-d', '--debug', action='store_true', help='Use the debug dictionary')

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    parameters_dict = {'population_size': [10_000],  # , '10000'],  # 10,000 is too much!!
                       'num_of_populations': [100],  # Num of samples at the end
                       'num_of_generations': [10_000],  # '1000000' is too big
                       'super_generation': [str(2 ** 31 - 1)], #, str(2 ** 31 - 1)],
                       'locus_len': [250],
                       # mean E. coli recombination length. took from: http://www.pnas.org/content/suppl/2014/10/09/1413272111.DCSupplemental/pnas.201413272SI.pdf
                       'num_of_loci': [5],  # 2 non-recombining, 1 selective, 5 non-selective
                       'num_of_non_recombining_loci': [2],
                       'selective_allele_timing_proportion': [0.99],
                       'p_a': [0.25],  # Observed frequency
                       'p_c': [0.25],  # Observed frequency
                       'p_g': [0.25],  # Observed frequency
                       'p_transition': [0.5],
                       'selection_intensity': [2], # '10', '50', '100', '1000'],
                       # there's no 'inf' in c++ so I use the largest int (2**31-1) to shut the superstep down
                       'mutation_rate_per_nucleotide_per_step': [1e-6],  # '0.0000001'],
                       'recombination_to_mutaiton_rate': [0.38],
                       # hundredth with respect to within_pop_recombination
                       'between_pops_recombination_rate_to_recombination_rate_per_bacterium_per_step': [0.05],
                       }

    if args.debug:
        parameters_dict = {'population_size': ['100'],  # , '10000'],  # 10,000 is too much!!
                           'num_of_populations': ['10'],  # Num of samples at the end
                           'num_of_generations': ['100'],  # '1000000' is too big
                           'super_generation': [str(2 ** 31 - 1)],  # '100000', str(2 ** 31 - 1)],
                           'locus_len': ['250'],
                           # mean E. coli recombination length. took from: http://www.pnas.org/content/suppl/2014/10/09/1413272111.DCSupplemental/pnas.201413272SI.pdf
                           'num_of_loci': ['5'],  # 2 non-recombining, 1 selective, 3 non-selective
                           'num_of_non_recombining_loci': ['2'],
                           'selective_allele_timing_proportion': [0.8],
                           'p_a': ['0.2421'],  # Observed frequency
                           'p_c': ['0.2444'],  # Observed frequency
                           'p_g': ['0.2691'],  # Observed frequency
                           'p_transition': ['0.5'],
                           'selection_intensity': ['5'],# '10', '50', '100', '1000'],
                           # there's no 'inf' in c++ so I use the largest int (2**31-1) to shut the superstep down
                           'mutation_rate_per_nucleotide_per_step': ['0.00001'],  # '0.0000001'],
                           'recombination_to_mutaiton_rate': ['0.38'],
                           # hundredth with respect to within_pop_recombination
                           'between_pops_recombination_rate_to_recombination_rate_per_bacterium_per_step': ['0.01'],
                           }

    simulate_loci(parameters_dict, args.num_of_simulations, args.window_size,
                  args.output_path, args.queue, args.simulations_per_job, args.debug)


