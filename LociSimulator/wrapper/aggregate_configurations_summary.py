import pandas as pd
import matplotlib.pyplot as plt
import os
import logging
from subprocess import check_output
from time import sleep
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def get_number_of_results(configuration_simulations_path):
    path_to_look_for = f'{configuration_simulations_path}/simulations/*/sweeps_analysis/done.txt'
    try:
        return int(check_output(f'ls -1 {path_to_look_for} | wc', shell=True).decode().split()[0])
    except:
        logger.info('ls command failed.. returning 0.')
        return 0


def aggregate_configuration_results(configuration_simulations_path, number_of_expected_results, sleeping_time=10):

    logger.info(f'\n\nWaiting for {number_of_expected_results} sweeps summary results (i.e., done.txt files)...')
    current_number_of_results = -1
    while current_number_of_results < number_of_expected_results:
        current_number_of_results = get_number_of_results(configuration_simulations_path)
        logger.info(f'{current_number_of_results} results are ready.')
        sleep(sleeping_time)

    logger.info(f'\n\n')

    # aggregate results
    configuration_summary = ''
    for path, folders, files in os.walk(configuration_simulations_path):
        for file in files:
            if file.endswith('sweeps_summary.txt'):
                sweeps_summary_file_path = os.path.join(path, file)
                logger.info(f'Adding {sweeps_summary_file_path}')
                with open(sweeps_summary_file_path) as f:
                    # the sweeps summary file contains only 2 rows: a header row and a stats row
                    if not configuration_summary:
                        # first time stats that are added
                        configuration_summary += f'{f.readline().rstrip()}\n'  # add header
                    else:
                        f.readline()  # skip header
                    configuration_summary += f'{f.readline().rstrip()}\n'

    # save aggregated results
    configuration_summary_text_path = f'{configuration_simulations_path}/configuration_summary.txt'
    with open(configuration_summary_text_path, 'w') as f:
        f.write(configuration_summary)
    logger.info(f'Configuration summary table is ready at:\n{configuration_summary_text_path}')

    df = pd.read_csv(configuration_summary_text_path)
    for column in df.columns:
        if column == 'msa_name' or 'above' in column:
            logger.info(f'Skipping {column} column...')
            continue

        fig = plt.figure()
        plt.title(f'{column}\n')
        plt.hist(df[column], bins=50)  # TODO: density=True
        fig.savefig(f'{configuration_simulations_path}/{column}.png', dpi=300, bbox_inches='tight')
        plt.close()
        fig = plt.figure()
        plt.title(f'{column}\n')
        df[column].plot.kde()
        fig.savefig(f'{configuration_simulations_path}/{column}_smooth.png', dpi=300, bbox_inches='tight')
        plt.close()

    logger.info(f'Plots are ready as well!')


if __name__ == '__main__':
        from sys import argv
        logger.info(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('configuration_simulations_path', help='A path to simulations results of a single configuration.',
                            type=lambda path: path.rstrip('/') if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('number_of_results', default=0, type=int,
                            help='waits until have enough results (in case of parallelizing each configuration)')

        args = parser.parse_args()

        aggregate_configuration_results(args.configuration_simulations_path, args.number_of_results)
