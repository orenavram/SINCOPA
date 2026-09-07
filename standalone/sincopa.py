import os
import logging
import Bio.SeqUtils
from time import sleep
from compute_homoplasy import compute_homoplasy
from compute_sweeps_score import compute_sweeps_score
from fix_msa import fix_msa
from adjust_tree_to_msa import fix_tree
from auxiliaries import *

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def verify_fasta_format(fasta_path):
    logger.info('Validating FASTA format')
    Bio.SeqUtils.IUPACData.ambiguous_dna_letters += 'U-'
    legal_chars = set(Bio.SeqUtils.IUPACData.ambiguous_dna_letters.lower() + Bio.SeqUtils.IUPACData.ambiguous_dna_letters)
    with open(fasta_path) as f:
        line_number = 0
        try:
            line = f.readline()
            line_number += 1
            if not line.startswith('>'):
                return f'Illegal FASTA format. First line starts with "{line[0]}" instead of ">".'
            previous_line_was_header = True
            putative_end_of_file = False
            curated_content = f'>{line[1:]}'.replace("|", "_")
            for line in f:
                line_number += 1
                line = line.strip()
                if not line:
                    if not putative_end_of_file:
                        putative_end_of_file = line_number
                    continue
                if putative_end_of_file:
                    return f'Illegal FASTA format. Line {putative_end_of_file} in MSA is empty.'
                if line.startswith('>'):
                    if previous_line_was_header:
                        return f'Illegal FASTA format. MSA contains an empty record. Both lines {line_number-1} and {line_number} start with ">".'
                    else:
                        previous_line_was_header = True
                        curated_content += f'>{line[1:]}\n'.replace("|", "_")
                        continue
                else:
                    previous_line_was_header = False
                    for c in line:
                        if c not in legal_chars:
                            return f'Illegal FASTA format. Line {line_number} contains illegal DNA character "{c}".'
                    curated_content += f'{line}\n'
        except UnicodeDecodeError as e:
            logger.info(e.args)
            line_number += 1
            return f'Illegal FASTA format. Line {line_number} contains non-ASCII character(s).'
    with open(fasta_path, 'w') as f:
        f.write(curated_content)


def verify_newick_format(tree_path):
    logger.info('Validating NEWICK format')
    pass


def verify_msa_is_consistent_with_tree(msa_path, tree_path):
    logger.info('Validating MSA is consistent with the species tree')
    tree_strains = get_tree_labels(tree_path)
    logger.info(f'Phylogenetic tree contains the following strains:\n{tree_strains}')
    with open(msa_path) as f:
        logger.info(f'Checking MSA...')
        for line in f:
            if line.startswith('>'):
                strain = line.lstrip('>').rstrip('\n')
                if strain not in tree_strains:
                    msg = f'{strain} appears in the MSA but not in the tree. ' \
                          f'Please make sure the tree contains all species in the MSA.'
                    logger.error(msg)
                    return msg
                else:
                    logger.info(f'{strain} appears in tree!')


def validate_input(msa_path, tree_path, error_path):
    logger.info('Validating input...')
    error_msg = verify_fasta_format(msa_path)
    if error_msg:
        fail(error_msg, error_path)
    error_msg = verify_newick_format(tree_path)
    if error_msg:
        fail(error_msg, error_path)
    error_msg = verify_msa_is_consistent_with_tree(msa_path, tree_path)
    if error_msg:
        fail(error_msg, error_path)


def fix_input(msa_path, tree_path, output_dir, tmp_dir):
    tree_name = os.path.split(tree_path)[-1]
    adjusted_tree = f'{output_dir}/{os.path.splitext(tree_name)[0]}_fixed{os.path.splitext(tree_name)[-1]}'
    fix_tree(msa_path, tree_path, tmp_dir, adjusted_tree)
    msa_name = os.path.split(msa_path)[-1]
    fixed_msa_path = f'{output_dir}/{os.path.splitext(msa_name)[0]}_fixed{os.path.splitext(msa_name)[-1]}'
    fix_msa(msa_path, fixed_msa_path)
    return fixed_msa_path, adjusted_tree


def sincopa(msa_path, tree_path, window_size, output_dir, tmp_dir):
    homplasy_path = f'{output_dir}/homoplasy.txt'
    control_file_path = os.path.join(tmp_dir, 'control.txt')
    sweeps_scores_path = f'{output_dir}/sweeps_scores.txt'
    sweeps_plot_path = f'{output_dir}/sweeps_scores.png'
    sweeps_summary_path = f'{output_dir}/sweeps_summary.txt'
    done_path = f'{output_dir}/done.txt'

    header = 'msa_name,max_score,number_of_sequences,centrality,msa_length,window_size,index_of_max,' \
             'mean_score,median_score,max_mean_division,max_median_division,relative_location_of_peak,' \
             'apd,pi,above95,above75,above50,above25,above05'
    with open(sweeps_summary_path, 'w') as f:
        f.write(f'{header}\n')

    compute_homoplasy(msa_path, tree_path, control_file_path, homplasy_path)

    compute_sweeps_score(msa_path, homplasy_path, sweeps_scores_path, sweeps_summary_path,
                         sweeps_plot_path, window_size, stats_writing_mode='a')

    with open(done_path, 'w'):
        pass


def main(msa_path, tree_path, window_size, output_dir_path, html_path=None):
    error_path = f'{output_dir_path}/error.txt'
    try:
        os.makedirs(output_dir_path, exist_ok=True)
        tmp_dir = f'{os.path.split(msa_path)[0]}/tmp'
        os.makedirs(tmp_dir, exist_ok=True)
        validate_input(msa_path, tree_path, error_path)
        msa_path, tree_path = fix_input(msa_path, tree_path, output_dir_path, tmp_dir)
        sincopa(msa_path, tree_path, window_size, output_dir_path, tmp_dir)
        logger.info('SUCCEEDED = True')
    except Exception as e:
        logger.info(f'SUCCEEDED = False')
        logger.error(str(e))
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_msa_path',
                        help='A path to a DNA MSA file to look for sweeps.',
                        type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
    parser.add_argument('input_tree_path',
                        help='A path to a background species tree that contains (at least) all the species in the '
                             'input MSA. The tree should be reconstructed by external data and not by the MSA provided.',
                        type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
    parser.add_argument('output_dir_path',
                        help='A path to a folder in which the sweeps analysis will be written.',
                        type=lambda path: path.rstrip('/'))
    parser.add_argument('--window_size', type=int, default=50,
                        help='The size of a window to which a score will be computed')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    main(args.input_msa_path, args.input_tree_path, args.window_size, args.output_dir_path)
