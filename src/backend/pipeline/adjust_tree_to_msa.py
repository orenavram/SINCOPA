import os
import subprocess
import re
import logging

from auxiliaries import get_tree_labels

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def remove_bootstrap_values(in_tree_path, out_tree_path):
    with open(in_tree_path) as f:
        tree_as_str = f.read()
    tree_as_str = re.sub('\)\d+:', '):', tree_as_str)
    with open(out_tree_path, 'w') as f:
        f.write(tree_as_str)


def prune_tree(msa_path, tree_to_prune_path, tmp_dir, output_tree_path,
               remove_taxa_script_path='/groups/pupko/orenavr2/src/removeTaxa'):

    logger.debug(f'Pruning tree for {msa_path}')
    msa_name = os.path.split(msa_path)[1]

    # get list of taxa in the full tree
    tree_taxa = get_tree_labels(tree_to_prune_path)

    # get list of taxa in msa
    with open(msa_path) as f:
        msa = f.read()

    msa_taxa = re.findall(r'>(\S+)\r?\n', msa)
    list_of_taxa_to_remove = [taxon for taxon in tree_taxa if taxon not in msa_taxa]

    # list of names to prune
    list_of_taxa_to_remove_path = f'{tmp_dir}/{msa_name}.txt'
    with open(list_of_taxa_to_remove_path, 'w') as f:
        f.write('\n'.join(list_of_taxa_to_remove))

    cmd_for_pruning = f'{remove_taxa_script_path} {tree_to_prune_path} {list_of_taxa_to_remove_path} {output_tree_path}'
    logger.info(f'Fetching pruning command:\n{cmd_for_pruning}')
    subprocess.run(cmd_for_pruning, shell=True)


def fix_tree(msa_path, tree_to_adjust, tmp_dir, output_tree_path):

    tree_without_bootstrap_path = tree_to_adjust+'.no_bootstrap'

    remove_bootstrap_values(tree_to_adjust, tree_without_bootstrap_path)

    prune_tree(msa_path, tree_without_bootstrap_path, tmp_dir, output_tree_path)

    logger.info(f'Tree was adjusted successfully and was saved to {output_tree_path}')

    return output_tree_path


if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('msa_path',
                            help='A path to an MSA file to which the tree should be adjusted',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('tree_to_adjust',
                            help='A path to a species tree that contains (at least) all the species in the input MSA',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('tmp_dir',
                            help='A path to a folder in which a txt file (with the same name as the msa_file) will be'
                                 'created. All the species that do not appear in the msa (and thus will be removed) '
                                 'will be written to the file that was created',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('output_tree_path',
                            help='A path to a file in which the pruned tree will be written',
                            type=lambda path: path if os.path.exists(os.path.split(path)[0]) else parser.error(
                                f'output folder {os.path.split(path)[0]} does not exist!'))
        args = parser.parse_args()

        fix_tree(args.msa_path, args.tree_to_adjust,
                 args.tmp_dir, args.output_tree_path)
