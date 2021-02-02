import os
import logging
import shutil
import Bio.SeqUtils
import CONSTANTS as CONSTS  # from /bioseq/sincopa/
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
                return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. First line in MSA starts with "{line[0]}" instead of ">".'
            previous_line_was_header = True
            putative_end_of_file = False
            curated_content = f'>{line[1:]}'.replace("|", "_")
            for line in f:
                line_number += 1
                line = line.strip()
                if not line:
                    if not putative_end_of_file: # ignore trailing empty lines
                        putative_end_of_file = line_number
                    continue
                if putative_end_of_file:  # non empty line after empty line
                    return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {putative_end_of_file} in MSA is empty.'
                if line.startswith('>'):
                    if previous_line_was_header:
                        return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. MSA contains an empty record. Both lines {line_number-1} and {line_number} start with ">".'
                    else:
                        previous_line_was_header = True
                        curated_content += f'>{line[1:]}\n'.replace("|", "_")
                        continue
                else:  # not a header
                    previous_line_was_header = False
                    for c in line:
                        if c not in legal_chars:
                            return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {line_number} in MSA contains illegal DNA character "{c}".'
                    curated_content += f'{line}\n'
        except UnicodeDecodeError as e:
            logger.info(e.args)
            line_number += 1  # the line that was failed to be read
            return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {line_number} in MSA contains one (or more) non <a href="https://en.wikipedia.org/wiki/ASCII" target="_blank">ascii</a> character(s).'
    # override the old file with the curated content
    with open(fasta_path, 'w') as f:
        f.write(curated_content)


def verify_newick_format(tree_path):
    logger.info('Validating NEWICK format')
    # TODO
    pass


def verify_msa_is_consistent_with_tree(msa_path, tree_path):
    logger.info('Validating MSA is consistent with the species tree')

    tree_strains = get_tree_labels(tree_path)
    logger.info(f'Phylogenetic tree contains the following strains:\n{tree_strains}')

    with open(msa_path) as f:
        logger.info(f'Checking MSA...')
        for line in f:
            if line.startswith('>'):
                # make sure each header appears in tree
                strain = line.lstrip('>').rstrip('\n')
                logger.debug(f'Checking {strain}...')
                if strain not in tree_strains:
                    msg = f'{strain} species appears in the input MSA but not in the phylogenetic tree. ' \
                        f'Please make sure the phylogenetic tree you provide contains (at least) all ' \
                        f'the species in the provided MSA.'
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

    # write summary header
    header = 'msa_name,max_score,number_of_sequences,centrality,msa_length,window_size,index_of_max,' \
             'mean_score,median_score,max_mean_division,max_median_division,relative_location_of_peak,' \
             'apd,pi,above95,above75,above50,above25,above05'
    with open(sweeps_summary_path, 'w') as f:
        f.write(f'{header}\n')

    # step 1
    compute_homoplasy(msa_path, tree_path, control_file_path, homplasy_path)
    sleep(CONSTS.RELOAD_INTERVAL)

    # step 2
    compute_sweeps_score(msa_path, homplasy_path, sweeps_scores_path, sweeps_summary_path,
                         sweeps_plot_path, window_size, stats_writing_mode='a')

    # make sure sweeps summary file is ready
    sleep(CONSTS.RELOAD_INTERVAL)
    with open(done_path, 'w'):
        pass


def main(msa_path, tree_path, window_size, output_dir_path, html_path):

    error_path = f'{output_dir_path}/error.txt'
    try:
        if html_path:
            run_number = initialize_html(CONSTS, output_dir_path, html_path)
            final_zip_path = f'{os.path.split(output_dir_path)[0]}/{CONSTS.WEBSERVER_NAME}_{run_number}'

        os.makedirs(output_dir_path, exist_ok=True)

        tmp_dir = f'{os.path.split(msa_path)[0]}/tmp'  # same folder as the input msa
        os.makedirs(tmp_dir, exist_ok=True)

        validate_input(msa_path, tree_path, error_path)

        msa_path, tree_path = fix_input(msa_path, tree_path, output_dir_path, tmp_dir)

        sincopa(msa_path, tree_path, window_size, output_dir_path, tmp_dir)

        if html_path:
            shutil.make_archive(final_zip_path, 'zip', output_dir_path)
            finalize_html(html_path, error_path, run_number)

    except Exception as e:
        logger.info(f'SUCCEEDED = False')
        if html_path:
            error_msg = e.args[-1]
            if os.path.exists(error_path):
                with open(error_path) as f:
                    error_msg = f.read()
            edit_failure_html(CONSTS, error_msg, html_path, run_number)
            add_closing_html_tags(html_path, CONSTS, run_number)


def finalize_html(html_path, error_path, run_number):
    succeeded = not os.path.exists(error_path)
    logger.info(f'SUCCEEDED = {succeeded}')
    if succeeded:
        edit_success_html(CONSTS, html_path, run_number)
    else:
        edit_failure_html(CONSTS, error_path, html_path, run_number)
    add_closing_html_tags(html_path, CONSTS, run_number)


def edit_success_html(CONSTS, html_path, run_number):
    update_html(html_path, 'RUNNING', 'FINISHED')

    append_to_html(html_path, f'''
            <div class="container" style="{CONSTS.CONTAINER_STYLE}" align='center'>
            <h3><b>
            <a href='{CONSTS.WEBSERVER_NAME}_{run_number}.zip' target='_blank'>Download zipped full results (textual & visual)</a>
            </b></h3>
            <br>
                <img src="outputs/sweeps_scores.png" id="sweeps_scores">
            </div>
''')


def edit_failure_html(CONSTS, error_msg, html_path, run_number):
    update_html(html_path, 'RUNNING', 'FAILED')
    append_to_html(html_path,
                   f'<div class="container" style="{CONSTS.CONTAINER_STYLE}" align="justify"><h3>\n'
                   f'<font color="red">{error_msg}</font></h3><br><br>'
                   f'Please make sure your input is OK and then try to re-run your job or '
                   f'<a href="mailto:{CONSTS.ADMIN_EMAIL}?subject={CONSTS.WEBSERVER_NAME}%20Run%20Number:%20{run_number}">'
                   f'contact us'
                   f'</a> '
                   f'for further information.<br>'
                   f'</div>\n')


def add_closing_html_tags(html_path, CONSTS, run_number):
    update_html(html_path, CONSTS.PROCESSING_MSG, '')  # remove "web server is now processing your request" message
    update_html(html_path, 'progress-bar-striped active', 'progress-bar-striped')  # stop_progress_bar

    append_to_html(html_path, f'''
            <hr>
                <h4 class=footer>
                    <p align='center'>Questions and comments are welcome! Please
                        <span class="admin_link"> 
                        <a href="mailto:{CONSTS.ADMIN_EMAIL}?subject={CONSTS.WEBSERVER_NAME.upper()}%20Run%20Number%20{run_number}">contact us</a> 
                        </span>
                    </p>
                </h4>
                <div id="bottom_links" align="center"><span class="bottom_link">
                <a href="{CONSTS.WEBSERVER_URL}" target="_blank">Home</a> 
                 | 
                <a href="{CONSTS.WEBSERVER_URL}/overview.html" target="_blank">Overview</a> 
            </span>
        </div>
        <br><br><br>
    </body>
</html>''')

    # have to be last thing that is done
    sleep(2*CONSTS.RELOAD_INTERVAL)
    update_html(html_path, CONSTS.RELOAD_TAGS, f'<!--{CONSTS.RELOAD_TAGS}-->')  # stop refresh



def initialize_html(CONSTS, output_dir_path, html_path):

    path_tokens = output_dir_path.split('/')
    # e.g., "/bioseq/data/results/sincopa/12345678/outputs"
    run_number = path_tokens[path_tokens.index(CONSTS.WEBSERVER_NAME) + 1]

    update_html(html_path, 'QUEUED', 'RUNNING')
    update_html(html_path, CONSTS.PROGRESS_BAR_ANCHOR, CONSTS.PROGRESS_BAR_TAG)  # add progress bar

    return run_number


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
                                 'input MSA. The tree should be reconstructed by en external data and not by the MSA '
                                 'provided.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('output_dir_path',
                            help='A path to a folder in which the sweeps analysis will be written.',
                            type=lambda path: path.rstrip('/'))
        parser.add_argument('--window_size', type=int, default=50,
                            help='The size of a window to which a score will be computed')
        parser.add_argument('--html_path', default=None,
                            help='A path to an html file that will be updated during the run.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))

        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

        args = parser.parse_args()

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        main(args.input_msa_path, args.input_tree_path, args.window_size, args.output_dir_path, args.html_path)


