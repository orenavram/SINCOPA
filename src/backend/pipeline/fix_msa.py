import os
import logging

from auxiliaries import load_header2sequences_dict

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def fix_ambiguous_columns(column, strain2sequence, legal_chars='ACGT-'):
    col = [strain2sequence[strain][column] for strain in strain2sequence
           if strain2sequence[strain][column] in legal_chars]

    major_allele = max(set(col), key=col.count)
    for strain in strain2sequence:
        if strain2sequence[strain][column] not in legal_chars:
            logger.info(f'Ambiguous character {strain2sequence[strain][column]} was detected. Replacing column #{column} of {strain} from {strain2sequence[strain][column]} to {major_allele}')
            strain2sequence[strain] = strain2sequence[strain][:column] + major_allele + strain2sequence[strain][column+1:]


def get_first_column_without_gap(strain2sequence, indexes):
    for col_index in indexes:
        for strain in strain2sequence:
            if strain2sequence[strain][col_index] == '-':
                break
        else:
            # no gap was found!
            return col_index


def trim_msa(strain2sequence, msa_length):
    first_left_col_without_gaps = get_first_column_without_gap(strain2sequence, range(msa_length))  # left side
    first_right_col_without_gaps = get_first_column_without_gap(strain2sequence, range(msa_length-1, -1, -1))  # right side

    for strain in strain2sequence:
        strain2sequence[strain] = strain2sequence[strain][first_left_col_without_gaps: 1+first_right_col_without_gaps]

    return 1+first_right_col_without_gaps - first_left_col_without_gaps


def fix_msa(msa_path, output_path, minimal_length=300):
    strain2sequence, msa_length = load_header2sequences_dict(msa_path, get_length=True, upper_sequence=True)

    for strain in strain2sequence:
        if len(strain2sequence[strain]) != msa_length:
            raise ValueError(f'Illegal MSA. Not all sequences are of the same length. E.g., {strain} sequence length is {len(strain2sequence[strain])} where others are of length {msa_length}.')

    for col in range(msa_length):
        fix_ambiguous_columns(col, strain2sequence)

    logger.info(f'Trimming {os.path.split(msa_path)[-1]}')
    logger.info(f'MSA length before trimming is {msa_length}bps')
    trimmed_msa_length = trim_msa(strain2sequence, msa_length)
    logger.info(f'MSA length after trimming is {trimmed_msa_length}bps (a total of {msa_length - trimmed_msa_length}bps were trimmed)')

    if trimmed_msa_length < minimal_length:
        logger.info(f'Discarding too short MSA (needs to be at least {minimal_length}bps wide).')
        return

    fixed_alignment = ''
    for strain in strain2sequence:
        fixed_alignment += '>' + strain + '\n' + strain2sequence[strain] + '\n'

    with open(output_path, 'w') as f:
        f.write(fixed_alignment)

    logger.info(f'MSA was fixed and stored successfully at {output_path}')


if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('input_msa_path',
                            help='A path to an MSA file that should be fixed. Fixation includes: '
                                 '(1) ambiguous chars removal and (2) msa edges trimming.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('output_msa_path',
                            help='A path to the fixed MSA',
                            type=lambda path: path if os.path.exists(os.path.split(path)[0]) else parser.error(
                                f'output folder {os.path.split(path)[0]} does not exist!'))
        parser.add_argument('--drop-shorter', default=300, help='A shorter MSA will be discarded.',
                            type=lambda x: int(x) if int(x) > 0 else parser.error(f'Minimal number of upstream basepairs should be positive!'))
        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

        args = parser.parse_args()

        fix_msa(args.input_msa_path, args.output_msa_path, args.drop_shorter)


