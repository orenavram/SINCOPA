import logging
import sys

from auxiliaries import load_header2sequences_dict

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')

import re
import os
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from sklearn import metrics
from numpy import mean, median, argmax
from collections import Counter
from itertools import combinations


def get_sequences_from_fasta(msa_path):
    # return a list of the sequences without inner "\n"s
    sequences = re.split('@', re.sub('\r?\n', '', re.sub('>.+\r?\n', '@', open(msa_path).read())))[1:]
    return sequences


def get_pairwise_distance(seq1, seq2):
    return float(sum([c1 != c2 for c1, c2 in zip(seq1, seq2)]))


def get_average_pairwise_distance_and_pi(header2sequence, msa_length):
    number_of_species = len(header2sequence)
    print(number_of_species)
    num_of_pairs = number_of_species * (number_of_species-1) / 2  # number_of_species choose 2
    total_relative_pairwise_distance = 0
    pi = 0

    sequence2count = {}
    for header in header2sequence:
        sequence2count[header] = sequence2count.get(header2sequence[header], 0) + 1

    sequences2frequency = {sequence: sequence2count[sequence]/number_of_species for sequence in sequence2count}

    sequence_pairs = set(combinations(sequence2count, 2))
    for sequence_pair in sequence_pairs:
        apd = get_pairwise_distance(*sequence_pair)/msa_length
        total_relative_pairwise_distance += apd
        pi += 2 * sequences2frequency[sequence_pair[0]] * sequences2frequency[sequence_pair[1]] * apd

    logger.info(f'pi is {pi}')
    logger.info(f'total relative apd is {total_relative_pairwise_distance / num_of_pairs}')
    return total_relative_pairwise_distance / num_of_pairs, pi


def get_alignment_columns_from_dict_msa(header2sequence, msa_length):
    msa_columns = []

    # iterate over msa columns
    for j in range(msa_length):
        # extract coulmn j
        column = ''.join([header2sequence[header][j] for header in header2sequence])
        msa_columns.append(column)

    return msa_columns, len(column)


def load_homoplasy(homoplasy_path):
    # return a list of ints indicating whether the corresponding columns in an alignment are homoplasy or not
    with open(homoplasy_path) as f:
        homoplasy = [int(x) for x in f.read().split()]

    return homoplasy


# calculate the Hamming similarity of two given columns
def calculate_hamming_similarity(col1, col2):

    chars1 = sorted(set(col1), key=col1.count, reverse=True)
    chars2 = sorted(set(col2), key=col2.count, reverse=True)

    for i in range(len(chars1)):
        col1 = col1.replace(chars1[i], str(i))

    for i in range(len(chars2)):
        col2 = col2.replace(chars2[i], str(i))

    return sum(c1 == c2 for c1, c2 in zip(col1, col2)) / len(col1)


# calculate the Eli similarity of two given columns
def calculate_greedy_association(col1, col2):
    res = 0
    char_to_num = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': 4}
    pairs_occurences = [[0] * 5, [0] * 5, [0] * 5, [0] * 5, [0] * 5]
    for i in range(len(col1)):
        pairs_occurences[char_to_num[col1[i]]][char_to_num[col2[i]]] += 1

    counts = Counter(col1)
    sortedCounts = [x[0] for x in sorted(counts.items(), key=lambda x: x[1], reverse=True)]
    for char in sortedCounts:
        maxOccurenceIndex = argmax(pairs_occurences[char_to_num[char]])
        res += pairs_occurences[char_to_num[char]][maxOccurenceIndex]

        for i in range(len(char_to_num.keys())):
            pairs_occurences[i][maxOccurenceIndex] = -1

    return res / len(col1)


# def calculate_greedy_association2(col1, col2, chars='ACGT-'):
#     res = 0
#     char2char2charcount = {char: dict.fromkeys(chars, 0) for char in chars}
#
#     print(char2char2charcount)
#     for ch1, ch2 in zip(col1, col2):
#         char2char2charcount[ch1][ch2] = char2char2charcount[ch1].get(ch2, 0) + 1
#
#     counts = Counter(col1)
#     sorted_chars_by_frequency = [char for char in sorted(counts, key=counts.get, reverse=True)]
#     for char in sorted_chars_by_frequency:
#         max_occurence_mate_pair = max(char2char2charcount[char], key=char2char2charcount[char].get)
#         res += char2char2charcount[char][max_occurence_mate_pair]
#
#         char2char2charcount[char] = {char: dict.fromkeys(chars, -1) for char in chars}
#
#     return res / len(col1)


def calculate_symmetric_greedy_association(col1, col2):
    # calculate a symmetric greedy association of two given columns
    # assert calculate_greedy_association(col1, col2) == calculate_greedy_association2(col1, col2)
    # assert calculate_greedy_association(col2, col1) == calculate_greedy_association2(col2, col1)
    return (calculate_greedy_association(col1, col2) + calculate_greedy_association(col2, col1)) / 2


def get_ith_window(iterable, i, window_size):
    return iterable[i: i + window_size]


def calculate_total_edge_contribution(edge_column, similarity_function,
                                      window_columns_without_edge_column, window_homoplasy_without_edge_column):

    # calculate the contribution of an edge column to the total window score
    edge_contributions = 0

    for column, homoplasy in zip(window_columns_without_edge_column, window_homoplasy_without_edge_column):
        if homoplasy:
            # otherwise, ignore this column
            edge_contributions += similarity_function(edge_column, column)

    return edge_contributions


def get_first_window_score(window_columns, window_homoplasy, similarity_function):
    # calculate a window score.
    # inefficient for second window (etc...) since the every window overlaps with the previous one
    window_score = 0

    homoplasious_window_columns = [column for i, column in enumerate(window_columns) if window_homoplasy[i]]
    for i in range(len(homoplasious_window_columns) - 1):
        for j in range(i + 1, len(homoplasious_window_columns)):
            window_score += similarity_function(homoplasious_window_columns[i], homoplasious_window_columns[j])

    return window_score


def get_next_window_score(prev_window_score, similarity_function,
                          prev_window_columns, prev_window_homoplasy,
                          window_columns, window_homoplasy):

    leftmost_edge_contribution = rightmost_edge_contribution = 0

    if prev_window_homoplasy[0]:
        # previous leftmost edge has homoplasy; otherwise, no contribution...
        edge_column = prev_window_columns[0]
        leftmost_edge_contribution = calculate_total_edge_contribution(edge_column, similarity_function,
                                                                       prev_window_columns[1:], prev_window_homoplasy[1:])

    if window_homoplasy[-1]:
        # current rightmost edge has homoplasy; otherwise, no contribution...
        edge_column = window_columns[-1]
        rightmost_edge_contribution = calculate_total_edge_contribution(edge_column, similarity_function,
                                                                        window_columns[:-1], window_homoplasy[:-1])

    return prev_window_score - leftmost_edge_contribution + rightmost_edge_contribution


def get_contig_scores(msa_columns, msa_length, homoplasy, window_size, similarity_function, normalization_factor):

    # setting similarity function
    if similarity_function == 'symmetric_greedy':
        similarity_function = calculate_symmetric_greedy_association
    elif similarity_function == 'greedy':
        similarity_function = calculate_greedy_association
    elif similarity_function == 'mutual':
        similarity_function = metrics.mutual_info_score
    else:
        similarity_function = calculate_hamming_similarity

    # calculates scores for every window in msa
    number_of_windows = msa_length - window_size + 1
    logger.info(f'Total number of windows is {number_of_windows}...')

    # scores vector initialization
    window_scores = [0] * (msa_length - window_size + 1)

    # first window initialization:
    logger.info(f'Calculating score for window #0 with {similarity_function.__name__} as similarity function')
    window_columns = get_ith_window(msa_columns, 0, window_size)
    window_homoplasy = get_ith_window(homoplasy, 0, window_size)
    window_scores[0] = get_first_window_score(window_columns, window_homoplasy, similarity_function)

    # sliding window over the next windows
    for i in range(1, number_of_windows):
        if i % 100 == 0:
            logger.info(f'Calculating score for window #{i}...')
        prev_window_columns = window_columns
        prev_window_homoplasy = window_homoplasy
        window_columns = get_ith_window(msa_columns, i, window_size)
        window_homoplasy = get_ith_window(homoplasy, i, window_size)

        window_scores[i] = get_next_window_score(window_scores[i - 1], similarity_function,
                                                  prev_window_columns, prev_window_homoplasy,
                                                  window_columns, window_homoplasy)

    # normalizing each score by the answer of possible pairs in a window
    # normalization must be done AFTER the all the scores were obtained since in order to compute each score
    # efficiently, current (unnormalized) score uses previous(unnormalized) score!
    normalized_window_scores = [score * normalization_factor for score in window_scores]
    return normalized_window_scores


def write_summary(msa_name, scores, window_size, header2sequence, msa_length, number_of_sequences,
                  meta_output_path, mode='w'):

    max_score = max(scores)
    index_of_max = scores.index(max_score) + 1  # start from 1 rather than 0
    mean_score = mean(scores)
    median_score = median(scores)
    max_mean_division = -1 if mean_score == 0 else max_score / mean_score
    max_median_division = -1 if median_score == 0 else max_score / median_score
    relative_location_of_peak = index_of_max / msa_length
    centrality = min(1 - relative_location_of_peak, relative_location_of_peak) * 2  # scaling from 0 to 1
    apd, pi = get_average_pairwise_distance_and_pi(header2sequence, msa_length)

    above05 = above25 = above50 = above75 = above95 = 0
    if max_score > 0.95:
        above05 = above25 = above50 = above75 = above95 = 1
    elif max_score > 0.75:
        above05 = above25 = above50 = above75 = 1
    elif max_score > 0.5:
        above05 = above25 = above50 = 1
    elif max_score > 0.25:
        above05 = above25 = 1
    elif max_score > 0.05:
        above05 = 1

    meta_data = [msa_name,max_score,number_of_sequences,centrality,msa_length,window_size,index_of_max,
                 mean_score,median_score,max_mean_division,max_median_division,relative_location_of_peak,
                 apd,pi,above95,above75,above50,above25,above05]

    with open(meta_output_path, mode) as f:
        f.write(','.join([meta_data[0]] + [f'{abs(score):.4f}' for score in meta_data[1:]]) + '\n')


def plot_scores(scores, msa_length, window_size, output_path, similarity_function,
                title='', x_label='Window #', y_label='SINCOPA Score', epsilon=0.05):

    # create plot
    plt.plot(range(msa_length - window_size + 1), scores)

    # set axes limits
    # plt.xlim(-1, msa_length - window_size + 1)
    plt.ylim(0 - epsilon, 1 + epsilon)

    # set labels
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if similarity_function == 'mutual':
        # Semi-log for mutual info. Irrelelvant for symmetric_greedy.
        plt.yscale('log')

    # save
    plt.savefig(output_path)
    plt.close()


def compute_sweeps_score(msa_path, homoplasy_path, scores_output_path, meta_output_path, plot_path, window_size,
                         similarity_function='symmetric_greedy', stats_writing_mode='w'):

    logger.info('Starting to compute_sweeps_score...')

    header2sequence, msa_length = load_header2sequences_dict(msa_path, get_length=True, upper_sequence=True)

    msa_columns, number_of_sequences = get_alignment_columns_from_dict_msa(header2sequence, msa_length)
    msa_name = os.path.split(msa_path)[-1]
    homoplasy = load_homoplasy(homoplasy_path)

    assert msa_length == len(homoplasy), 'MSA length is inconsistent with homoplasy length ' \
                                         '(they should be both of the same length).'

    scores = get_contig_scores(msa_columns, msa_length, homoplasy, window_size, similarity_function,
                               normalization_factor=2/(window_size * (window_size-1)))  # 1 / window_size choose 2

    # write scores to file
    with open(scores_output_path, 'w') as f:
        f.write('\n'.join([f'{abs(score):.4f}' for score in scores]) + '\n')

    write_summary(msa_name, scores, window_size, header2sequence, msa_length, number_of_sequences,
                  meta_output_path, stats_writing_mode)

    plot_scores(scores, msa_length, window_size, plot_path, similarity_function,
                title=f'{msa_name}\n(across {number_of_sequences} sequences)')


if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('msa_path',
                            help='A path to an MSA file to compute homoplasy',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('homoplasy_path', help='A path to the homoplasy computation')
                            # homoplasy path does not exist in the case where a species tree is provided
        parser.add_argument('scores_output_path',
                            help='A path in which the window scores will be written to',
                            type=lambda path: path if os.path.exists(os.path.split(path)[0]) else parser.error(
                                f'output folder {os.path.split(path)[0]} does not exist!'))
        parser.add_argument('meta_output_path',
                            help='A path in which the meta output will be written to',
                            type=lambda path: path if os.path.exists(os.path.split(path)[0]) else parser.error(
                                f'output folder {os.path.split(path)[0]} does not exist!'))
        parser.add_argument('plot_path',
                            help='A path in which the scores disdribution will be plotted to',
                            type=lambda path: path if os.path.exists(os.path.split(path)[0]) else parser.error(
                                f'output folder {os.path.split(path)[0]} does not exist!'))
        parser.add_argument('window_size', type=int,
                            help='The size of a window to which a score will be computed')
        parser.add_argument('--similarity_function', default='symmetric_greedy',
                            help='The type of similarity to be computed between msa columns',
                            choices=['symmetric_greedy', 'greedy', 'mutual', 'hamming'])

        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

        args = parser.parse_args()

        compute_sweeps_score(args.msa_path, args.homoplasy_path, args.scores_output_path, args.meta_output_path,
                             args.plot_path, args.window_size, args.similarity_function)


