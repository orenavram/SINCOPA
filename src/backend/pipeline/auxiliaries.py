from Bio import Phylo


def get_tree_labels(tree_to_prune_path):
    return [node.name for node in Phylo.read(tree_to_prune_path, 'newick').get_terminals()]


def fail(error_msg, error_file_path):
    with open(error_file_path, 'w') as error_f:
        error_f.write(error_msg + '\n')
    raise Exception(error_msg)


def update_html(html_path, src, dst):
    # The initial file exists (generate by the cgi) so we can read and parse it.
    with open(html_path) as f:
        html_content = f.read()
    html_content = html_content.replace(src, dst)
    with open(html_path, 'w') as f:
        f.write(html_content)


def append_to_html(html_path, new_content):
    with open(html_path) as f:
        html_content = f.read()
    html_content += new_content
    with open(html_path, 'w') as f:
        f.write(html_content)


def load_header2sequences_dict(fasta_path, get_length=False, upper_sequence=False):
    header_to_sequence_dict = {}
    seq_length = 0

    with open(fasta_path) as f:
        header = f.readline().lstrip('>').rstrip()
        sequence = ''
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                seq_length = len(sequence)
                if upper_sequence:
                    header_to_sequence_dict[header] = sequence.upper()
                else:
                    # leave untouched
                    header_to_sequence_dict[header] = sequence
                header = line.lstrip('>')
                sequence = ''
            else:
                sequence += line

        # don't forget last record!!
        if sequence != '':

            if upper_sequence:
                header_to_sequence_dict[header] = sequence.upper()
            else:
                # leave untouched
                header_to_sequence_dict[header] = sequence

    if get_length:
        return header_to_sequence_dict, seq_length
    else:
        return header_to_sequence_dict



