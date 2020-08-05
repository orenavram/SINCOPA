import pandas as pd
import os
import wget
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

for path, folders, files in os.walk('/Users/Oren/Dropbox/Projects/sincopa/data/'):
    for file in files:
        if file != 'assembly_summary.txt':
            continue

        species = os.path.split(path)[-1]
        file_path = os.path.join(path, file)
        out_path = os.path.join(path, 'assembly_summary_complete_genome.txt')
        #
        # path = '/Users/Oren/Dropbox/Projects/sincopa/genomic_data/Chlamydia_trachomatis'
        # file_path = f'{path}/assembly_summary.txt'
        # out_path = f'{path}/assembly_summary_complete_genome.txt'
        meta_df = pd.read_csv(file_path, skiprows=1, sep='\t')
        meta_df = meta_df[meta_df.assembly_level == 'Complete Genome']  # discard instances of incomplete genomes
        meta_df.to_csv(out_path, sep='\t')
        logger.info(f'\n\nGetting {species} genomes')
        path = f'{path}/genomes'
        os.makedirs(path, exist_ok=True)
        logger.info('The following paths were stored:')
        for i, ftp_path in enumerate(meta_df.ftp_path):
            local_path = f'{path}/{os.path.split(ftp_path)[-1]}_genomic.fna.gz'
            wget.download(f'{ftp_path}/{os.path.split(ftp_path)[-1]}_genomic.fna.gz', local_path)
            logger.info(local_path)
            # if i == 10:
            #     break




