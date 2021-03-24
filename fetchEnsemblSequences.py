"""
Author: Andrew Harris
Python 3.8
Version: 
Last update: 
Summary: 
"""
import argparse
import logging
from pathlib import Path

## Dependencies
import ensembl_rest
import pandas as pd

########################### Logging Init ###########################
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s: %(message)s')
Path("logs/").mkdir(parents=True, exist_ok=True)
file_handler = logging.FileHandler('logs/fetchEnsemblSequences.log')
stream_handler = logging.StreamHandler()
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.addHandler(stream_handler)


########################## Helper Functions ##########################



########################### Main Function ###########################
def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        action='store',
        # required=True,
        default='ensembl_gene_lookup_input.tsv',
        help="Ensembl gene ID's",
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        action='store',
        default='./ensembl_sequences',
        help='',
    )
    args = parser.parse_args()
    
    # --- Input Argparse Variables ---
    INPUT= Path(args.input)
    OUTPUT= Path(args.output)

    # Make output directory
    OUTPUT.mkdir(parents=True, exist_ok=True)
    # Collect geneID's
    # File headers == ['gene_id', 'common_name']
    genesToFind = pd.read_csv(INPUT, sep='\t')

    for row in genesToFind.itertuples(index=False):
        geneID, commonName = row
        outFileName = OUTPUT / f"{commonName}_{geneID}.fasta"
        fasta_seq = ensembl_rest.sequence_id(geneID, headers={'content-type': 'text/x-fasta'})
        with open(outFileName, 'w') as oh:
            oh.write(fasta_seq)
        logger.info(f"CommonName: {commonName}\nGeneID: {geneID}\nSeqFile: {outFileName}")
        logger.info('----------------------------------------------------------------')


    return

if __name__ == "__main__":
    main()

