"""
Author: Andrew Harris
Python 3.8 

Dependency install commands
conda install ensembl-rest or pip install ensembl-rest
conda install pandas or pip install pandas
"""
import argparse
import logging
import os
from pathlib import Path
import textwrap

## Dependencies
import ensembl_rest
import pandas as pd

########################### Main Function ###########################
def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        action='store',
        required=True,
        # default='ensembl_gene_lookup_input.tsv',
        help="Excel/tsv file with gene name, common name, GeneID, and ProteinID",
    )
    parser.add_argument(
        '-g',
        '--gene',
        type=str,
        action='store',
        # default='ensembl_gene_lookup_input.tsv',
        help="Ensembl gene name to look up",
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        action='store',
        # default='./ensembl_sequences',
        help='',
    )
    args = parser.parse_args()
    
    # --- Input Argparse Variables ---
    INPUT= Path(args.input)
    GENE = str(args.gene)
    OUTPUT= Path(args.output)
    print(GENE, OUTPUT)
    GENE_OUTPUT = OUTPUT / 'GeneID' / GENE
    PROTEIN_OUTPUT = OUTPUT / 'ProteinID' / GENE
    GENE_LOG_OUTPUT = OUTPUT / 'gene_logs/'
    PROTEIN_LOG_OUTPUT = OUTPUT / 'protein_logs/'
    GENE_LOG_OUTPUT_FILE = GENE_LOG_OUTPUT / f'{GENE}_fetchEnsemblSequences_GeneID.log'
    PROTEIN_LOG_OUTPUT_FILE = PROTEIN_LOG_OUTPUT / f'{GENE}_fetchEnsemblSequences_ProteinID.log'

    # Clear log file if it exist
    if os.path.exists(GENE_LOG_OUTPUT_FILE):
        os.remove(GENE_LOG_OUTPUT_FILE)
    if os.path.exists(PROTEIN_LOG_OUTPUT_FILE):
        os.remove(PROTEIN_LOG_OUTPUT_FILE)

    # Make output directory
    GENE_OUTPUT.mkdir(parents=True, exist_ok=True)
    GENE_LOG_OUTPUT.mkdir(parents=True, exist_ok=True)
    PROTEIN_OUTPUT.mkdir(parents=True, exist_ok=True)
    PROTEIN_LOG_OUTPUT.mkdir(parents=True, exist_ok=True)

    # Set up logging output files + handlers + streams
    formatter = logging.Formatter('%(message)s')
    # Gene logging
    gene_logger = logging.getLogger(GENE_LOG_OUTPUT_FILE.as_posix())
    gene_logger.setLevel(logging.INFO)
    gene_stream_handler = logging.StreamHandler()
    gene_file_handler = logging.FileHandler(GENE_LOG_OUTPUT_FILE)
    gene_file_handler.setFormatter(formatter)
    gene_logger.addHandler(gene_file_handler)
    gene_logger.addHandler(gene_stream_handler)
    # Protein logging
    protein_logger = logging.getLogger(PROTEIN_LOG_OUTPUT_FILE.as_posix())
    protein_logger.setLevel(logging.INFO)
    protein_stream_handler = logging.StreamHandler()
    protein_file_handler = logging.FileHandler(PROTEIN_LOG_OUTPUT_FILE)
    protein_file_handler.setFormatter(formatter)
    protein_logger.addHandler(protein_file_handler)
    protein_logger.addHandler(protein_stream_handler)

    # Load input file into df + fill na values with str(NULL)
    df = pd.read_csv(INPUT, sep='\t')
    df = df.fillna('NULL')
    # Fetch information + output + log
    for row in df.itertuples(index=False):
        if 'NULL' in row.Gene:
            geneName = 'NULL_Novel'
        else:
            geneName = row.Gene
        # Row attributes
        commonName = row.CommonName
        sampleOrder = row.Order
        geneID = row.GeneID
        proteinID = row.ProteinID
        # Output filenames
        geneOutputFileName = GENE_OUTPUT / f"{geneName}_{commonName}_{sampleOrder}_{geneID}.fasta"
        proteinOutputFileName = PROTEIN_OUTPUT / f"{geneName}_{commonName}_{sampleOrder}_{proteinID}.fasta"
        # Fetch gene sequences
        try:
            gene_fasta_results = ensembl_rest.sequence_id(geneID)
            with open(geneOutputFileName, 'w') as oh:
                oh.write(f">{geneOutputFileName.stem}\n")
                oh.write(textwrap.fill(gene_fasta_results['seq'], width=80))
            gene_logger.info('----------------------------------------------------------------')
            gene_logger.info(f"CommonName: {commonName}\nGeneID: {geneID}\nSeqFile: {geneOutputFileName}\nDescription: {gene_fasta_results['desc']}")
        except:
            gene_logger.info('----------------------------------------------------------------')
            gene_logger.info(f'Could not find info for {geneID}')
            pass
        # Fetch protein sequences
        try:
            protein_fasta_results = ensembl_rest.sequence_id(proteinID)
            with open(proteinOutputFileName, 'w') as oh:
                oh.write(f">{proteinOutputFileName.stem}\n")
                oh.write(textwrap.fill(protein_fasta_results['seq'], width=80))
            protein_logger.info('----------------------------------------------------------------')
            protein_logger.info(f"CommonName: {commonName}\nProteinID: {proteinID}\nSeqFile: {proteinOutputFileName}\nDescription: {protein_fasta_results['desc']}")
        except:
            protein_logger.info('----------------------------------------------------------------')
            protein_logger.info(f'Could not find info for {proteinID}')
            pass
        continue
    return

if __name__ == "__main__":
    main()
    
