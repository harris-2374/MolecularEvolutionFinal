"""
Author: Andrew Harris
Python 3.8
Version: 0.0.1
Last update: 2-23-2021
"""
import argparse
import logging
import os
from pathlib import Path

## Dependencies
import pandas as pd

"""
    Outline:
        1. Load files into pandas dataframes
        2. Filter dataframe
        3. Group data by 'gene_stable_id'
        
    Filters:
        1. Percent identity >= 50%
        2. 'homology_type' == 'ortholog_one2many' or 'ortholog_one2one'
        3. 
"""

########################### Logging Init ###########################
# If log file exist, erase it
if os.path.exists('logs/speciesHomologyParser.log'):
    os.remove('logs/speciesHomologyParser.log')

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s: %(message)s')
Path("logs/").mkdir(parents=True, exist_ok=True)
file_handler = logging.FileHandler('logs/speciesHomologyParser.log')
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
        required=True,
        help="Homology file path",
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        action='store',
        default='./',
        help='',
    )
    parser.add_argument(
        '-s',
        '--species',
        type=str,
        action='store',
        default=None,
        help='Species of interest list (single species name per line -- basic txt file)',
    )
    parser.add_argument(
        '-p',
        type=int,
        action='store',
        default=50,
        help='Percent identity cutoff',
    )
    args = parser.parse_args()
    
    # --- Input Argparse Variables ---
    INPUT = Path(args.input)
    OUTPUT = Path(args.output)
    SPECIES = Path(args.species)
    PERCENTID = int(args.p)
    logger.info(f"File Name: {INPUT.name}")

    # Step 1: Load into dataframe
    df = pd.read_csv(INPUT, sep='\t')  # sep='\t' tells pandas it is a tab delimited file (tsv)

    # Step 2: Filter DataFrame
    length_before_filter = len(df)
    
    one2many = df[(df['identity'] >= PERCENTID) & (df['homology_type'] == 'ortholog_one2many')]
    one2one = df[(df['identity'] >= PERCENTID) & (df['homology_type'] == 'ortholog_one2one')]
    
    filtered_df = pd.concat([one2one, one2many])
    
    logger.info(f"Intial number of entries: {length_before_filter:,}")
    logger.info(f"Remaining number of entries after filter: {len(filtered_df):,} ({len(filtered_df)-length_before_filter:,})")

    # Step 3: Group data by 'gene_stable_id'
    grouped_df = filtered_df.groupby(by=['gene_stable_id'])
    logger.info(f"Number of genes: {len(grouped_df):,}")
    
    # Step 4: Run through each gene_stable_id
    tally = 0  # Testing
    for gene, data in grouped_df:
        print("\n===========================================")
        print(f"Gene name: {gene}")
        print(f"Number of entries: {len(data)}")
        print(data)
        print(data.homology_species)
        print("===========================================\n")
        
        # Testing
        if tally == 25:
            break
        tally += 1

    return

if __name__ == "__main__":
    main()

