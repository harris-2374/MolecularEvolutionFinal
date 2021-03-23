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
    SOI = Species of interest

    Outline for species homology file:
        1. Load files into pandas dataframeframes
        2. Filter dataframeframe (see filter section below)
        3. Group dataframe by 'gene_stable_id'
        4. Search for genes from genes-of-interest file
        5. Output information with or w/out all speices (only keep species-of-interest if SOI file provided)

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
        '-g',
        '--genes',
        type=str,
        action='store',
        default=None,
        help='Gene of interest',
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
    SPECIES = args.species
    GENES = Path(args.genes)
    PERCENTID = int(args.p)
    logger.info(f"File Name: {INPUT.name}")

    OUTPUT.mkdir(parents=True, exist_ok=True)

    # Step 1: Load dataframe into dataframeframes
    # sep='\t' tells pandas it is a tab delimited file (tsv)
    df = pd.read_csv(INPUT, sep='\t')
    geneDF = pd.read_csv(GENES, names=['geneID', 'geneName'], sep='\t')
    logger.info(f'Number of genes to search for: {len(geneDF)}')
    try:
        species_path = Path(SPECIES)
        if 'csv' in species_path.name:
            speciesDF = pd.read_csv(species_path, sep='\t')
        elif 'xls' in species_path.name:
            speciesDF = pd.read_excel(species_path, engine='openpyxl')
    except TypeError:
        speciesDF = None

    # Step 2: Filter dataframeFrame
    length_before_filter = len(df)
    
    one2many = df[(df['identity'] >= PERCENTID) & (df['homology_type'] == 'ortholog_one2many')]
    one2one = df[(df['identity'] >= PERCENTID) & (df['homology_type'] == 'ortholog_one2one')]
    
    filtered_df = pd.concat([one2one, one2many])
    
    logger.info(f"Intial number of entries: {length_before_filter:,}")
    logger.info(f"Remaining number of entries after filter: {len(filtered_df):,} ({len(filtered_df)-length_before_filter:,})")

    # Step 3: Group dataframe by 'gene_stable_id'
    grouped_df = filtered_df.groupby(by=['gene_stable_id'])
    logger.info(f"Number of genes: {len(grouped_df):,}")
    
    # Step 4: Find genes of interest and output dataframeframe
    # to a file.
    for gene, dataframe in grouped_df:
        if gene in geneDF['geneID'].to_list():
            geneINFO = geneDF[geneDF['geneID'] == gene]  # Remove all dataframeframe entries besides the current gene of interest
            geneINFO = geneINFO.reset_index(drop=True)   # Reset index so gene is at index position 0
            assert len(geneINFO) == 1  # Make sure gene is not added to list twice
            geneName = geneINFO['geneName'][0]
            outFilename = OUTPUT / f'{geneName}_{gene}_results.tsv'
            if not SPECIES:
                # Output dataframe with all dataframe
                logger.info("===========================================")
                logger.info(f"Gene name: {gene}")
                logger.info(f"Number of entries: {len(dataframe)}")
                logger.info(f"Number of unique species: {len(dataframe['homology_species'].unique())}")
                logger.info(dataframe)
                logger.info("===========================================")
                # Write results to file
                dataframe.to_csv(outFilename, sep='\t', index=False)
            else:
                # Remove non-species of interest dataframe
                filteredDF = pd.DataFrame()
                uniqueNames = dataframe['homology_species'].drop_duplicates()
                namesToRemove = pd.concat([uniqueNames, speciesDF['scientific_name']])
                namesToRemove = namesToRemove.drop_duplicates(keep=False)
                for s in namesToRemove:
                    if filteredDF.empty:
                        filteredDF = dataframe[dataframe.homology_species != s]
                    else:
                        filteredDF = filteredDF[filteredDF.homology_species != s]
                # Now add in common names
                conversionDict = {i:j for i,j in zip(speciesDF['scientific_name'], speciesDF['common_name'])}
                filteredDF['common_name'] = filteredDF['homology_species'].apply(lambda x: conversionDict[x])
                logger.info("===========================================")
                logger.info(f"Gene name: {gene}")
                logger.info(f"Number of unique species: {len(dataframe['homology_species'].unique())}")
                logger.info(f"Number of entries after filtering: {len(filteredDF)}")
                logger.info(filteredDF)
                logger.info("===========================================")
                filteredDF.to_csv(outFilename, sep='\t', index=False)
    # Exit script by removing output directory if it is empty
    # If output directory is not empty, just exit the script
    outputFiles = [f for f in OUTPUT.iterdir()]
    if len(outputFiles) > 0:
        return
    else:
        OUTPUT.rmdir()
        return

if __name__ == "__main__":
    main()

