"""
Author: Andrew Harris
Python 3.8
Version: 0.0.1
Last update: 3/27/2021
"""
import argparse
import logging
from pathlib import Path

## Dependencies
import pandas as pd

########################### Logging Init ###########################
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s: %(message)s')
Path("logs/").mkdir(parents=True, exist_ok=True)
file_handler = logging.FileHandler('logs/parseActinSupplementFile.log')
stream_handler = logging.StreamHandler()
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.addHandler(stream_handler)


########################## Helper Functions ##########################



########################### Main Function ###########################
def main():
    # parser = argparse.ArgumentParser(description='')
    # parser.add_argument(
    #     '-i',
    #     '--input',
    #     type=str,
    #     action='store',
    #     required=True,
    #     help="",
    # )
    # parser.add_argument(
    #     '-o',
    #     '--output',
    #     type=str,
    #     action='store',
    #     default='./',
    #     help='',
    # )
    # args = parser.parse_args()
    
    # --- Input Argparse Variables ---
    # INPUT = Path(args.input)
    # OUTPUT= Path(args.output)
    INPUT = Path("ActinGeneFamilyGeneNames.xlsx")
    OUTPUT = Path("genesToLookUp.txt")
    

    # --- Load file into Pandas dataframe ---
    df = pd.read_excel(INPUT, engine='openpyxl')

    # --- Collect unique gene names ---
    # Most genes names have a species identifier at
    # the end (gene_species). Just keep gene name. 
    
    geneNames = []
    for row in df.itertuples(index=False):
        if '_' in row.gene1:
            gene1 = row.gene1[:-4]
        else:
            gene1 = row.gene1
            
        if '_' in row.gene2:
            gene2 = row.gene2[:-4]
        else:
            gene2 = row.gene2

        geneNames.append(gene1)
        geneNames.append(gene2)
        continue
    finalGeneList = pd.Series(list(set(geneNames)))
    finalGeneList.to_csv(OUTPUT, index=False, header=False)
    logger.info(f"Number of output genes:{len(finalGeneList)}")


    return

if __name__ == "__main__":
    main()

