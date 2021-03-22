"""
Author: Andrew Harris
Python 3.8
Version: 0.0.1
Date: 2/14/2021
"""
import argparse
from pathlib import Path

import pandas as pd

########################## Helper Functions ##########################
def writeTreeFile(data, filename):
    with open(filename, 'w') as oh:
        oh.write(f"{data}\n")
    return

def writeExcelSeqFile(data, filename):
    df = pd.DataFrame()
    index = 0
    for l in data:
        df.loc[index] = l
        index += 1
        continue
    df.to_excel(filename, index=False)
    return


def getFileChunks(INPUT, fileChunkOutput):
    """This function will return a dictionary with keys chunk_1 through chunk-n,
    key value is a list of lines from the file """
    currentChunk = 1
    with open(INPUT) as fh:
        geneGroups = fh.read().split("//\n")
        for count, group in enumerate(geneGroups, 1):
            if count % 10 == 0:
                print(f"-- {count:,}/{len(geneGroups):,} --")
            seqlines = [l for l in group.split("\n") if l != '']
            currChunkDir = fileChunkOutput / f"chunk_{currentChunk}"
            currChunkDir.mkdir(parents=True, exist_ok=True)  # Make output directory
            currChunkSeqFile = currChunkDir / f"chunk_{currentChunk}_SEQ.tsv"
            currChunkTreeFile = currChunkDir / f"chunk_{currentChunk}_Newick.tree"
            for i, l in enumerate(reversed(seqlines), 1):
                if not l: # Skips blank lines
                    continue
                elif ';' in l:  # Write out newick file
                    writeTreeFile(l, currChunkTreeFile)
                    continue
                elif "SEQ" in l:
                    seqdata = [l.split(" ") for l in seqlines[:len(seqlines)-i]]
                    avglen = sum([len(i) for i in seqdata]) / len(seqdata)
                    print(avglen, max([len(i) for i in seqdata]), min([len(i) for i in seqdata]))

                    for l in seqdata:
                        if len(l) > 9:
                            print(l)

                    try:
                        df = pd.DataFrame(data=seqdata, columns=range(len(seqdata[0])))
                    except ValueError:
                        print(seqdata)
                        raise ValueError
                    # df.to_csv(currChunkSeqFile, sep="\t", index=False)
                    break
                # else:
                #     break
            currentChunk += 1
            continue
    return


########################### Main Function ###########################
def main():
    # These are the flags the user will use to tie input
    # values like file paths and parameter values.
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        action='store',
        required=True,
        help="",
    )
    parser.add_argument(
        '-g',
        '--INPUT_GENE',
        type=str,
        action='store',
        required=True,
        help="Comma separated list of INPUT_GENEs to search for [Replace spaces with underscores]",
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        action='store',
        default='./ensemblGeneFamily/',
        help='',
    )
    args = parser.parse_args()
    
    # --- Input Argparse Variables ---
    # We need a way to tie the inputs provided
    # by the user to something we can use and 
    # manipulate in our code.
    INPUT = Path(args.input)
    INPUT_GENE = args.INPUT_GENE
    OUTPUT = Path(args.output)
    OUTPUT.mkdir(parents=True, exist_ok=True)

    """
    Input file structure: (Compara.102.protein_default.nh)
        - 'SEQ' data
        - 'DATA' header
        - Newick phyloINPUT_GENEtic tree
        - '//' to indicate break to next set of data
        - Blank line
        - Next set of 'SEQ' data...
    """
    # Step 1: Parse file into chunks
    fileChunkOutput = OUTPUT / f'{INPUT_GENE}'
    currentChunk = 1
    with open(INPUT) as fh:
        geneGroups = fh.read().split("//\n")
        for count, group in enumerate(geneGroups, 1):
            seqlines = [l for l in group.split("\n") if l != '']
            for i, l in enumerate(reversed(seqlines), 1):
                if not l: # Skips blank lines
                    continue
                elif ';' in l:
                    tree = l
                    continue
                elif "SEQ" in l:
                    seqdata = [l.split(" ") for l in seqlines[:len(seqlines)-i]]
                    # This section fixes an issue where INPUT_GENEs have
                    # spaces in their names at the end of the row
                    # This joins the name by replacing the space
                    # with an underscore 
                    for l in seqdata:
                        if len(l) > 9:
                            l[8:len(l)+1] = ['_'.join(l[8:len(l)+1])] 
                    try:
                        df = pd.DataFrame(data=seqdata, columns=['SEQ', 'Species', 'ProteinID', 'Chromosome', 'Start', 'Stop', 'gain-loss?', 'GeneID', 'Gene'])
                    except:
                        print(seqdata)
                    df['Gene'] = df['Gene'].apply(lambda x: 'NULL' if not x else x) # Replace None with "NULL", None throws error
                    nullCount = df['Gene'].to_list().count('NULL')
                    for n, gene in enumerate(df['Gene']):
                        if INPUT_GENE.lower() in gene.lower():
                            currChunkDir = fileChunkOutput / f"chunk_{currentChunk}"
                            currChunkDir.mkdir(parents=True, exist_ok=True)  # Make output directory
                            currChunkSeqFile = currChunkDir / f"chunk_{currentChunk}_SEQ.tsv"
                            currChunkTreeFile = currChunkDir / f"chunk_{currentChunk}_Newick.tree"
                            df.to_csv(currChunkSeqFile, sep="\t", index=False)
                            writeTreeFile(tree, currChunkTreeFile)
                            currentChunk += 1
                            break
                    break
            if count % 1000 == 0:
                print(f"-- {count:,}/{len(geneGroups):,} --")
            continue
    return

if __name__ == "__main__":
    main()