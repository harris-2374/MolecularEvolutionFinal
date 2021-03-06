"""
Author: Andrew Harris
Python 3.8
Version: 0.0.1
Date: 4/17/2021
"""
import argparse
from pathlib import Path

import pandas as pd
from ete3 import Tree
from ete3.coretype.tree import TreeError

########################## Helper Functions ##########################
def writeTreeFile(data, filename):
    with open(filename, 'w') as oh:
        oh.write(f"{data}\n")
    return


def writeNullOutput(nullResult, df):
    with open(nullResult, 'w') as oh:
        oh.write('Resulting data only contained null values\n')
        oh.write(f"{df}")
        return


def getFileChunks(INPUT, fileChunkOutput):
    """This function will return a dictionary with keys chunk_1 through chunk-n,
    key value is a list of lines from the file """
    currentChunk = 1
    with open(INPUT) as fh:
        geneFamily = fh.read().split("//\n")
        for count, group in enumerate(geneFamily, 1):
            if count % 10 == 0:
                print(f"-- {count:,}/{len(geneFamily):,} --")
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


def make_scientific_name_tree(tree, df):
    tree_str = str(tree)
    for pid, sn in zip(df['ProteinID'].to_list(), df['Species'].to_list()):
        tree_str = tree_str.replace(pid, sn)
        continue
    return tree_str


def make_gene_name_tree(tree, df):
    tree_str = str(tree)
    for pid, gene in zip(df['ProteinID'].to_list(), df['Gene'].to_list()):
        tree_str = tree_str.replace(pid, gene)
        continue
    return tree_str


def make_common_name_tree(tree, sdf):
    tree_str = str(tree)
    for sn, cn, in zip(sdf['scientific_name'], sdf['common_name']):
        tree_str = tree_str.replace(sn, cn)
        continue
    return tree_str


def drop_non_species_of_interest(df, speciesDF):
    indexToDrop = []
    for n, species in zip(df.index, df['Species']):
        if species not in speciesDF['scientific_name'].to_list():
            indexToDrop.append(n)
        else:
            continue
    df = df.drop(index=indexToDrop)
    return df

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
        '-o',
        '--output',
        type=str,
        action='store',
        default='./ensemblGeneFamilyResults/',
        help='',
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
        '-s',
        '--species',
        type=str,
        action='store',
        default=None,
        help='Species of interest list (single species name per line -- basic txt file)',
    )
    args = parser.parse_args()
    
    # --- Input Argparse Variables ---
    # We need a way to tie the inputs provided
    # by the user to something we can use and 
    # manipulate in our code.
    INPUT = Path(args.input)
    OUTPUT = Path(args.output)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    SPECIES = args.species
    INPUT_GENES = args.INPUT_GENE

    # Load in species file
    # I have it written to take an excel file
    # or a csv file ¯\_(ツ)_/¯ 
    try:
        species_path = Path(SPECIES)
        if 'csv' in species_path.name:
            speciesDF = pd.read_csv(species_path, sep='\t')
        elif 'xls' in species_path.name:
            speciesDF = pd.read_excel(species_path, engine='openpyxl')
    except TypeError:
        speciesDF = pd.DataFrame()

    # Load in gene name file and set into list
    genesToLookUp = [g.strip() for g in open(INPUT_GENES).readlines()]

    """
    Input file structure: (Compara.102.protein_default.nh)
        - 'SEQ' data
        - 'DATA' header
        - Newick phyloINPUT_GENEtic tree
        - '//' to indicate break to next set of data
        - Blank line
        - Next set of 'SEQ' data...
    """
    for geneOfInterest in genesToLookUp:
        print(f'--- Collecting results for {geneOfInterest} ---')
        # Step 1: Parse file into chunks
        fileChunkOutput = OUTPUT / f'{geneOfInterest}'
        currentChunk = 1
        with open(INPUT) as fh:
            geneFamily = fh.read().split("//\n")
            for group in geneFamily:
                seqlines = [l for l in group.split("\n") if l != '']
                for i, l in enumerate(reversed(seqlines), 1):
                    if not l: # Skips blank lines if not caught
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
                            df = df.drop(labels=['SEQ', 'gain-loss?', 'Start', 'Stop', 'Chromosome'], axis=1)
                            CommonNameConvertDict = {sn:cn for sn,cn in zip(speciesDF['scientific_name'], speciesDF['common_name'])}
                            CommonNameConvertDict['None'] = None
                            df['CommonName'] = [CommonNameConvertDict.get(str(g)) for g in df['Species']]
                            OrderConvertDict = {sn:order for sn,order in zip(speciesDF['scientific_name'], speciesDF['order'])}
                            df['Order'] = [OrderConvertDict.get(str(g)) for g in df['Species']]
                            df['Chunk'] = [f'{geneOfInterest}_chunk_{currentChunk}']*len(df)
                        except:
                            # print(seqdata)
                            pass
                        df['Gene'] = df['Gene'].apply(lambda x: 'NULL' if not x else x) # Replace None with "NULL", None throws error
                        for gene in df['Gene']:
                            if str(geneOfInterest).lower() == str(gene).lower():
                                print(f"Data found for {geneOfInterest}")
                                # Output File Names
                                speciesCountOutput = OUTPUT / 'speciesCounts'
                                speciesCountOutput.mkdir(parents=True, exist_ok=True)  # Create output directory
                                currChunkDir = fileChunkOutput / f"chunk_{currentChunk}"
                                currChunkDir.mkdir(parents=True, exist_ok=True)  # Create output directory
                                currChunkSeqFile = currChunkDir / f"chunk_{currentChunk}_SEQ.tsv"
                                currChunkPidTreeFile = currChunkDir / f"chunk_{currentChunk}_ProteinID_Newick.tree"
                                currChunkSciNameTreeFile = currChunkDir / f"chunk_{currentChunk}_ScientificName_Newick.tree"
                                currChunkCommonNameTreeFile = currChunkDir / f"chunk_{currentChunk}_CommonName_Newick.tree"
                                currChunkGeneTreeFile = currChunkDir / f"chunk_{currentChunk}_GeneName_Newick.tree"
                                CountSummaryOutputFileName = speciesCountOutput / f'{geneOfInterest}_copy_number_summary_chunk_{currentChunk}.txt'
                                orderGeneCount = speciesCountOutput / f'{geneOfInterest}_OrderCounts_chunk_{currentChunk}.txt'
                                nullResult = currChunkDir / 'null_result.txt'
                                malformedTree = currChunkDir / 'malformed_tree.txt'
                                # Filter out non-species of interest entries
                                if speciesDF.empty:
                                    df.to_csv(currChunkSeqFile, sep="\t", index=False)
                                    tree = Tree(tree)
                                    tree.prune(df['ProteinID'].to_list(), preserve_branch_length=True)
                                    # Convert tree to scientific names
                                    scientificNameTree = make_scientific_name_tree(Tree(tree), df)
                                    geneNameTree = make_gene_name_tree(Tree(tree), df)

                                    # Organize data into columns [geneID, Speceies, Order, SpeciesCopyNum, OrderCopyNum]
                                    species_count_dict = {f:df['Species'].to_list().count(f) for f in df['Species'].unique()}
                                    df['SpeciesCopyNumber'] = [species_count_dict[f] for f in df['Species']]
                                    
                                    OrderCounts = {f:df['Order'].to_list().count(f) for f in df['Order'].unique()}
                                    df['OrderCopyNumber'] = [OrderCounts[f] for f in df['Order']]
                                    df = df[['Gene', 'Species', 'CommonName', 'Order', 'SpeciesCopyNumber', 'OrderCopyNumber', 'Chunk', 'GeneID', 'ProteinID']]  # Reorder columns
                                    # Output files
                                    writeTreeFile(Tree(tree).prune(df['ProteinID'].to_list(), preserve_branch_length=True), currChunkPidTreeFile)
                                    writeTreeFile(scientificNameTree, currChunkSciNameTreeFile)
                                    writeTreeFile(geneNameTree, currChunkGeneTreeFile)
                                    df.to_csv(CountSummaryOutputFileName, sep='\t', index=False)
                                    df.to_csv(currChunkSeqFile, sep="\t", index=False)
                                    currentChunk += 1
                                    break
                                else:
                                    # Remove species that are not in species of interest file
                                    df = drop_non_species_of_interest(df, speciesDF)
                                    # Sort by Order
                                    df = df.sort_values(by='Order')
                                    # If all species have NULL as gene, output null result file
                                    if (len(df['Gene'].unique()) == 1) and (df['Gene'].unique()[0] == 'NULL'):
                                        writeNullOutput(nullResult, df)
                                        break
                                    # This checks to make sure the newick tree is valid,
                                    # if not then it will return a file saying the tree
                                    # is malformed 
                                    try:
                                        tree = Tree(tree)
                                        tree.prune(df['ProteinID'].to_list(), preserve_branch_length=True)
                                    except TreeError:
                                        writeTreeFile(tree, malformedTree)
                                        currentChunk += 1
                                        break
                                    
                                    # Convert tree leaves to scientific, common, and gene names
                                    scientificNameTree = make_scientific_name_tree(tree.write(), df)
                                    CommonNameTree = make_common_name_tree(scientificNameTree, speciesDF)
                                    geneNameTree = make_gene_name_tree(tree.write(), df)

                                    # Organize data into columns [geneID, Speceies, Order, SpeciesCopyNum, OrderCopyNum]
                                    species_count_dict = {f:df['Species'].to_list().count(f) for f in df['Species'].unique()}
                                    df['SpeciesCopyNumber'] = [species_count_dict[f] for f in df['Species']]
                                    
                                    OrderCounts = {f:df['Order'].to_list().count(f) for f in df['Order'].unique()}
                                    df['OrderCopyNumber'] = [OrderCounts[f] for f in df['Order']]

                                    # Reorder columns
                                    df = df[['Gene', 'Species', 'CommonName', 'Order', 'SpeciesCopyNumber', 'OrderCopyNumber', 'Chunk', 'GeneID', 'ProteinID']]

                                    # Output all files
                                    writeTreeFile(tree.write(), currChunkPidTreeFile)
                                    writeTreeFile(scientificNameTree, currChunkSciNameTreeFile)
                                    writeTreeFile(CommonNameTree, currChunkCommonNameTreeFile)
                                    writeTreeFile(geneNameTree, currChunkGeneTreeFile)
                                    df.to_csv(CountSummaryOutputFileName, sep='\t', index=False)

                                    # Output null file if no data present
                                    if df.empty:
                                        with open(nullResult, 'w') as oh:
                                            oh.write('No data available')
                                            currentChunk += 1
                                            break
                                    else:
                                        df.to_csv(currChunkSeqFile, sep="\t", index=False)
                                        currentChunk += 1
                                        break
                            else:
                                continue
                        break
                # File tracker
                # if count % 10000 == 0:
                #     print(f"-- {count:,}/{len(geneFamily):,} --")
                continue
    # This section of code will output all species count data into a single file
    countDir = OUTPUT / 'speciesCounts'
    countFiles = [f for f in countDir.iterdir() if f.is_file()]
    cumulative_df = pd.DataFrame(columns=['Gene', 'Species', 'CommonName', 'Order', 'SpeciesCopyNumber', 'OrderCopyNumber', 'Chunk', 'GeneID'])
    cumulativeFileName = OUTPUT / 'CumulativeSpeciesCounts.tsv'

    for f in countFiles:
        fdf = pd.read_csv(f, sep='\t')
        cumulative_df = pd.concat([cumulative_df, fdf])
        continue

    cumulative_df.to_csv(cumulativeFileName, sep='\t', index=False)
    
    return

if __name__ == "__main__":
    main()