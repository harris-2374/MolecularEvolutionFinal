import argparse
from pathlib import Path
from Bio import AlignIO
import pandas as pd
"""
This script will take in multi-alginment fasta file and will locate positions 
where a single sample has a change in comparison to all other samples. 

The script will ignore a position if >50% of the samples have "-" or missing data
If the script identifies the unique amino acid to be "-" it will ignore it

The user will ultimately need to prune the output dataframe to remove hits for 
non-species of interst. (primates for our case)
"""

def main():
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
    args = parser.parse_args()
    INPUT = Path(args.input)
    OUTPUT = Path(args.output)    
    
    aln = AlignIO.read(open(INPUT), "fasta")
    diff_df = pd.DataFrame(columns=['seqName', 'position', 'uniqueAA', 'allSampleBases'])
    df_idx = 0
    for i in range(len(aln[0].seq)):
        pos_data = aln[:, i]
        sample_names = [n.id for n in aln[:, i:i]]
        if (pos_data.count("-")/len(pos_data)) > 0.5:
            continue
        elif len(list(set(pos_data))) > 2:
            continue
        else:
            aa_to_find = [a for a in pos_data if (pos_data.count(a) == 1) and (pos_data != '-')]
            if len(aa_to_find) == 0:
                continue
            aa_to_find_pos = pos_data.rindex(aa_to_find[0])
            if pos_data[aa_to_find_pos] == '-':
                continue
            elif 'primate' in sample_names[aa_to_find_pos]:
                continue
            print(f"Unique AA found! -- PositionInAlinment:{i} -- ChangedAminoAcid: {pos_data[aa_to_find_pos]} -- SampleName:{sample_names[aa_to_find_pos]}")
            diff_df.at[df_idx, 'seqName'] = sample_names[aa_to_find_pos]
            diff_df.at[df_idx, 'position'] = i
            diff_df.at[df_idx, 'uniqueAA'] = pos_data[aa_to_find_pos]
            diff_df.at[df_idx, 'allSampleBases'] = pos_data
            df_idx += 1
            continue
    if len(diff_df) == 0:
        OUTPUT = OUTPUT.parents[0] / f"{OUTPUT.stem}_NO_RESULTS.xlsx"

    print(diff_df)
    diff_df.to_excel(OUTPUT, index=False)
    return

if __name__ == '__main__':
    main()


