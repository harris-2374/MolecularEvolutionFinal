"""
Author: Andrew Harris
Python 3.8
Version: 0.0.1
Last update: 2/14/2021
Summary: This script will count the number of genes per node and
add it to the Newick file.
"""
import argparse
from pathlib import Path

## Dependencies
from ete3 import Tree


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
        default='./',
        help='',
    )
    args = parser.parse_args()
    
    # --- Input Argparse Variables ---
    """
    We need a way to tie the inputs provided
    by the user to something we can use and 
    manipulate in our code.
    """
    INPUT= Path(args.input)
    OUTPUT= Path(args.output)

    """
    This snipit does a few things at once.
    First, it opens the INPUT file -- open() --
    Since I know that this Newick file only has
    one line that is the tree, I can access the first
    line of the file by adding [0] to the end of readlines(), 
    The .strip() part just removs the "\n" at the end each line
    """
    newick = open(INPUT).readlines()[0].strip()

    """
    This is a class imported from ete3. 
    We put our Newick tree into the Tree class,
    this allows us to traverse our 
    tree in a programatic way
    """
    tree = Tree(newick)  
    
    """
    Move through the tree node-by-node
    postorderly and add the number of leafs (genes)
    that are attached to the current node to it's features.
    This changes our standard Newick file to an extended
    form of Newick called New Hampshire eXtended format (NHX).
    """
    for node in tree.traverse("postorder"):
        # print(node.write())
        # print(len(node.get_leaves()))
        node.add_feature('numGenes', len(node.get_leaves()))
        continue
    
    """
    We then write the updated Newick file to a new file,
    explicitly telling the write function to include our
    calculated feature 'numGenes' in the Newick tree. 
    """
    tree.write(outfile='nodeGeneCount_NHX.tree', features=['numGenes'])
    
    """
    This just prints the ascii version of the tree.
    It will be pretty big, so may not be userful. 
    """
    print(tree)

    """
    If you want to look at this in a user interface,
    remove the "#" in front of tree.show() then run 
    in your terminal with the virtual enviornment on,
    'pip install PyQT5'. 
    """
    #tree.show()    
    return

if __name__ == "__main__":
    main()
