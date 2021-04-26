# MolecularEvolutionFinal
 Scripts for GENE648 Molecular Evolution final project
 
 ## ensemblGeneFamilyHomologyFileParser.py
    Summary:
        This script takes in an Ensembl gene family emf file and returns
        information for a specific genes. The resulting information contains
        protein and gene ID's along with other metadata.If a species
        of interest file is provided, the script will remove all information
        from any species not found in the lookup file.
    
    Input:
        1. Ensembl emf file (required)
        2. Gene of interest (required)
        3. Species of interest lookup file (optional) 
            - Tab delimited columns of scientific_name, common_name, and enseml_ID

 ## fetchEnsemblSequences.py
    Summary:
        This script utilizes the python Ensembl rest api to fetch fasta 
        formatted nucleotide sequences and writes them to fasta file 
        in a provided output directory. The input file consist of two tab 
        delimeted columns where the first column is the Ensembl gene ID and 
        the other is a species identifier. 

    Input:
        1. Gene ID and species identifier file (required)

 ## identifyVariantsOfInterest.py
    Summary:
        This script takes in a multi-alignment fasta file (clustal output in our case)
        and it identifies bases where a single carnivore has a different amino acid compared
        to all other species.
    
    Filtration Steps:
        1. Ignores positions where greater than 50% of the samples have missing data ("-")
        2. Ignores pisitions where more than two different amino acids are present
        3. Ignores results where the variant of interest in is a primate (study specific -- can be changed)
    
    Input:
        Multi-alignment fasta file
    
    Output:
        Excel file with fasta header of sample with variant of interest, the position in the alignment, 
        the altered amino acid and a string of all bases from each sample to verify script has worked correctly.  

 ## speciesHomologyParser.py (Ultimately not used in study)
    Summary:
        This script takes in an Ensembl species homology tsv file and returns
        information for specific gene ID's. If a species
        of interest file is provided, the script will remove all information
        from any species not found in the lookup file.

    Input:
        1. Ensembl species holomology tsv file (required)
        2. Genes of interest lookup file (required)
        3. Species of interest lookup file (optional) 
            - Tab delimited columns of scientific_name, common_name, and enseml_ID
