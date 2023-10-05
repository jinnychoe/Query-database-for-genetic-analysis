Query Database for Genetic Analysis
=============================================
This repository contains a set of Python scripts designed to retrieve gene information, sequences, and species with homologous genes. This script has been tested on Linux systems and can be used to analyze gene homology.

Purpose
------------------------------------
The objective of this program is to analyze gene homology using publicly available biological databases. The pipeline involves the following steps:

1. Retrieving the Ensembl ID: The script retrieves the Ensembl ID of the gene of interest using the MyGene database.

2. Obtaining the gene sequence: Using the Ensembl ID, the script retrieves the DNA sequence of the gene from the Ensembl database.

3. Storing the sequence data: The gene sequence is saved in a FASTA file named `{gene_name}_sequence.fasta`.

4. Finding the longest ORF: The longest open reading frame (ORF) is searched for in the gene sequence.

5. Translating the ORF: The longest ORF is translated in its amino acid sequence. This amino acid sequence is then appended to the FASTA file.

6. Homologous species search: The script searches for homologous genes in other species using the MyGene database. Taxonomic IDs are returned.

7. Species name lookup: The taxonomic ID is used to determine the species name through the Entrez database.

8. Writing Homology List: The list of species with homologous genes is written to a file named `{gene_name}_homology_list.txt`.

Execution in Ubuntu/Linux:
------------------------------------
1. Ensure that the `query.py` script is located in the working directory.

2. Open a terminal and navigate to the working directory.

3. Install dependencies.

pip install biopython

pip install requests

4. Provide executable permissions to the script:

chmod +x query.py

5. Execute the script with the following command:  
 
python3 query.py
 
6. The following outputs are generated:  
 
`{gene_name}_sequence.fasta` contains the gene sequence and the translated amino acid sequence of the longest ORF.

`{gene_name}_homology_list.txt` contains the list of species names that have genes homologous to your gene of interest.
