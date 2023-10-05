from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
import requests

# Function to get Ensembl ID for a given gene name
def get_ensembl_id(gene_name):

    # Define MyGene server URL and endpoint
    mygene_server = 'http://mygene.info/v3' 
    mygene_endpoint = '/query'
    
    # Set parameters for query to MyGene
    mygene_params = {
        'q': f'symbol:{gene_name}',
        'species': 'human',
        'fields': 'ensembl.gene',
        'size': 1
    }
    print(f'\nQuery for gene: {gene_name}')
    
    # Send request to MyGene to retrieve Ensembl ID 
    mygene_response = requests.get(mygene_server + mygene_endpoint, params=mygene_params)
    mygene_output = mygene_response.json()
    print (f'\nSending {gene_name} to Ensembl for query...')

    # Check if Ensembl ID retrieved successfully
    if mygene_output.get('hits'):
        ensembl_id = mygene_output['hits'][0]['ensembl']['gene']
        print(f'\nEnsembl ID:', ensembl_id) 
        return ensembl_id
    else:
        print(f"Ensembl ID not found for gene: {gene_name}") 
        return None

# Function to get sequence for a given Ensembl ID     
def get_sequence(ensembl_id):
    # Define Ensembl server URL and endpoint for retrieving sequence data
    ensembl_server = 'http://rest.ensembl.org'
    ensembl_endpoint = f'/sequence/id/{ensembl_id}'
    
    # Send request to Ensembl to retrieve sequence data
    seq_response = requests.get(ensembl_server + ensembl_endpoint, headers={'Content-Type': 'application/json'})
    
    # Parse response as JSON
    seq_data = seq_response.json()
    
    # Return sequence data
    print(f'\nGetting sequence from Ensembl...') 
    return seq_data['seq'] 

# Function to find longest open reading frame (ORF) in a given sequence    
def get_longest_orf(sequence):
    # Define start codon and stop codons
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    
    # Initialize variables to store longest ORF and current ORF
    longest_orf = ""
    current_orf = ""
    
    # Initialize index to traverse sequence
    i = 0
    while i < len(sequence):
    
        # Extract codon from sequence
        codon = sequence[i:i+3]
        
        # Check if codon is start codon
        if codon == start_codon:
            current_orf = codon
            i += 3
            
            # Continue to build current ORF
            while i < len(sequence):
                codon = sequence[i:i+3]
                current_orf += codon
                i += 3
                
                # Check if codon is stop codon
                if codon in stop_codons:
                
                    # Update longest ORF if current ORF is longer
                    if len(current_orf) > len(longest_orf):
                        longest_orf = current_orf
                    break
        else:
            i += 3
    
    return longest_orf 

# Function to retrieve homologous species for given gene name
def get_homologous_species(gene_name):
    # Define MyGene server URL and endpoint
    mygene_server = 'http://mygene.info/v3'
    mygene_endpoint = '/query'
    
    # Set parameters for query to MyGene to retrieve homologous species
    mygene_params = {
        'q': f'symbol:{gene_name}',
        'fields': 'homologene',
        'size': 1
    }
    
    # Send request to MyGene to search for homologous genes in other species
    mygene_response = requests.get(mygene_server + mygene_endpoint, params=mygene_params)
    mygene_output = mygene_response.json()
    print (f'\nSearching for homologous genes in other species...')

    # Check if homologous genes were found for given gene name
    if mygene_output['hits'][0]['homologene']['genes']:
    
        # Extract homologene data and create list of species
        homologene_data = mygene_output['hits'][0]['homologene']
        species_list = [entry[0] for entry in homologene_data['genes']]
        return species_list
    else:
        # Return empty list if no homologous genes were found
        return []
        
# Function to get species name using taxonomic ID        
def get_species_name(taxid):
    # Set your email for Entrez to identify user
    Entrez.email = "youremail@server.com"  
    
    # Fetch taxonomy information for given taxonomic ID in XML format
    handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    
    # Read and parse the XML records
    records = Entrez.read(handle)
    handle.close()

    # Check if taxonomy records were successfully retrieved
    if records:
        # Extract the scientific name of the species
        species_name = records[0]['ScientificName']
        return species_name
    else:
        # Return "Unknown" if no taxonomy records were found
        return "Unknown"

# Function to write species name to file        
def write_species_list_to_file(species_list, filename):
    # Open file for writing
    with open(filename, 'w') as file:
        # Iterate through the list of species IDs
        for species_id in species_list:
            # Retrieve the species name 
            species_name = get_species_name(species_id)  
            # Write the species ID and name to the file
            file.write(f'{species_id}\t{species_name}\n')
    

def main():
    # Define gene name
    gene_name = 'MC1R'
    
    # Get Ensembl ID for given gene name
    ensembl_id = get_ensembl_id(gene_name)
    
    if ensembl_id:  
        # Get DNA sequence for Ensembl ID
        sequence = get_sequence(ensembl_id)
        
        # Find longest open reading frame (ORF) in sequence
        longest_orf = get_longest_orf(sequence)
      
        # Define filename for FASTA file
        fasta_filename = f'{gene_name}_sequence.fasta'

        with open(fasta_filename, 'w') as fasta_file:
            # Write gene name and Ensembl ID as sequence header
            fasta_file.write(f'>{gene_name}_{ensembl_id}_sequence\n')
            
            # Write DNA sequence
            fasta_file.write(f'{sequence}\n')
            print(f'\n{gene_name} Ensembl sequence recorded in {gene_name}_sequence.fasta')

            if longest_orf:
                # Write the longest ORF to the FASTA file
                fasta_file.write(f'>{gene_name}_{ensembl_id}_longest_orf\n')
                fasta_file.write(f'{longest_orf}\n') 
                print(f'\nLongest ORF sequence recorded in {fasta_filename}')

                # Translate longest ORF to amino acid sequence
                protein_seq = Seq(longest_orf).translate()
                print(f'\nTranslating longest open reading frame of {gene_name} to amino acid sequence...')

                # Write protein sequence to FASTA file
                fasta_file.write(f'>{gene_name}_{ensembl_id}_protein_sequence\n')
                fasta_file.write(f'{protein_seq}\n') 
                print(f'\n{gene_name} amino acid sequence recorded in {gene_name}_sequence.fasta')

            else: 
                print("No ORF sequence found.")
    else:
        print("Ensembl ID not found for the given gene name.")

    # Get list of homologous species for gene
    homologous_species = get_homologous_species(gene_name)

    if homologous_species:
        print(f'\nHomologous genes found in other species')
        species_list = []
        for taxid in homologous_species:
            # Iterate through list of taxids and get species names
            species_name = get_species_name(taxid)
            species_list.append((taxid, species_name))

        # Define filename for species list
        species_list_filename = f'{gene_name}_homology_list.txt'
        with open(species_list_filename, 'w') as file:
            # Open file for writing and write species names
            for taxid, species_name in species_list:
                file.write(f'{species_name}\n')
        print(f'\nSpecies names recorded in {species_list_filename}\n')
    else:
        print("No homologous species found.")


if __name__ == "__main__":
    main()
