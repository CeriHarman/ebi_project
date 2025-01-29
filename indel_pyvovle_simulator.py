import os
import numpy as np
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse
import pyvolve
import random
import glob
from collections import Counter
from collections import defaultdict

# INDELIBLE FUNCTIONS

def create_NB_control(p0, p1, w0, w1, w2, tree, length, kappa, indel_rate, randomseed, indel_q=0.35):
    try:
        with open(os.path.join(os.getcwd(), "control.txt"), "w") as file:
            file.write("[TYPE] CODON 1\n")
            file.write("[SETTINGS]\n[printrates] TRUE\n[randomseed] " + str(randomseed) + "\n[output] FASTA\n")
            file.write("[MODEL] modelname\n[submodel]\n" + str(kappa) + "\n" + str(p0) + " " + str(p1) + "\n" + str(w0) + " " + str(w1) + " " + str(w2))
            file.write("\n[indelrate] " + str(indel_rate) + "\n[indelmodel] NB " + str(indel_q) + " 1")
            file.write("\n[TREE] treename  \n" + tree)
            file.write("\n[PARTITIONS] partitionname [treename modelname " + str(length) + "]")
            file.write("\n[EVOLVE] partitionname 1 dna\n")
       # print("Control file successfully written")
    except Exception as e:
        print(f"Error writing control file: {e}")


def run_indelible(indelible, randomseed, tree, length, indel_rate, kappa, p0, p1, w0, w1, w2):
    try:
        create_NB_control(p0, p1, w0, w1, w2, tree, length, kappa, indel_rate, randomseed)
        os.system(f"{os.path.join(os.getcwd(), indelible)} > indelible_log.txt")
        #print("Indelible simulation run successfully.")
    except Exception as e:
        print(f"Error running Indelible: {e}")



### USED TO GET SEQ LENGTH FROM INDELIBLE FOR PYVOLVE
def load_sequences(dna_rates_file): #load the dna rates file - because it contains the number of codon sites (pyvovle also works with codon sites) and counts the length = number of sites
    with open(dna_rates_file, 'r') as f: 
        site_rates = np.array([int(line.split()[1]) for line in f.readlines()[10:]]) # first 10 lines are a header/ so remove this and access second column = where the rates are in the file
    return len(site_rates)

# PYVOVLE FUNCTIONS

#list of codons for matrices to use
def codon_list(): 
    """provides a list of 61 codons minus stop codons - for the custom matrix"""
    return [
        "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
        "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
        "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
        "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"
    ]

###get char to check this
def is_synonymous(codon1, codon2):
    """check if two codons are synonymous, and returns TRUE if a nucleotide
    differs between them"""
    if len(codon1) != 3 or len(codon2) != 3:
        return False
    changes = sum(1 for a, b in zip(codon1, codon2) if a != b)
    return changes == 1

def create_matrix(omega, base_rate=0.01):
    """created substitution rate matrix with base subsitution rates, rates of 
    non/synonymous can be adjusted with omega"""
    codons = codon_list()
    matrix = np.zeros((61, 61))
    for i in range(61):
        for j in range(61):
            if i != j:
                if is_synonymous(codons[i], codons[j]):
                    matrix[i, j] = base_rate
                else:
                    matrix[i, j] = base_rate * omega
                if j > i:
                    matrix[i,j] = matrix[i,j]*2
        matrix[i, i] = -np.sum(matrix[i, :])  # Row sums to zero
    return matrix

#### ^ CHECK THIS 

#run pyvolve, takes original use tree sequence, omega values are set in main = need to change, and sequence length comes from load_sequences / site rates ffrom indelible
def run_pyvolve(tree_file, sequence_length, omega_values):
    """runs evolutionary sequence simulations using pyvolve, takes in phylogenic tree, for each
    omega - creates substitution matrix and model, returns fasta file for each matrix/simualtion"""
    my_tree = pyvolve.read_tree(file=tree_file)
    codons = codon_list()

    for idx, omega in enumerate(omega_values, start=1):
        matrix = create_matrix(omega) #custom matrices using omega values inputted in MAIN ()
        custom_model = pyvolve.Model("custom", {"matrix": matrix, "code": codons})  #makes the matrices into an actual model to use 
        my_partition = pyvolve.Partition(models=custom_model, size=sequence_length) # creates a partition - at least one is needed to 'evolve' pyvolve
        my_evolver = pyvolve.Evolver(partitions=my_partition, tree=my_tree) # takes the partition e.g. what model and length from above and tree - to run pyvolve

        output_prefix = f"simulation_omega_{idx}"
        my_evolver(seqfile=f"{output_prefix}.fas", ratefile=f"{output_prefix}_rates.txt", infofile=f"{output_prefix}_info.txt") # line that actually runs pyvolve and gives the different out puts e.g. simulation_omega_1/2/3.fas
    #print("Pyvolve simulations complete!")





# MERGE PYVOLVE SIMULATIONS BASED ON SITE RATES FROM INDELIBLE

#takes the three simualtions from pyvolve and merged them based on the site rates from indelible, producing a fasta with three different omegas/site rates
def merged_sequences(dna_rates_file, fasta_files, output_file):
    with open(dna_rates_file, 'r') as f:  # might be able to remove this and just use load_sequences - so an make more streamline
        site_rates = np.array([int(line.split()[1]) for line in f.readlines()[10:]]) #opens file and parses header, looks at column 2 = site rates e.g. 12012102100011222
    
    sequences_dict = [{} for _ in range(len(fasta_files))]  

    for i, fasta_file in enumerate(fasta_files): # open the pyvovle simulated data
        with open(fasta_file) as f:
            for taxon_id, sequence in SimpleFastaParser(f): # separate the taxons and sequences into a dictionary , to keep the structure
                sequences_dict[i][taxon_id] = sequence
    
    taxon_ids = set(sequences_dict[0].keys())
    if not all(set(d.keys()) == taxon_ids for d in sequences_dict):
        raise ValueError("mismatch in taxon IDs across input FASTA files.") # make sure they are the same 
    
    merged_sequences = {taxon_id: [] for taxon_id in taxon_ids} 

    for taxon_id in taxon_ids: # for each taxon in the dictionary - look ath the sequences and merge - 
        sequences = [np.array(list(sequences_dict[i][taxon_id])) for i in range(len(fasta_files))]
        #site_source = site_rates 
        merged_sequence = np.concatenate([sequences[source][idx *3:(idx + 1 ) * 3] for idx, source in enumerate(site_rates)]) # merging based on site_source = site rates, just changes #site_source to site_rates
        merged_sequences[taxon_id] = "".join(merged_sequence) #join

        with open(output_file, 'w') as f:
            for taxon_id, sequence in merged_sequences.items():
                f.write(f">{taxon_id}\n{sequence}\n")    #save as new file to 'merged_sequences"
        
        #print(f"merged seqs written to {output_file}")



        #############
# dna true from indelible - to get the pattern from, merged_file from function above - file we are applying indels to, output file = output with - as indels, indel_rm = indels are there but removed
def apply_indels(dna_true_file, merged_file, output_file, indel_rm_file):

    dna_true_dict = {}  #  start dictionary from true file - inorder to strip for pattern, separate taxons and sequences
    with open(dna_true_file) as f:
        for taxon_id, sequence in SimpleFastaParser(f):
            dna_true_dict[taxon_id] = sequence

    merged_dict = {}  #start dictionary for merged - separate taxond and sequences
    with open(merged_file) as f:
        for taxon_id, sequence in SimpleFastaParser(f):
            merged_dict[taxon_id] = sequence

    if set(dna_true_dict.keys()) != set(merged_dict.keys()):   #make sure the same length
        raise ValueError("mismatch in taxon IDs between DNA true and merged files")

    updated_sequences = {}   #new sequence with indels
    indels_removed_sequences = {} # seq with indels but removed

    for taxon_id in dna_true_dict:
        true_sequence = np.array(list(dna_true_dict[taxon_id]))  
        merged_sequence = np.array(list(merged_dict[taxon_id]))

        gap_mask = (true_sequence == '-')
        num_non_gaps = np.sum(~gap_mask)
        #print(f"Taxon: {taxon_id}, Length of true sequence: {len(true_sequence)}, "
        #    f"Number of non-gap bases: {num_non_gaps}, Length of merged sequence: {len(merged_sequence)}")

        if len(merged_sequence) < num_non_gaps:
            raise ValueError(
                f"Merged sequence for taxon {taxon_id} is too short. "
                f"Expected at least {num_non_gaps} bases but got {len(merged_sequence)}."
            )

        new_sequence = np.full_like(true_sequence, '-', dtype='<U1')
        new_sequence[~gap_mask] = merged_sequence[:num_non_gaps]
        updated_sequences[taxon_id] = ''.join(new_sequence)
       # print(f"Taxon: {taxon_id}, Length of new {len(new_sequence)}")

        indels_removed_sequences[taxon_id] = ''.join(new_sequence[~gap_mask])


    with open(output_file, "w") as f:
        for taxon_id, sequence in updated_sequences.items():
            f.write(f">{taxon_id}\n{sequence}\n")
    print(f"files with indels as '-' saved to {output_file}")

    with open(indel_rm_file, "w") as f:
        for taxon_id, sequence in indels_removed_sequences.items():
            f.write(f">{taxon_id}\n{sequence}\n")
    print(f"files with indels but no '-' saved to {indel_rm_file}")



def main():
    parser = argparse.ArgumentParser(description="Run Indelible and Pyvolve simulations.")
    parser.add_argument('--tree-file', type=str, required=True, help="Path to the tree file.")
    parser.add_argument('--output-dir', type=str, required=True, help="Output directory.")
    args = parser.parse_args()

    randomseed = 12345
    length = 500
    indel_rate = 0.1
    kappa = 2.0
    p0 = 0.3
    p1 = 0.3
    w0 = 0.1
    w1 = 0.8
    w2 = 2.5
    omega_values = [0.2, 1.0, 2.0]
    indelible = "indelible"

    with open(args.tree_file, 'r') as tree_file:
        tree = tree_file.read().strip()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    os.chdir(args.output_dir)

    run_indelible(indelible, randomseed, tree, length, indel_rate, kappa, p0, p1, w0, w1, w2)
    indelible_rates_file = os.path.join(args.output_dir, "dna_RATES.txt")
    if not os.path.exists(indelible_rates_file):
        print(f"Error: {indelible_rates_file} not found.")
        return

    dna_rates_file = os.path.join(args.output_dir, "dna_RATES.txt")
    if not os.path.exists(indelible_rates_file):
        print(f"Error: {indelible_rates_file} not found.")
        return

    sequence_length = load_sequences(dna_rates_file)
    run_pyvolve(args.tree_file, sequence_length, omega_values)

    simulation_fasta_files = glob.glob("simulation_omega_*.fas")
    for simulation_fasta in simulation_fasta_files:
        output_fasta = simulation_fasta.replace(".fas", "_with_rates.fas")

    fasta_files = [os.path.join(args.output_dir, f"simulation_omega_{i}.fas") for i in range(1, 4)]
    output_file = os.path.join(args.output_dir, "merged_sequences.fas")

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    merged_sequences(dna_rates_file, fasta_files, output_file)
   # print(f"Merged sequences have been written to {output_file}")

    dna_true_file = os.path.join(args.output_dir, "dna_TRUE.fas")
    merged_file = os.path.join(args.output_dir, "merged_sequences.fas")
    output_file = os.path.join(args.output_dir, "merged_indels.fas")
    indel_rm_file = os.path.join(args.output_dir, "indels_rm.fas")

    apply_indels(dna_true_file, merged_file, output_file, indel_rm_file)

    os.system("rm merged_sequences.fas")

    # Use glob to find and remove all matching files
    for file in glob.glob("simulation_omega_*"):
        os.system(f"rm {file}")

if __name__ == "__main__":
    main()
