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



#####   INDELIBLE   #####


# makes the indelible control file
def create_NB_control(p0, p1, w0, w1, w2, tree, length, kappa, indel_rate, randomseed, indel_q=0.35):
    try:
        with open("control.txt", "w") as file:
            file.write("[TYPE] CODON 1\n")
            file.write("[SETTINGS]\n[printrates] TRUE\n[randomseed] " + str(randomseed) + "\n[output] FASTA\n")
            file.write("[MODEL] modelname\n[submodel]\n" + str(kappa) + "\n" + str(p0) + " " + str(p1) + "\n" + str(w0) + " " + str(w1) + " " + str(w2))
            file.write("\n[indelrate] " + str(indel_rate) + "\n[indelmodel] NB " + str(indel_q) + " 1")
            file.write("\n[TREE] treename  \n" + tree)
            file.write("\n[PARTITIONS] partitionname [treename modelname " + str(length) + "]")
            file.write("\n[EVOLVE] partitionname 1 dna\n")
        print("Control file successfully written")
    except Exception as e:
        print(f"Error writing control file: {e}")

#run indelible
def run_indelible(indelible, randomseed, tree, length, indel_rate, kappa, p0, p1, w0, w1, w2):
    """run indelible simulation, calls create_NB_control, executes indelible , redirecting
    output to a log file indelible_log.txt"""
    try:
        create_NB_control(p0, p1, w0, w1, w2, tree, length, kappa, indel_rate, randomseed)
        os.system(f"{indelible} > indelible_log.txt")
        print("Indelible simulation run successfully.")
    except Exception as e:
        print(f"Error running Indelible: {e}")



####    PYVOLVE    ####

#codon list for custom matrices
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
        matrix[i, i] = -np.sum(matrix[i, :])  # Row sums to zero
    return matrix

#### ^ CHECK THIS 


def run_pyvolve(tree_file, sequence_length, omega_values):
    """runs evolutionary sequence simulations using pyvolve, takes in phylogenic tree, for each
    omega - creates substitution matrix and model, returns fasta file for each matrix/simualtion"""
    my_tree = pyvolve.read_tree(file=tree_file)
    codons = codon_list()

    for idx, omega in enumerate(omega_values, start=1):
        matrix = create_matrix(omega)
        custom_model = pyvolve.Model("custom", {"matrix": matrix, "code": codons})
        my_partition = pyvolve.Partition(models=custom_model, size=sequence_length)
        my_evolver = pyvolve.Evolver(partitions=my_partition, tree=my_tree)

        output_prefix = f"simulation_omega_{idx}"
        #print(f"Running simulation for omega {omega}...")   ####print statement remove
        my_evolver(seqfile=f"{output_prefix}.fas", ratefile=f"{output_prefix}_rates.txt", infofile=f"{output_prefix}_info.txt")
    print("Pyvolve simulations complete!")


###dont this this is actually used anymore - can be removed potentially
def load_sequences(dna_rates_file):
    with open(dna_rates_file, 'r') as f:
        site_rates = np.array([int(line.split()[1]) for line in f.readlines()[10:]])
    return len(site_rates)



#SITE RATES AND MERGE############################################################


from Bio import SeqIO

####this is repeated
def parse_site_rates(dna_rates_file):
    """opens the dna_RATES file and extracts the second column = rates"""
    site_rates = []
    with open(dna_rates_file, 'r') as f:
        lines = f.readlines()

        # skipping first 10 lines bc thats a header
        for line in lines[10:]:
            columns = line.split()
            if len(columns) > 1:  
                site_rates.append(int(columns[1])) 
    
    return site_rates, len(site_rates)



#merge seqeunces
# with open(fasta_files) as handle:
    #for values in SimpleFastaParser(handle):
        #print(values)
#output== (a. GATAG, b. TAGATGC, c. ACTAGCT ...)

#######################################################################
#merged sequences



def merged_sequences(dna_rates_file, fasta_files, output_file):
    with open(dna_rates_file, 'r') as f:
        site_rates = np.array([int(line.split()[1]) for line in f.readlines()[10:]])
    
    sequences_dict = [{} for _ in range(len(fasta_files))]

    for i, fasta_file in enumerate(fasta_files):
        with open(fasta_file) as f:
            for taxon_id, sequence in SimpleFastaParser(f):
                sequences_dict[i][taxon_id] = sequence
    
    taxon_ids = set(sequences_dict[0].keys())
    if not all(set(d.keys()) == taxon_ids for d in sequences_dict):
        raise ValueError("mismatch in taxon IDs across input FASTA files.")
    
    merged_sequences = {taxon_id: [] for taxon_id in taxon_ids}

    for taxon_id in taxon_ids:
        sequences = [np.array(list(sequences_dict[i][taxon_id])) for i in range(len(fasta_files))]
        site_source = site_rates #- 1 #based on indexing for the actual positions/ 0 based indexing e.g. 012345 not 12345
        merged_sequence = np.concatenate([sequences[source][idx *3:(idx + 1 ) * 3] for idx, source in enumerate(site_source)])
        merged_sequences[taxon_id] = "".join(merged_sequence)

        with open(output_file, 'w') as f:
            for taxon_id, sequence in merged_sequences.items():
                f.write(f">{taxon_id}\n{sequence}\n")
        
        print(f"merged seqs written to {output_file}")



        #############

def apply_indels(dna_true_file, merged_file, output_file, indel_rm_file):

    dna_true_dict = {}
    with open(dna_true_file) as f:
        for taxon_id, sequence in SimpleFastaParser(f):
            dna_true_dict[taxon_id] = sequence

    merged_dict = {}
    with open(merged_file) as f:
        for taxon_id, sequence in SimpleFastaParser(f):
            merged_dict[taxon_id] = sequence

    if set(dna_true_dict.keys()) != set(merged_dict.keys()):
        raise ValueError("mismatch in taxon IDs between DNA true and merged files")

    updated_sequences = {}
    indels_removed_sequences = {}

    for taxon_id in dna_true_dict:
        true_sequence = np.array(list(dna_true_dict[taxon_id]))
        merged_sequence = np.array(list(merged_dict[taxon_id]))

        gap_mask = (true_sequence == '-')
        num_non_gaps = np.sum(~gap_mask)
        print(f"Taxon: {taxon_id}, Length of true sequence: {len(true_sequence)}, "
            f"Number of non-gap bases: {num_non_gaps}, Length of merged sequence: {len(merged_sequence)}")

        if len(merged_sequence) < num_non_gaps:
            raise ValueError(
                f"Merged sequence for taxon {taxon_id} is too short. "
                f"Expected at least {num_non_gaps} bases but got {len(merged_sequence)}."
            )

        new_sequence = np.full_like(true_sequence, '-', dtype='<U1')
        new_sequence[~gap_mask] = merged_sequence[:num_non_gaps]
        updated_sequences[taxon_id] = ''.join(new_sequence)
        print(f"Taxon: {taxon_id}, Length of new {len(new_sequence)}")

        indels_removed_sequences[taxon_id] = ''.join(new_sequence[~gap_mask])


    with open(output_file, "w") as f:
        for taxon_id, sequence in updated_sequences.items():
            f.write(f">{taxon_id}\n{sequence}\n")
    print(f"Indels applied and saved to {output_file}")

    with open(indel_rm_file, "w") as f:
        for taxon_id, sequence in indels_removed_sequences.items():
            f.write(f">{taxon_id}\n{sequence}\n")
    print(f"indel-removed sequences saved to {indel_rm_file}")




##need the other output file without the indel hyphens , but indels


####    MAIN    ####
def main():
    parser = argparse.ArgumentParser(description="Run Indelible and Pyvolve simulations.")
    parser.add_argument('--tree-file', type=str, required=True, help="Path to the tree file.")
    parser.add_argument('--output-dir', type=str, required=True, help="Output directory.")
    args = parser.parse_args()

    # Parameters
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
    indelible = "/mnt/c/users/cerih/ebi/ebi_project/indelible"

    with open(args.tree_file, 'r') as tree_file:
        tree = tree_file.read().strip()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    os.chdir(args.output_dir)

#indelible

    run_indelible(indelible, randomseed, tree, length, indel_rate, kappa, p0, p1, w0, w1, w2)
    indelible_rates_file = os.path.join(args.output_dir, "dna_RATES.txt")
    if not os.path.exists(indelible_rates_file):
        print(f"Error: {indelible_rates_file} not found.")
        return


    dna_rates_file = "/mnt/c/users/cerih/ebi/ebi_project/dna_RATES.txt"
    if not os.path.exists(indelible_rates_file):
        print(f"Error: {indelible_rates_file} not found.")
        return

    

## sequence length
    sequence_length = load_sequences(dna_rates_file)
#run pyvolve
    run_pyvolve(args.tree_file, sequence_length, omega_values)

    simulation_fasta_files = glob.glob("simulation_omega_*.fas")
    for simulation_fasta in simulation_fasta_files:
        output_fasta = simulation_fasta.replace(".fas", "_with_rates.fas")



   #merge seq/site rates
    fasta_files = [
        "/mnt/c/users/cerih/ebi/ebi_project/simulation_omega_1.fas",
        "/mnt/c/users/cerih/ebi/ebi_project/simulation_omega_2.fas",
        "/mnt/c/users/cerih/ebi/ebi_project/simulation_omega_3.fas"
    ]
    output_file = "/mnt/c/users/cerih/ebi/ebi_project/merged_sequences.fas"

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    merged_sequences(dna_rates_file, fasta_files, output_file)

    print(f"Merged sequences have been written to {output_file}")

    #INDEL PATTERN 
    dna_true_file = "/mnt/c/users/cerih/ebi/ebi_project/dna_TRUE.fas"
    merged_file = "/mnt/c/users/cerih/ebi/ebi_project/merged_sequences.fas"
    output_file = "/mnt/c/users/cerih/ebi/ebi_project/merged_indels.fas"
    indel_rm_file = "/mnt/c/users/cerih/ebi/ebi_project/indels_rm.fas"

    apply_indels(dna_true_file, merged_file, output_file, indel_rm_file)

if __name__ == "__main__":
    main()
