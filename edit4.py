#merge_test
import os
import numpy as np
from Bio import SeqIO
import argparse
import pyvolve
import random
import glob
from collections import Counter


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
    try:
        create_NB_control(p0, p1, w0, w1, w2, tree, length, kappa, indel_rate, randomseed)
        os.system(f"{indelible} > indelible_log.txt")
        print("Indelible simulation run successfully.")
    except Exception as e:
        print(f"Error running Indelible: {e}")



####    PYVOLVE    ####

#codon list for custom matrices
def codon_list():
    """Returns a list of all codons."""
    return [
        "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
        "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
        "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
        "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"
    ]

###get char to check this
def is_synonymous(codon1, codon2):
    if len(codon1) != 3 or len(codon2) != 3:
        return False
    changes = sum(1 for a, b in zip(codon1, codon2) if a != b)
    return changes == 1

def create_matrix(omega, base_rate=0.01):
    """Creates a substitution rate matrix for codons."""
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
    """Runs Pyvolve simulations."""
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

def load_sequences(fasta_file):
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences


####    SITE RATE FUNCTIONS    ####

def get_site_rates(rate_file_path):
    """Extracts site rates from a file."""
    site_rates = []
    with open(rate_file_path, 'r') as file:
        for idx, line in enumerate(file):
            if idx < 9:  # Skip header lines
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    rate = float(parts[1])
                    site_rates.append(rate)
                except ValueError as e:
                    print(f"Skipping invalid line: {line.strip()} ({e})")
    return site_rates

def apply_site_rates(fasta_file, site_rates, output_file):
    sequences = np.array([list(seq) for seq in load_sequences(fasta_file)])
    site_rates = np.array(site_rates)
    sequence_length = len(sequences[0])
    if len(site_rates) != sequence_length:
        site_rates = np.tile(site_rates, sequence_length // len(site_rates) + 1)[:sequence_length]

    modified_sequences = []
    for idx, seq in enumerate(sequences):
        mutated_seq = np.copy(seq)
        mutation_mask = site_rates > 1.0
        mutated_seq[mutation_mask] = np.vectorize(lambda _: random.choice(['A', 'C', 'G', 'T']))(mutated_seq[mutation_mask])
        modified_sequences.append(f">taxon{idx+1}\n{''.join(mutated_seq)}\n")

    with open(output_file, "w") as f:
        f.writelines(modified_sequences)
    print(f"Modified sequences written to {output_file}")


####    MERGING    ####

def read_fasta_to_array(file_path):

    headers = []
    sequences = []
    with open(file_path, 'r') as file:
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                headers.append(line)
                if sequence:
                    sequences.append(list(sequence))
                sequence = []
            else:
                sequence.extend(line)
        if sequence:
            sequences.append(list(sequence))
    return headers, np.array(sequences, dtype="str")

def seq_together(sequences):
    max_length = max(map(len, sequences))
    padded = np.full((len(sequences), max_length), "-", dtype="str")
    
    for i, seq in enumerate(sequences):
        padded[i, :len(seq)] = list(seq)
    
    merged = []
    for col in padded.T:
        most_common = Counter(col).most_common()
        if len(most_common) > 1 and most_common[0][1] == most_common[1][1]:
            most_common_char = random.choice([char for char, _ in most_common])
        else:
            most_common_char = most_common[0][0]
        
        merged.append(most_common_char)
    
    return ''.join(merged)

def merge_fasta_with_taxon(input_files, output_file):
    all_headers = []
    all_sequences = []

    for file_path in input_files:
        headers, sequences = read_fasta_to_array(file_path)
        all_headers.extend(headers)
        all_sequences.extend(sequences)
    
    merged_sequences = []
    taxon_dict = {}
    
    for header, seq in zip(all_headers, all_sequences):
        taxon = header.split()[0]  # Assuming taxon is part of the header before a space
        if taxon not in taxon_dict:
            taxon_dict[taxon] = []
        taxon_dict[taxon].append(seq)

    for taxon, sequences in taxon_dict.items():
        merged_seq = seq_together(sequences)
        merged_sequences.append(f">{taxon}\n{merged_seq}\n")
    
    with open(output_file, 'w') as f:
        f.writelines(merged_sequences)

    print(f"Merged sequences with taxon structure saved to '{output_file}'.")



####    INDEL PATTERN   ####

def apply_indels(dna_true_file, merged_file, output_file):
    dna_true_records = list(SeqIO.parse(dna_true_file, "fasta"))
    merged_records = list(SeqIO.parse(merged_file, "fasta"))

    if len(dna_true_records) != len(merged_records):
        print("Warning: The number of taxa in dna_TRUE.fas and merged_sequences.fas do not match")
        return

    updated_records = []

    for dna_true_record, merged_record in zip(dna_true_records, merged_records):
        indel_pattern = str(dna_true_record.seq)
        binary_vector = [1 if char == '-' else 0 for char in indel_pattern]
        merged_sequence = str(merged_record.seq)

        new_sequence = []
        merged_index = 0
        for is_gap in binary_vector:
            if is_gap: 
                new_sequence.append("-")
            else:  
                if merged_index < len(merged_sequence):
                    new_sequence.append(merged_sequence[merged_index])
                    merged_index += 1

        new_sequence.extend(merged_sequence[merged_index:])
        updated_records.append(f">{merged_record.id}\n{''.join(new_sequence)}")

    with open(output_file, "w") as output_handle:
        output_handle.write("\n".join(updated_records))

    print(f"Indel patterns applied and saved to {output_file}")


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

    #INDELIBLE 
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

    #SITE RATES
    site_rates = get_site_rates(indelible_rates_file)

    #PYVOLVE
    run_pyvolve(args.tree_file, length, omega_values)

    simulation_fasta_files = glob.glob("simulation_omega_*.fas")
    for simulation_fasta in simulation_fasta_files:
        output_fasta = simulation_fasta.replace(".fas", "_with_rates.fas")
        apply_site_rates(simulation_fasta, site_rates, output_fasta)

    #MERGE

    input_files = [
        "/mnt/c/users/cerih/ebi/ebi_project/simulation_omega_1_with_rates.fas",
        "/mnt/c/users/cerih/ebi/ebi_project/simulation_omega_2_with_rates.fas",
        "/mnt/c/users/cerih/ebi/ebi_project/simulation_omega_3_with_rates.fas"
    ]
    
    output_file = "/mnt/c/users/cerih/ebi/ebi_project/merged_sequences.fas"  

    merge_fasta_with_taxon(input_files, output_file)

    #INDEL PATTERN 
    dna_true_file = "/mnt/c/users/cerih/ebi/ebi_project/dna_TRUE.fas"
    merged_file = "/mnt/c/users/cerih/ebi/ebi_project/merged_sequences.fas"
    output_file = "/mnt/c/users/cerih/ebi/ebi_project/merged_indels.fas"

    apply_indels(dna_true_file, merged_file, output_file)

    #fin

if __name__ == "__main__":
    main()
