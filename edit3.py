
#the main file that is all the 4 other sections, #indelible task # pyvolve_task # vectorising # merging_data 
# aim to be a simple .py file that can run with any tree file as input and will return an organised fasta file

#indelible_task
import os
import numpy as np
from Bio import SeqIO
import argparse
import pyvolve
import random
import glob


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
        print(f"Control file successfully written")
    except Exception as e:
        print(f"Error writing control file: {e}")


def run_indelible(indelible, randomseed, tree, length, indel_rate, kappa, p0, p1, w0, w1, w2):
    try:
        create_NB_control(p0, p1, w0, w1, w2, tree, length, kappa, indel_rate, randomseed)
        os.system(f"{indelible} > indelible_log.txt")
        print("Indelible simulation run successfully.")
    except Exception as e:
        print(f"Error running Indelible: {e}")


def seq_length(fasta_file):
    """
    Count the total length of sequences in a FASTA file.
    Args:
        fasta_file (str): Path to the FASTA file.
    Returns:
        int: Total sequence length (sum of all sequence lengths).
    """
    total_length = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        total_length += len(record.seq)
    return total_length
print("total length{total_length}")

#for the code to make a custom matrix
def codon_list():
    """pre-defined list of 61 codons (excluding stop codons)."""
    return [
        "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
        "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
        "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
        "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"
    ]


#for omega>
def is_synonymous(codon1, codon2):
    """Determine if two codons are synonymous."""
    if len(codon1) != 3 or len(codon2) != 3:
        return False
    changes = sum(1 for a, b in zip(codon1, codon2) if a != b)
    return changes == 1

#create the matrices based on omega input
def create_matrix(omega, base_rate=0.01):
    """create a matrix for a given omega."""
    codons = codon_list()
    matrix = np.zeros((61, 61))
    for i in range(61):
        for j in range(61):
            if i != j:
                if is_synonymous(codons[i], codons[j]):
                    matrix[i, j] = base_rate
                else:
                    matrix[i, j] = base_rate * omega
        matrix[i, i] = -np.sum(matrix[i, :])  # row sums to zero
    return matrix


#####
def load_sequences(fasta_file):
    """Load sequences from a FASTA file into a list of strings."""
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences



def get_site_rates(rate_file_path):
    site_rates = []  # Initialize an empty list to store site rates

    with open(rate_file_path, 'r') as file:
        for idx, line in enumerate(file):
            # Skip the first 9 lines (header and unnecessary rows)
            if idx < 9:
                continue

            # Split the line into parts
            parts = line.strip().split()

            # Skip the header line if it exists on line 10
            if idx == 9 and parts == ["Site", "Class", "Partition", "Inserted?"]:
                continue

            # Ensure the line has enough parts to extract the rate
            if len(parts) >= 2:
                try:
                    # Convert the second column to an integer and append to the list
                    rate = float(parts[1])  # Use float to handle any decimal values
                    site_rates.append(rate)
                except ValueError as e:
                    # Handle invalid rate values gracefully
                    print(f"Skipping invalid line: {line.strip()} ({e})")

    return site_rates

def apply_site_rates(fasta_file, site_rates, output_file):
    # Load sequences into a list of strings
    sequences = load_sequences(fasta_file)
    sequences = np.array([list(seq) for seq in sequences])  # Convert to 2D numpy array

    # Convert site rates to numpy array
    site_rates = np.array(site_rates)

    # Ensure site_rates and sequence lengths match
    sequence_length = len(sequences[0])
    if len(site_rates) != sequence_length:
        # If site rates are fewer, repeat them to match sequence length
        repetitions = sequence_length // len(site_rates) + 1  # +1 to ensure we cover the full length
        site_rates = np.tile(site_rates, repetitions)[:sequence_length]  # Slice to match exact sequence length

    # Apply site rates to each sequence
    modified_sequences = []

    for idx, seq in enumerate(sequences):
        mutated_seq = np.copy(seq)  # Copy the sequence to mutate

        # create a mutation mask: positions with site rate > 1.0 should be mutated
        mutation_mask = site_rates > 1.0  # Use site rates to decide where mutations should occur

        # apply mutations based on the mutation mask
        mutated_seq[mutation_mask] = np.vectorize(lambda _: random.choice(['A', 'C', 'G', 'T']))(mutated_seq[mutation_mask])

        # ad the modified sequence to the list with the correct header format
        modified_sequences.append(f">taxon{idx+1}\n{''.join(mutated_seq)}\n")

    with open(output_file, "w") as f:
        f.writelines(modified_sequences)
    print(f"Modified sequences written to {output_file}")



###
#running pyvolve but taking the seq length from the function/ so that it can use the number of sites from indelible but /3 so that it can accommodate for the x3 in indelible
def run_pyvolve(tree_file, sequence_length, omega_values):
    """run Pyvolve simulations for multiple omega values."""
    my_tree = pyvolve.read_tree(file=tree_file)
    codons = codon_list()

    for idx, omega in enumerate(omega_values, start=1):
        matrix = create_matrix(omega)
        custom_model = pyvolve.Model("custom", {"matrix": matrix, "code": codons})
        my_partition = pyvolve.Partition(models=custom_model, size=sequence_length)
        my_evolver = pyvolve.Evolver(partitions=my_partition, tree=my_tree)

        output_prefix = f"simulation_omega_{idx}"
        print(f"Running simulation for omega {omega}...")
        my_evolver(seqfile=f"{output_prefix}.fas", ratefile=f"{output_prefix}_rates.txt", infofile=f"{output_prefix}_info.txt")
    print("Pyvolve simulations complete!")



def main():
    parser = argparse.ArgumentParser(description="Run Indelible and Pyvolve simulations")
    parser.add_argument('--tree-file', type=str, required=True, help="Path to the tree file you want to use")
    parser.add_argument('--output-dir', type=str, required=True, help="Directory where you want the output files to go")
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
    omega_values = [0.2, 1.0, 2.0]  # List of omega values for simulation
    indelible = "/mnt/c/users/cerih/ebi/ebi_project/indelible"

    with open(args.tree_file, 'r') as tree_file:
        tree = tree_file.read().strip()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    os.chdir(args.output_dir)

    # 1 - Run Indelible
    run_indelible(indelible, randomseed, tree, length, indel_rate, kappa, p0, p1, w0, w1, w2)
    print("Indelible simulation complete.")

    indelible_rates_file = os.path.join(args.output_dir, "dna_RATES.txt")

    if not os.path.exists(indelible_rates_file):
        print(f"Error: {indelible_rates_file} not found in the specified output directory.")
        return

    # 2 - Get site rates from Indelible output
    site_rates = get_site_rates(indelible_rates_file)
    print(f"Site rates extracted: {site_rates}")

    # 3 - Run Pyvolve simulations
    run_pyvolve(args.tree_file, length, omega_values)

    # 4 - Apply site rates to Pyvolve-generated simulation files
    simulation_fasta_files = glob.glob("simulation_omega_*.fas")

    if not simulation_fasta_files:
        print("Error: No FASTA simulation files found in the output directory.")
        return

    for simulation_fasta in simulation_fasta_files:
        output_fasta = simulation_fasta.replace(".fas", "_with_rates.fas")
        apply_site_rates(simulation_fasta, site_rates, output_fasta)

    print("Site rates applied successfully to all simulations.")

if __name__ == "__main__":
    main()
