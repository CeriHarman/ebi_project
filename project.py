import os
import numpy as np
from Bio import SeqIO
import argparse
import pyvolve


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
    omega_values = [0.2, 1.0, 2.0] ###omegaaaa values - need to intergrate this into the main matrix
    indelible = "/mnt/c/users/cerih/ebi/ebi_project/indelible"

    with open(args.tree_file, 'r') as tree_file:
        tree = tree_file.read().strip()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    os.chdir(args.output_dir)

    # 1 -  run Indelible
    run_indelible(indelible, randomseed, tree, length, indel_rate, kappa, p0, p1, w0, w1, w2)
    print("Indelible simulation complete.")

    # 2 - get sequence length from Indelible output
    fasta_file = "dna.fas"  # Adjust if Indelible uses a different output filename
    sequence_length = seq_length(fasta_file)
    print(f"Total sequence length from Indelible: {sequence_length}")

    # 4 - run Pyvolve simulations
    run_pyvolve(args.tree_file, sequence_length, omega_values)

    #5 add site rates
    # 6 merge indeible simulations
    #7 add indels
    #8 make copy and remove ---



if __name__ == "__main__":
    main()
