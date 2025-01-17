import sys
import numpy as np
import matplotlib.pyplot as plt

# transition probabilities
P_INBRED_TO_OUTBRED = 1 / (1.5 * 10**6)
P_OUTBRED_TO_INBRED = 1 / (4 * 10**6)

# error rate and epsilon for numerical stability
e = 1 / 1000
EPSILON = 1e-10  

def emission_prob(genotype, p, q, inbred_state):
    """calculate emission probability based on genotype and state."""
    if genotype == '0/0' or genotype == '1/1':  # homozygous
        if inbred_state == 'inbred':
            return max(1 - e, EPSILON)  # inbred state probability
        else:  # outbred state
            return max(1 - 2 * p * q, EPSILON)  # outbred state probability
    elif genotype == '0/1' or genotype == '1/0':  # heterozygous
        if inbred_state == 'inbred':
            return max(e, EPSILON)  # inbred state probability
        else:  # outbred state
            return max(2 * p * q, EPSILON)  # outbred state probability
    return EPSILON  # return a small value for undefined genotypes

def viterbi(genotypes, allele_freqs):
    """perform viterbi algorithm to determine the most probable states."""
    n = len(genotypes)  # number of genotypes
    states = ['inbred', 'outbred']  # possible states
    
    # initialize dynamic programming table and backpointer table
    v = np.full((2, n), -np.inf) 
    backpointer = np.zeros((2, n), dtype=int) 

    p, q = allele_freqs[0], 1 - allele_freqs[0]  # allele frequencies
    # initialize the first column of the viterbi table
    v[0, 0] = np.log(emission_prob(genotypes[0], p, q, 'inbred')) + np.log(0.5)
    v[1, 0] = np.log(emission_prob(genotypes[0], p, q, 'outbred')) + np.log(0.5)

    # fill the viterbi table for each genotype
    for t in range(1, n):
        p, q = allele_freqs[t], 1 - allele_freqs[t]  # allele frequencies for the current position
        for s in range(2):
            if s == 0:  # inbred state
                trans_probs = [v[0, t-1] + np.log(1 - P_INBRED_TO_OUTBRED), v[1, t-1] + np.log(P_OUTBRED_TO_INBRED)]
                v[0, t] = np.log(emission_prob(genotypes[t], p, q, 'inbred')) + max(trans_probs)
                backpointer[0, t] = np.argmax(trans_probs)
            else:  # outbred state
                trans_probs = [v[0, t-1] + np.log(P_INBRED_TO_OUTBRED), v[1, t-1] + np.log(1 - P_OUTBRED_TO_INBRED)]
                v[1, t] = np.log(emission_prob(genotypes[t], p, q, 'outbred')) + max(trans_probs)
                backpointer[1, t] = np.argmax(trans_probs)

    # traceback to find the best path
    best_path = []
    last_state = np.argmax(v[:, -1])
    best_path.append(states[last_state])

    for t in range(n - 1, 0, -1):
        last_state = backpointer[last_state, t]
        best_path.append(states[last_state])

    best_path.reverse()  # reverse the path since we traced it backwards
    return best_path

def parse_vcf(file_path):
    """parse vcf file and extract genotypes for each individual by chromosome."""
    genotypes_by_individual = {}
    
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("#CHROM"):
                headers = line.strip().split("\t")
                individuals = headers[9:]  # get the list of individuals
                for individual in individuals:
                    genotypes_by_individual[individual] = {}
                continue
            if line.startswith("#"):
                continue  # skip comment lines
            
            fields = line.strip().split("\t")
            chrom, pos, _, ref, alt, _, _, info, _, *sample_genotypes = fields
            pos = int(pos)

            # extract allele frequency from the INFO field
            info_dict = {item.split("=")[0]: item.split("=")[1] for item in info.split(";") if "=" in item}
            p = float(info_dict.get("AF", 0.5))  # default allele frequency is 0.5 if not provided

            for i, genotype_info in enumerate(sample_genotypes):
                genotype = genotype_info.split(":")[0].replace("|", "/")  # extract genotype
                individual = individuals[i]
                
                if chrom not in genotypes_by_individual[individual]:
                    genotypes_by_individual[individual][chrom] = {"positions": [], "genotypes": [], "p_values": []}
                
                # store the genotype and associated position and allele frequency
                genotypes_by_individual[individual][chrom]["positions"].append(pos)
                genotypes_by_individual[individual][chrom]["genotypes"].append(genotype)
                genotypes_by_individual[individual][chrom]["p_values"].append(p)

    return genotypes_by_individual

def plot_inbred_outbred(genotypes_by_individual):
    """plot inbred and outbred regions for each individual."""
    # create a figure for the plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # iterate through each individual and plot their inbred/outbred regions
    y_pos = 0  # initial y position for plotting each individual
    for individual, chromosomes in genotypes_by_individual.items():
        for chrom, data in chromosomes.items():
            positions = data["positions"]
            genotypes = data["genotypes"]
            p_values = data["p_values"]

            # run viterbi to get the best path (inbred/outbred states)
            best_path = viterbi(genotypes, p_values)

            # plot each region for the current individual
            for i in range(1, len(best_path)):
                if best_path[i - 1] == 'inbred':
                    ax.plot([positions[i - 1], positions[i]], [y_pos, y_pos], color='blue', lw=6)
                else:
                    ax.plot([positions[i - 1], positions[i]], [y_pos, y_pos], color='red', lw=6)

            y_pos += 1  # move to the next individual

    # add labels and title
    ax.set_xlabel('Genomic Position')
    ax.set_ylabel('Individuals')
    ax.set_title('Inbred vs Outbred Regions Across Genomic Positions')
    plt.tight_layout()
    plt.show()

def main(vcf_file):
    """main function to parse vcf file and plot inbred vs outbred regions."""
    genotypes_by_individual = parse_vcf(vcf_file)
    print("individual\tstart\tstop")

    for individual, chromosomes in genotypes_by_individual.items():
        for chrom, data in chromosomes.items():
            positions = data["positions"]
            genotypes = data["genotypes"]
            p_values = data["p_values"]
            
            # run viterbi algorithm to get the best path for this individual
            best_path = viterbi(genotypes, p_values)

            # extract inbred regions
            inbred_regions = []
            current_start = None
            for i, state in enumerate(best_path):
                if state == 'inbred':
                    if current_start is None:
                        current_start = positions[i]
                else:
                    if current_start is not None:
                        inbred_regions.append((current_start, positions[i - 1]))
                        current_start = None
            if current_start is not None:
                inbred_regions.append((current_start, positions[-1]))

            # merge consecutive inbred regions if they are adjacent or close
            merged_regions = []
            for start, stop in inbred_regions:
                if merged_regions and start <= merged_regions[-1][1] + 1:
                    merged_regions[-1] = (merged_regions[-1][0], stop)
                else:
                    merged_regions.append((start, stop))

            for start, stop in merged_regions:
                print(f"{individual}\t{start}\t{stop}")

    # plot the inbred/outbred regions after processing
    plot_inbred_outbred(genotypes_by_individual)


if __name__ == "__main__":
    # check if the user provided the correct number of arguments (vcf file)
    if len(sys.argv) != 2:
        print("usage: python viterbi_algorithm.py input.vcf")
        sys.exit(1)

    # get the input vcf file from command-line arguments
    vcf_file = sys.argv[1]
    # run the main function
    main(vcf_file)
