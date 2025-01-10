import sys
import numpy as np

P_INBRED_TO_OUTBRED = 1 / (1.5 * 10**6)
P_OUTBRED_TO_INBRED = 1 / (4 * 10**6)

e = 1 / 1000
EPSILON = 1e-10  

def emission_prob(genotype, p, q, inbred_state):
    if genotype == '0/0' or genotype == '1/1':  # Homozygous
        if inbred_state == 'inbred':
            return max(1 - e, EPSILON)
        else:  # outbred
            return max(1 - 2 * p * q, EPSILON)
    elif genotype == '0/1' or genotype == '1/0':  # Heterozygous
        if inbred_state == 'inbred':
            return max(e, EPSILON)
        else:  # outbred
            return max(2 * p * q, EPSILON)
    return EPSILON

def viterbi(genotypes, allele_freqs):
    n = len(genotypes)
    states = ['inbred', 'outbred']
    
    v = np.full((2, n), -np.inf) 
    backpointer = np.zeros((2, n), dtype=int) 

    p, q = allele_freqs[0], 1 - allele_freqs[0]
    v[0, 0] = np.log(emission_prob(genotypes[0], p, q, 'inbred')) + np.log(0.5)
    v[1, 0] = np.log(emission_prob(genotypes[0], p, q, 'outbred')) + np.log(0.5)

    for t in range(1, n):
        p, q = allele_freqs[t], 1 - allele_freqs[t]
        for s in range(2):
            if s == 0:  # Inbred state
                trans_probs = [v[0, t-1] + np.log(1 - P_INBRED_TO_OUTBRED), v[1, t-1] + np.log(P_OUTBRED_TO_INBRED)]
                v[0, t] = np.log(emission_prob(genotypes[t], p, q, 'inbred')) + max(trans_probs)
                backpointer[0, t] = np.argmax(trans_probs)
            else:  # Outbred state
                trans_probs = [v[0, t-1] + np.log(P_INBRED_TO_OUTBRED), v[1, t-1] + np.log(1 - P_OUTBRED_TO_INBRED)]
                v[1, t] = np.log(emission_prob(genotypes[t], p, q, 'outbred')) + max(trans_probs)
                backpointer[1, t] = np.argmax(trans_probs)

    best_path = []
    last_state = np.argmax(v[:, -1])
    best_path.append(states[last_state])

    for t in range(n - 1, 0, -1):
        last_state = backpointer[last_state, t]
        best_path.append(states[last_state])

    best_path.reverse()
    return best_path

def parse_vcf(file_path):
    """Parse VCF file and extract genotypes for each individual by chromosome."""
    genotypes_by_individual = {}
    
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("#CHROM"):
                headers = line.strip().split("\t")
                individuals = headers[9:]  
                for individual in individuals:
                    genotypes_by_individual[individual] = {}
                continue
            if line.startswith("#"):
                continue  
            
            fields = line.strip().split("\t")
            chrom, pos, _, ref, alt, _, _, info, _, *sample_genotypes = fields
            pos = int(pos)

            info_dict = {item.split("=")[0]: item.split("=")[1] for item in info.split(";") if "=" in item}
            p = float(info_dict.get("AF", 0.5))  

            for i, genotype_info in enumerate(sample_genotypes):
                genotype = genotype_info.split(":")[0].replace("|", "/")  
                individual = individuals[i]
                
                if chrom not in genotypes_by_individual[individual]:
                    genotypes_by_individual[individual][chrom] = {"positions": [], "genotypes": [], "p_values": []}
                
                genotypes_by_individual[individual][chrom]["positions"].append(pos)
                genotypes_by_individual[individual][chrom]["genotypes"].append(genotype)
                genotypes_by_individual[individual][chrom]["p_values"].append(p)

    return genotypes_by_individual

def main(vcf_file):
    genotypes_by_individual = parse_vcf(vcf_file)
    print("individual\tstart\tstop")

    for individual, chromosomes in genotypes_by_individual.items():
        for chrom, data in chromosomes.items():
            positions = data["positions"]
            genotypes = data["genotypes"]
            p_values = data["p_values"]
            
            best_path = viterbi(genotypes, p_values)

            # Extract inbred regions
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

            # Merge consecutive inbred regions if they are adjacent or close
            merged_regions = []
            for start, stop in inbred_regions:
                if merged_regions and start <= merged_regions[-1][1] + 1:
                    merged_regions[-1] = (merged_regions[-1][0], stop)
                else:
                    merged_regions.append((start, stop))

            for start, stop in merged_regions:
                print(f"{individual}\t{start}\t{stop}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python submission.py input.vcf")
        sys.exit(1)

    vcf_file = sys.argv[1]
    main(vcf_file)
