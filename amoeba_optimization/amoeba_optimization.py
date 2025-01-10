import sys
import numpy as np

ERROR_RATE = 1 / 1000

def parse_vcf(vcf_file):
    """Parse the VCF file and return records."""
    with open(vcf_file, 'r') as f:
        return [line.strip().split('\t') for line in f if not line.startswith('#')]

def emission_probabilities(genotype, p, q):
    """Calculate emission probabilities for inbred and outbred states."""
    if genotype == "0/0" or genotype == "1/1":  # Homozygous
        return 1 - ERROR_RATE, 1 - 2 * p * q
    elif genotype == "0/1" or genotype == "1/0":  # Heterozygous
        return ERROR_RATE, 2 * p * q
    return None, None 

def forward_algorithm(vcf_data, trans_inbred_outbred, trans_outbred_inbred):
    """Compute total likelihood using the Forward Algorithm."""
    total_likelihood = 0
    for record in vcf_data:
        _, _, _, _, *genotypes = record
        p = 0.5  
        q = 1 - p

        f_inbred, f_outbred = 1.0, 1.0 

        for genotype in genotypes:
            e_inbred, e_outbred = emission_probabilities(genotype, p, q)
            if e_inbred is None: 
                continue

            next_inbred = f_inbred * (1 - trans_inbred_outbred) * e_inbred + \
                          f_outbred * trans_outbred_inbred * e_inbred
            next_outbred = f_outbred * (1 - trans_outbred_inbred) * e_outbred + \
                           f_inbred * trans_inbred_outbred * e_outbred

            f_inbred, f_outbred = next_inbred, next_outbred

        total_likelihood += np.log(f_inbred + f_outbred) 

    return total_likelihood

def amoeba_optimization(vcf_data, initial_params, tolerance=1e-6, max_iterations=100):
    """Optimize transition rates using the Amoeba algorithm."""
    def evaluate(params):
        """Evaluate likelihood with given transition rates."""
        return -forward_algorithm(vcf_data, params[1], params[0])

    simplex = [np.array(initial_params)]
    for i in range(len(initial_params)):
        perturbed = initial_params[:]
        perturbed[i] *= 1.1 
        simplex.append(np.array(perturbed))

    simplex = np.array(simplex)
    values = [evaluate(point) for point in simplex]

    for _ in range(max_iterations):
        order = np.argsort(values)
        simplex = simplex[order]
        values = [values[i] for i in order]

        centroid = np.mean(simplex[:-1], axis=0)

        worst = simplex[-1]
        reflection = centroid + (centroid - worst)
        reflection_value = evaluate(reflection)

        if reflection_value < values[0]:  
            expansion = centroid + 2 * (reflection - centroid)
            expansion_value = evaluate(expansion)
            if expansion_value < reflection_value:
                simplex[-1] = expansion
                values[-1] = expansion_value
            else:
                simplex[-1] = reflection
                values[-1] = reflection_value
        elif reflection_value < values[-2]:  
            simplex[-1] = reflection
            values[-1] = reflection_value
        else:  
            contraction = centroid + 0.5 * (worst - centroid)
            contraction_value = evaluate(contraction)
            if contraction_value < values[-1]:
                simplex[-1] = contraction
                values[-1] = contraction_value
            else: 
                best = simplex[0]
                simplex = [best + 0.5 * (point - best) for point in simplex]
                values = [evaluate(point) for point in simplex]

        if np.max(np.abs(np.array(values) - values[0])) < tolerance:
            break

    return simplex[0]

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python FirstName_LastName_Tier2.py input.vcf")
        sys.exit(1)

    vcf_file = sys.argv[1]
    vcf_data = parse_vcf(vcf_file)

    initial_params = [1 / (4 * 10**6), 1 / (1.5 * 10**6)]

    optimized_rates = amoeba_optimization(vcf_data, initial_params)

    print(f"P(transition outbred>inbred): {optimized_rates[0]:.10f}")
    print(f"P(transition inbred>outbred): {optimized_rates[1]:.10f}")
