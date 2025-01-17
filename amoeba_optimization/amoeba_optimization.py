import sys
import numpy as np

# error rate for the emission probabilities
ERROR_RATE = 1 / 1000

def parse_vcf(vcf_file):
    """parse the vcf file and return records as a list of lists."""
    with open(vcf_file, 'r') as f:
        # read each line, skipping header lines, and split by tab to get genotype information
        return [line.strip().split('\t') for line in f if not line.startswith('#')]

def emission_probabilities(genotype, p, q):
    """calculate emission probabilities for inbred and outbred states based on the genotype."""
    # check if the genotype is homozygous (either "0/0" or "1/1")
    if genotype == "0/0" or genotype == "1/1":  # homozygous
        # return probabilities for both inbred and outbred states
        return 1 - ERROR_RATE, 1 - 2 * p * q
    # check if the genotype is heterozygous (either "0/1" or "1/0")
    elif genotype == "0/1" or genotype == "1/0":  # heterozygous
        # return probabilities for both inbred and outbred states
        return ERROR_RATE, 2 * p * q
    return None, None  # return None if genotype is not recognized

def forward_algorithm(vcf_data, trans_inbred_outbred, trans_outbred_inbred):
    """compute total likelihood using the forward algorithm for hidden markov models."""
    total_likelihood = 0  # initialize the total likelihood
    # loop through each record (row) in the vcf data
    for record in vcf_data:
        _, _, _, _, *genotypes = record  # extract the genotype information
        p = 0.5  # initial probability for inbred state
        q = 1 - p  # probability for outbred state

        # initialize forward variables for inbred and outbred states
        f_inbred, f_outbred = 1.0, 1.0

        # process each genotype in the record
        for genotype in genotypes:
            # calculate emission probabilities for inbred and outbred states
            e_inbred, e_outbred = emission_probabilities(genotype, p, q)
            if e_inbred is None:  # skip if the genotype is not valid
                continue

            # calculate the next state probabilities for inbred and outbred
            next_inbred = f_inbred * (1 - trans_inbred_outbred) * e_inbred + \
                          f_outbred * trans_outbred_inbred * e_inbred
            next_outbred = f_outbred * (1 - trans_outbred_inbred) * e_outbred + \
                           f_inbred * trans_inbred_outbred * e_outbred

            # update the forward variables
            f_inbred, f_outbred = next_inbred, next_outbred

        # add the log likelihood of the final states to the total likelihood
        total_likelihood += np.log(f_inbred + f_outbred)

    return total_likelihood

def amoeba_optimization(vcf_data, initial_params, tolerance=1e-6, max_iterations=100):
    """optimize transition rates using the amoeba (nelder-mead) algorithm."""
    
    def evaluate(params):
        """evaluate the likelihood of the current parameters."""
        # negate the result of the forward algorithm to minimize the negative log-likelihood
        return -forward_algorithm(vcf_data, params[1], params[0])

    # initialize the simplex (starting points for the optimization)
    simplex = [np.array(initial_params)]
    
    # generate the perturbed points for the simplex
    for i in range(len(initial_params)):
        perturbed = initial_params[:]
        perturbed[i] *= 1.1  # increase the i-th parameter by 10%
        simplex.append(np.array(perturbed))

    # convert the simplex to a numpy array and evaluate the initial values
    simplex = np.array(simplex)
    values = [evaluate(point) for point in simplex]

    # perform the optimization using the amoeba algorithm (nelder-mead)
    for _ in range(max_iterations):
        # sort the simplex by function values
        order = np.argsort(values)
        simplex = simplex[order]
        values = [values[i] for i in order]

        # calculate the centroid (mean of the best points)
        centroid = np.mean(simplex[:-1], axis=0)

        # reflection step
        worst = simplex[-1]
        reflection = centroid + (centroid - worst)
        reflection_value = evaluate(reflection)

        # check if reflection is better than the best point
        if reflection_value < values[0]:
            # expansion step
            expansion = centroid + 2 * (reflection - centroid)
            expansion_value = evaluate(expansion)
            if expansion_value < reflection_value:
                simplex[-1] = expansion
                values[-1] = expansion_value
            else:
                simplex[-1] = reflection
                values[-1] = reflection_value
        elif reflection_value < values[-2]:
            # accept reflection if it is not worse than the second-to-worst point
            simplex[-1] = reflection
            values[-1] = reflection_value
        else:
            # contraction step
            contraction = centroid + 0.5 * (worst - centroid)
            contraction_value = evaluate(contraction)
            if contraction_value < values[-1]:
                simplex[-1] = contraction
                values[-1] = contraction_value
            else:
                # shrink the simplex towards the best point
                best = simplex[0]
                simplex = [best + 0.5 * (point - best) for point in simplex]
                values = [evaluate(point) for point in simplex]

        # check if the optimization has converged (tolerance condition)
        if np.max(np.abs(np.array(values) - values[0])) < tolerance:
            break

    return simplex[0]

if __name__ == "__main__":
    # check if the user provided the correct number of arguments (vcf file)
    if len(sys.argv) != 2:
        print("usage: python amoebaoptimization.py <input>.vcf")
        sys.exit(1)

    # get the input vcf file from command-line arguments
    vcf_file = sys.argv[1]
    # parse the vcf file to get the data
    vcf_data = parse_vcf(vcf_file)

    # define initial parameters for the transition rates
    initial_params = [1 / (4 * 10**6), 1 / (1.5 * 10**6)]

    # perform optimization of the transition rates using the amoeba algorithm
    optimized_rates = amoeba_optimization(vcf_data, initial_params)

    # print the optimized transition rates
    print(f"p(transition outbred>inbred): {optimized_rates[0]:.10f}")
    print(f"p(transition inbred>outbred): {optimized_rates[1]:.10f}")
