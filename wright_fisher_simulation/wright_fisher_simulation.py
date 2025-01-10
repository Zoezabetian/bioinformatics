import argparse
import numpy as np

def simulate_wright_fisher(allele_freq, pop_size, fitness):
    """Simulate Wright-Fisher dynamics with selection in a haploid population until fixation or loss."""
    freq = allele_freq
    generations = 0
    
    while 0 < freq < 1:
        # calculate effective fitness-adjusted frequency
        p_selected = freq * fitness / (freq * fitness + (1 - freq))
        
        # update allele frequency for the next generation
        freq = np.random.binomial(pop_size, p_selected) / pop_size
        generations += 1
        
    # return the outcome and number of generations taken to fix/lose
    return freq, generations

def main(args):
    fix_times = []
    loss_times = []
    
    for _ in range(args.replicates):
        outcome, generations = simulate_wright_fisher(args.allele_freq, args.pop_size, args.fitness)
        
        if outcome == 1:  # allele fixed
            fix_times.append(generations)
        else:  # allele lost
            loss_times.append(generations)
    
    # calculate mean and variance for fixation and loss times
    if fix_times:
        mean_fix = np.mean(fix_times)
        var_fix = np.var(fix_times)
        print(f"Allele was fixed in {mean_fix:.2f}. Variance: {var_fix:.2f}")
    
    if loss_times:
        mean_loss = np.mean(loss_times)
        var_loss = np.var(loss_times)
        print(f"Allele was lost in {mean_loss:.2f}. Variance: {var_loss:.2f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate Wright-Fisher allele fixation dynamics with selection.")
    parser.add_argument("--allele_freq", type=float, required=True, help="Initial frequency of the allele (between 0 and 1).")
    parser.add_argument("--pop_size", type=int, required=True, help="Population size (integer).")
    parser.add_argument("--fitness", type=float, required=True, help="Relative fitness of the allele.")
    parser.add_argument("--replicates", type=int, required=True, help="Number of simulation replicates.")
    args = parser.parse_args()
    
    main(args)
