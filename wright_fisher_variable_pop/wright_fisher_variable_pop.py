import argparse
import numpy as np

def load_population_schedule(tsv_path):
    """Load the population size schedule from a TSV file manually without pandas."""
    schedule = {}
    with open(tsv_path, 'r') as file:
        next(file) 
        for line in file:
            generation, popsize = map(int, line.strip().split('\t'))
            schedule[generation] = popsize
    return schedule

def simulate_fixation_loss_variable_pop(allele_freq, fitness, schedule):
    generations = 0
    fixed = False
    lost = False
    pop_size = schedule.get(0, 100)  # initial population size from generation 0, default 100 if not in schedule

    while not fixed and not lost:
        # update population size if the current generation matches a change point in the schedule
        if generations in schedule:
            pop_size = schedule[generations]

        # simulate allele inheritance with selection
        allele_count = np.random.binomial(pop_size, allele_freq * fitness / (allele_freq * fitness + (1 - allele_freq)))
        allele_freq = allele_count / pop_size

        generations += 1
        if allele_freq == 1:
            fixed = True
        elif allele_freq == 0:
            lost = True

    return generations, fixed

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--allele_freq", type=float, required=True)
    parser.add_argument("--fitness", type=float, required=True)
    parser.add_argument("--replicates", type=int, required=True)
    parser.add_argument("--pop_size_file", type=str, required=True)
    args = parser.parse_args()

    # load population schedule
    schedule = load_population_schedule(args.pop_size_file)  # updated reference

    # run replicates
    fixation_times = []
    loss_times = []
    for _ in range(args.replicates):
        generations, fixed = simulate_fixation_loss_variable_pop(args.allele_freq, args.fitness, schedule)
        if fixed:
            fixation_times.append(generations)
        else:
            loss_times.append(generations)

    if fixation_times:
        mean_fixation = np.mean(fixation_times)
        var_fixation = np.var(fixation_times)
        print(f"Allele was fixed in {mean_fixation:.2f} generations. Variance: {var_fixation:.2f}")
    if loss_times:
        mean_loss = np.mean(loss_times)
        var_loss = np.var(loss_times)
        print(f"Allele was lost in {mean_loss:.2f} generations. Variance: {var_loss:.2f}")

if __name__ == "__main__":
    main()
