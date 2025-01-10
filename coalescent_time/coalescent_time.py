import argparse
import numpy as np

def simulate_coalescent_time(pop_size, sample_size, target_event):
    """
    Simulate the time to reach a specific coalescent event in a population.
    """
    num_lineages = sample_size
    total_time = 0
    event_count = 0
    
    while num_lineages > 1 and event_count < target_event:
        # calculate the rate for the next coalescent event
        rate = 2 * pop_size / (num_lineages * (num_lineages - 1))
        
        # draw time until next coalescent event from exponential distribution
        time_to_next_event = np.random.exponential(rate)
        total_time += time_to_next_event
        event_count += 1
        
        # reduce the number of lineages by 1 due to coalescence
        num_lineages -= 1
    
    return total_time

def main(args):
    coalescent_times = []
    
    for _ in range(args.replicates):
        time_to_eighth_event = simulate_coalescent_time(args.pop_size, args.sample_size, target_event=8)
        coalescent_times.append(time_to_eighth_event)
    
    # calculate mean and variance for the eighth coalescent event
    mean_time = np.mean(coalescent_times)
    variance_time = np.var(coalescent_times)
    
    # output the results
    print(f"Time to eighth coalescent event: {mean_time:.2f}. Variance: {variance_time:.2f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate coalescent times for a Wright-Fisher population.")
    parser.add_argument("--pop_size", type=int, required=True, help="Effective population size.")
    parser.add_argument("--sample_size", type=int, required=True, help="Sample size (number of lineages).")
    parser.add_argument("--replicates", type=int, required=True, help="Number of simulation replicates.")
    args = parser.parse_args()
    
    main(args)
