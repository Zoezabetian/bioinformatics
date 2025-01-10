import sys
import random

# Interval Node
class IntervalNode:
    def __init__(self, interval):
        self.interval = interval  # (chromosome, start, end)
        self.max = interval[2]  # default max is the end of the interval
        self.left = None  # left child
        self.right = None  # right child

# Interval Tree
class IntervalTree:
    def __init__(self):
        self.root = None  # root node

    def build_balanced_tree(self, intervals):
        '''
        Build a balanced interval tree. Called by insert_intervals. Returns the root node of the tree.
        '''
        if not intervals:  # base case for empty list
            return None
        median = len(intervals) // 2  # find median
        root = IntervalNode(intervals[median])  # create root node
        root.left = self.build_balanced_tree(intervals[:median])  # left subtree
        root.right = self.build_balanced_tree(intervals[median + 1:])  # right subtree
        root.max = max(root.interval[2],
                       root.left.max if root.left else float('-inf'),
                       root.right.max if root.right else float('-inf'))
        return root

    def insert_intervals(self, intervals):
        '''
        Insert intervals into the interval tree. Called by __init__. No return value.'''
        merged_intervals = merge_intervals(intervals) # make sure set is merged before building
        self.root = self.build_balanced_tree(merged_intervals)

    def do_overlap(self, interval1, interval2):
        '''
        Check if two intervals overlap. Called by search_overlapping. Returns True if they overlap, False otherwise.'''
        return (interval1[0] == interval2[0] and
        interval1[1] < interval2[2] and
        interval2[1] < interval1[2])

    def search_overlapping(self, root, query_interval, result):
        ''' 
        Search for overlapping intervals in the interval tree recursively. Called by find_overlaps.
        '''
        if root is None: 
            return
        if self.do_overlap(root.interval, query_interval):
            result.append(root.interval)
        # only search left subtree if there's a chance of overlap
        if root.left is not None and root.left.max > query_interval[1]: # 
            self.search_overlapping(root.left, query_interval, result)
        # Search right subtree if the current interval's max is greater than the query start
        # ensures that we only search the right subtree if the starting position of the current interval 
        # is less than the end position of the query interval
        if root.right is not None and root.interval[1] < query_interval[2]: 
            self.search_overlapping(root.right, query_interval, result)

    def find_overlaps(self, query_interval):
        '''
        Find overlapping intervals in the interval tree. Returns a list of overlapping intervals.'''
        result = []
        self.search_overlapping(self.root, query_interval, result)
        return result

# Helper functions
def calculate_overlap_with_tree(set_a, interval_tree):
    '''
    Calculate the overlap between set_a and the intervals in the interval tree. Called by permutation_test. Returns the total overlap.'''
    overlap = 0
    merged_intervals = merge_intervals(set_a)
    
    for query_interval in merged_intervals:
        overlapping_intervals = interval_tree.find_overlaps(query_interval)
        
        for _, overlap_start, overlap_end in overlapping_intervals:
            overlap_start_pos = max(query_interval[1], overlap_start) # find the start position of the overlap
            overlap_end_pos = min(query_interval[2], overlap_end) # find the end position of the overlap
            overlap += overlap_end_pos - overlap_start_pos

    return overlap

def randomize_ranges(set_a, chromosome_length):
    '''
    Randomize the start positions of the intervals in set_a. Called by permutation_test. Returns a list of randomized intervals.'''
    randomized_set_a = []
    for chrom, start, end in set_a:
        interval_length = end - start
        new_start = random.randint(0, chromosome_length - interval_length)
        randomized_set_a.append((chrom, new_start, new_start + interval_length))
    return randomized_set_a

def merge_intervals(intervals):
    '''
    Merge overlapping intervals. Called by permutation_test and insert_intervals. Returns a list of merged intervals.'''
    if not intervals:
        return []
    intervals.sort(key=lambda x: (x[0], x[1]))
    merged = []  
    current = intervals[0]
    for next_interval in intervals[1:]:
        if current[0] == next_interval[0] and current[2] >= next_interval[1]:
            current = (current[0], current[1], max(current[2], next_interval[2]))
        else:
            merged.append(current)
            current = next_interval
    merged.append(current)
    return merged  

def permutation_test(set_a, set_b, chrom_lengths, num_permutations=10000):
    '''
    Main function to perform permutation test. Returns observed overlap and p-value.
    '''
    interval_tree = IntervalTree()
    interval_tree.insert_intervals(set_b)
    observed_overlap = calculate_overlap_with_tree(set_a, interval_tree)
    permuted_overlaps = [0] * num_permutations  # using a list for permuted overlaps
    for i in range(num_permutations):
        shuffled_a = randomize_ranges(set_a, chrom_lengths)
        permuted_overlaps[i] = calculate_overlap_with_tree(shuffled_a, interval_tree)
    p_value = (sum(x >= observed_overlap for x in permuted_overlaps)) / (num_permutations)
    return observed_overlap, p_value

def load_ranges(filename):
    '''
    Load intervals from a file. Returns a list of intervals.'''
    ranges = []  
    with open(filename, 'r') as file:
        for line in file:
            if line.strip():
                parts = line.strip().split('\t')
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                ranges.append((chrom, start, end))  
    return ranges

def parse_fai_file(fai_file):
    '''
    Parse the .fai file to get the chromosome length. Returns the length of the chromosome.'''
    with open(fai_file, 'r') as file:
        for line in file:
            _, length, *_ = line.strip().split()
            return int(length)

def main():    
    if len(sys.argv) != 4 and len(sys.argv) != 5:
        print("Usage: python submission.py path/to/SetA.bed path/to/SetB.bed path/to/genome.fa.fai [num_permutations]")  
        sys.exit(1) 
    set_a_path = sys.argv[1] 
    set_b_path = sys.argv[2]
    fai_path = sys.argv[3]
    num_permutations = int(sys.argv[4]) if len(sys.argv) == 5 else 10000 # set default for num_permutations if not provided
    
    chromosome_length = parse_fai_file(fai_path)
    set_a = load_ranges(set_a_path)
    set_b = load_ranges(set_b_path)
    
    observed_overlap, p_value = permutation_test(set_a, set_b, chromosome_length, num_permutations)
    print(f'Number of overlapping bases observed: {observed_overlap}, p value: {p_value:.4f}')

if __name__ == "__main__":  
    main()
