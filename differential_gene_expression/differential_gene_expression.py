import os
import sys
import statistics
import math

def read_and_normalize_files(directory):
    """
    reads and normalizes gene expression data from all csv files in a given directory and returns as a dictionary.
    """
    
    gene_data_dict = {}  # dict of lists: key = each gene, value = list of all that gene's normalized expression values from all files
    
    for filename in os.listdir(directory):
        if filename.endswith('.csv'):
            filepath = os.path.join(directory, filename)  # get full path of file
            total_expression = 0  # for expression sum
            expressions_list = []  # for (gene, expression) tuples

            with open(filepath, 'r') as file:
                next(file)  # skip header
                for line in file:
                    gene, expression = line.strip().split(',')
                    expression = float(expression)
                    expressions_list.append((gene, expression))  # store gene & expression as tuple
                    total_expression += expression  # add to total expression
            
            # normalize
            for gene, expression in expressions_list:  # iterate over each tuple in expressions list for the given file
                normalized_expression = expression / total_expression  # formula for normalization
                if gene not in gene_data_dict: 
                    gene_data_dict[gene] = []  # initialize empty list to populate with normalized expression values for that gene
                gene_data_dict[gene].append(normalized_expression)  # add normalized expression value from current file to end of list

    return gene_data_dict  # return dict of lists containing normalized expressions for each gene ({gene: [normalized_expressions])

def get_mean_median_tuple(gene_data):
    """
    computes the mean and median normalized expressions for each gene and returns as a tuple of means and medians dictionaries.
    """
    
    means_dict = {}  # key = gene, value = mean of normalized expression values
    medians_dict = {}  # key = gene, value = median of normalized expression values

    for gene, expressions in gene_data.items():
        means_dict[gene] = statistics.mean(expressions)  # calculate mean of expressions for each gene and store in means dict
        medians_dict[gene] = statistics.median(expressions)  # calculate median of expressions for each gene and store in medians dict

    return means_dict, medians_dict  # return tuple: ({gene: mean}, {gene: median})

def compute_log2_fold_change(control_means, treatment_means):
    """
    computes the log2 fold change between control and treatment means for each gene.
    """
    
    log2_fold_changes_dict = {}  # empty dict to store log2 fold changes

    for gene in control_means:
        if gene in treatment_means:  # check if gene is present in both control and treatment means
            control_mean = control_means[gene]  # get mean for control group
            treatment_mean = treatment_means[gene]  # get mean for treatment group
            
            # calculate log2 fold change
            log2_fold_change = math.log2(treatment_mean / control_mean) if control_mean > 0 else float('inf')  # formula: take log2() of the ratio of treatment to control means
            log2_fold_changes_dict[gene] = log2_fold_change  # store in dict for each gene

    return log2_fold_changes_dict  # return dict: ({gene: log2_fold_change})

def write_results_to_file(log2_fold_changes_dict, control_means_dict, control_medians_dict, treatment_means_dict, treatment_medians_dict):
    """
    writes the results to a file. each line contains one gene's data.
    """
    
    output_strings_list = [
        "gene\tmean_normalized_control_expression\tmedian_normalized_control_expression\tmean_normalized_treatment_expression\tmedian_normalized_treatment_expression\tlogFoldChange"
    ]
    
    # make a list of tuples from the log2_fold_changes_dict
    sorted_genes = sorted(log2_fold_changes_dict.items(), key=lambda x: x[1])  # sort by the log fold change values
    
    # loop over the sorted genes and append each gene's data to output_strings_list
    for gene, log_fold_change in sorted_genes:
        output_strings_list.append(  
            f"{gene}\t{control_means_dict[gene]}\t{control_medians_dict[gene]}\t{treatment_means_dict[gene]}\t{treatment_medians_dict[gene]}\t{log_fold_change}"
        )

    # write output to file
    with open("output1.txt", "w") as outfile:
        outfile.write("\n".join(output_strings_list))

def main(control_dir, treatment_dir):
    """
    main function to read gene expression data, compute statistics, and write results to a file.
    """
    
    # read and normalize gene expression data
    control_data_dict = read_and_normalize_files(control_dir)  # get control gene data dict
    treatment_data_dict = read_and_normalize_files(treatment_dir)  # get treatment gene data dict
    
    # compute statistics for output
    control_means_dict, control_medians_dict = get_mean_median_tuple(control_data_dict)  # get (means dict, medians dict) for control
    treatment_means_dict, treatment_medians_dict = get_mean_median_tuple(treatment_data_dict)  # get (means dict, medians dict) for treatment
    log2_fold_changes_dict = compute_log2_fold_change(control_means_dict, treatment_means_dict)  # get log2 fold change dict

    # write output to file
    write_results_to_file(log2_fold_changes_dict, control_means_dict, control_medians_dict, treatment_means_dict, treatment_medians_dict)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python your_assignment.py path/to/controls path/to/treatments")  
        sys.exit(1) 
        
    control_directory = sys.argv[1] 
    treatment_directory = sys.argv[2]  

    main(control_directory, treatment_directory)
