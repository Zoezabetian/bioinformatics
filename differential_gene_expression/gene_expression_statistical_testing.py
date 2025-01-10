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
            
            for gene, expression in expressions_list:  # iterate over each tuple in expressions list for the given file
                normalized_expression = expression / total_expression  # formula for normalization
                if gene not in gene_data_dict:
                    gene_data_dict[gene] = []  # initialize empty list to populate with normalized expression values for that gene
                gene_data_dict[gene].append(normalized_expression)  # add normalized expression value from current file to end of list

    return gene_data_dict  # dict of lists ({gene: [normalized_expressions]})

def get_mean_median_tuple(gene_data):
    """
    computes the mean and median normalized expressions for each gene and returns as a tuple of means and medians dictionaries.
    """
    
    means_dict = {}  # key = gene, value = mean of normalized expression values
    medians_dict = {}  # key = gene, value = median of normalized expression values

    for gene, expressions in gene_data.items():
        means_dict[gene] = statistics.mean(expressions)  # calculate mean of expressions for each gene and store in means dict
        medians_dict[gene] = statistics.median(expressions)  # calculate median of expressions for each gene and store in medians dict

    return means_dict, medians_dict  # return tuple of 2 dicts: ({gene: mean}, {gene: median})

def compute_log2_fold_change(control_means, treatment_means):
    """
    computes the log2 fold change between control and treatment means for each gene.
    """
    
    log2_fold_changes_dict = {}  # empty dict to store log2 fold changes

    for gene in control_means:
        if gene in treatment_means:
            control_mean = control_means[gene]  # get means for each group
            treatment_mean = treatment_means[gene]  
            
            log2_fold_change = math.log2(treatment_mean / control_mean) if control_mean > 0 else float('inf')
            log2_fold_changes_dict[gene] = log2_fold_change  # store in dict for each gene

    return log2_fold_changes_dict  # dict of log2 fold changes ({gene: log2_fold_change}).

def mann_whitney_u_test(control_expr_list, treatment_expr_list):
    """
    perform Mann-Whitney U test and return p-value for a gene's expression values in control and treatment groups.
    """
    
    num_control_obs = len(control_expr_list)
    num_treatment_obs = len(treatment_expr_list)
 
    combined_expr_list = control_expr_list + treatment_expr_list  # combine lists
    ranks = {expr: rank for rank, expr in enumerate(sorted(combined_expr_list), start=1)}  # assign ranks by sorting and enumerating
    
    rank_sum = sum(ranks[expr] for expr in control_expr_list)  # calculate the rank sum for the control group
    control_U = rank_sum - (num_control_obs * (num_control_obs + 1)) / 2  # U statistic for control
    treatment_U = num_treatment_obs * num_treatment_obs - control_U  # U statistic for treatment
    U = min(control_U, treatment_U)  # Mann-Whitney U statistic

    # calculate p-value using approximation
    expected_U = (num_control_obs * num_treatment_obs) / 2  # expected U under the null hypothesis
    std_U = math.sqrt((num_control_obs * num_treatment_obs * (num_control_obs + num_treatment_obs + 1)) / 12)  # standard deviation of U
    z_score = (U - expected_U) / std_U if std_U > 0 else 0
    p_value = 2 * (1 - statistics.NormalDist().cdf(abs(z_score)))  # approximate using standard normal distribution

    return p_value

def main(control_dir, treatment_dir):
    """
    main function to read gene expression data, compute statistics, and write results to a file.
    """
    
    control_data_dict = read_and_normalize_files(control_dir)
    treatment_data_dict = read_and_normalize_files(treatment_dir)
    
    control_means_dict, control_medians_dict = get_mean_median_tuple(control_data_dict)
    treatment_means_dict, treatment_medians_dict = get_mean_median_tuple(treatment_data_dict)
    log2_fold_changes_dict = compute_log2_fold_change(control_means_dict, treatment_means_dict)

    results = []

    for gene in log2_fold_changes_dict:
        mean_control = control_means_dict[gene]
        median_control = control_medians_dict[gene]
        mean_treatment = treatment_means_dict[gene]
        median_treatment = treatment_medians_dict[gene]
        log2_fold_change = log2_fold_changes_dict[gene]
        
        control_expressions_list = control_data_dict[gene] if gene in control_data_dict else []
        treatment_expressions_list = treatment_data_dict[gene] if gene in treatment_data_dict else []
        p_value = mann_whitney_u_test(control_expressions_list, treatment_expressions_list)

        results.append((gene, mean_control, median_control, mean_treatment, median_treatment, log2_fold_change, p_value))
    
    results.sort(key=lambda x: x[-1])  # sort results by p-value

    output_lines = ["gene\tmean_normalized_control_expression\tmedian_normalized_control_expression\tmean_normalized_treatment_expression\tmedian_normalized_treatment_expression\tlogFoldChange\tp_value"]

    for result in results:
        output_lines.append("\t".join(map(str, result))) 
    
    with open("output2.txt", "w") as outfile:
        outfile.write("\n".join(output_lines))  

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python your_assignment.py path/to/controls path/to/treatments") 
        sys.exit(1) 

    control_directory = sys.argv[1] 
    treatment_directory = sys.argv[2]  

    main(control_directory, treatment_directory)
