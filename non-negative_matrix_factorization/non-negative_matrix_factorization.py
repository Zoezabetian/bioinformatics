import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
import matplotlib.pyplot as plt


def load_data(dog_X_path, dog_clades_path):
    '''Load the SNP data and clades'''
    dog_X = np.load(dog_X_path)  # dog matrix has each row as a dog and each column as a SNP feature
    dog_clades = np.load(dog_clades_path, allow_pickle=True)  # object arrays need allow_pickle=True
    return dog_X, dog_clades


def apply_nmf(dog_X, n_components=5, max_iter=500):
    ''' The apply_nmf function breaks down each dog’s SNP data into a smaller set of main patterns, called components, using Non-negative 
        Matrix Factorization (NMF). This process produces two matrices, W and H, where the original data is approximately equal to W times H.
        ** Matrix H: the main patterns found in the data, with each row showing how much a specific SNP feature contributes to a component.
        ** Matrix W: each dog as a unique mix of components, with values indicating how much each component contributes to that dog’s profile.
        In short, NMF reveals each dog’s genetic data as a combination of shared patterns, helping to uncover hidden relationships.
    ''' 
    
    model = NMF(n_components=n_components, init='random', random_state=42, max_iter=max_iter)
    W = model.fit_transform(dog_X) # W: each row corresponds to a dog, and each column represents a component. 
                                   # for each dog, the sum of proportions across all components gives an overall representation of that dog.
                                   # (proportion is how much that component contributes to the dog's genetic profile)
    H = model.components_ # H: each row represents one component (5 rows), and each column represents a set of SNP features (extracted from the data by NMF)
    return W, H 


def normalize_for_visual(W):
    ''' Normalize the W matrix to represent proportions, (so the components for each individual sum to 1).
        Allows each dog’s representation in W to be viewed as a distribution across components, 
        making it easier to interpret relative importance across components for each clade
    '''
    W_normalized = W / W.sum(axis=1, keepdims=True) # divide each element in a row by the row’s sum ensures that each row’s total is 1, converting component values to proportions
    return W_normalized


def sort_by_dominant_cluster(W_normalized, dog_clades):
    ''' Determine the dominant (largest value) cluster for each dog.
        Sort the dogs by dominant cluster and proportion (highest to lowest).
    '''
    data_df = pd.DataFrame(W_normalized, columns=[f'Component {i+1}' for i in range(W_normalized.shape[1])]) # convert to DataFrame
    data_df['Clade'] = dog_clades # add clade information
    data_df['Dominant_Component'] = data_df.iloc[:, :-1].idxmax(axis=1)  # finds the column (component) with the highest value for each row (dog) (axis=1 to search across columns)
    data_df['Dominant_Proportion'] = data_df.iloc[:, :-2].max(axis=1)    #  maximum value in each row, which corresponds to the proportion of the dominant component

    sorted_df = data_df.sort_values(by=['Dominant_Component', 'Dominant_Proportion'], ascending=[True, False])   # sorts the dogs first by the dominant component and then by the proportion of that dominant component (in descending order)
    sorted_df = sorted_df.drop(columns=['Dominant_Component', 'Dominant_Proportion'])   # removes the temporary columns which were used only for sorting
    return sorted_df


def plot_sorted_dogs(W_normalized_sorted, output_file='NMF_Dogs.png'):
    ''' Create a stacked plot where each bar represents a dog, and y-axis represents the proportion to each cluster
        Visualize the bar plot
    '''
    components = [W_normalized_sorted.iloc[:, i].values for i in range(W_normalized_sorted.shape[1] - 1)]  # :, i selects all rows and the i-th column for i in range of all columns except the last two
    x = np.arange(W_normalized_sorted.shape[0])  # x-axis represents individual dogs
    fig, ax = plt.subplots(figsize=(15, 8))
    ax.stackplot(x, *components, labels=[f'Component {i+1}' for i in range(len(components))], alpha=0.7) # *components unpacks the list of arrays into individual arguments

    ax.set_xlabel('Dogs (Sorted by Dominant Cluster and Proportion)', fontsize=12)
    ax.set_ylabel('Proportion of Components (Normalized)', fontsize=12)
    ax.set_title('Proportion of NMF Components for Each Dog (Sorted)', fontsize=14)
    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


def identify_dominant_component(W_normalized, dog_clades):
    '''Identify the dominant component in Basenji and Wolf samples
    '''
    data_df = pd.DataFrame(W_normalized, columns=[f'Component {i+1}' for i in range(W_normalized.shape[1])]) # convert to DataFrame
    data_df['Clade'] = dog_clades # add clade information
    basenji_samples = data_df[data_df['Clade'] == '**Basenji'] # get Basenji samples
    wolf_samples = data_df[data_df['Clade'] == 'Wolf'] # get Wolf samples
    
    basenji_dominant_components = basenji_samples.drop(columns=['Clade']).idxmax(axis=1) # find dominant component for each Basenji sample
    wolf_dominant_components = wolf_samples.drop(columns=['Clade']).idxmax(axis=1) # find dominant component for Wolf samples
    
    basenji_dominant = basenji_dominant_components.mode()[0]  # get the most common dominant component from all Basenji samples
    wolf_dominant = wolf_dominant_components.mode()[0] # get the most common dominant component from all Wolf samples
    
    print(f"The dominant component in Basenji samples is: {basenji_dominant}")
    print(f"The dominant component in Wolf samples is: {wolf_dominant}")
    if basenji_dominant == wolf_dominant:
        print(f"Both Basenji and Wolf samples have the same dominant component: {basenji_dominant}")
    else:
        print("The dominant components differ between Basenji and Wolf samples.")


if __name__ == "__main__":
    dog_X_path = sys.argv[1]
    dog_clades_path = sys.argv[2]
    
    # Load matrix where each row is a dog and each column an SNP feature, and load clade labels
    dog_X, dog_clades = load_data(dog_X_path, dog_clades_path)
    
    # Apply NMF with 5 components to decompose the dataset
    # W: each row = dog, each column = component, sum of proportions across all components gives an overall representation of that dog
    # H: each row = component, each column = set of SNP features, indicates how each SNP feature contributes to each component
    # Matrix H could help if you wanted to interpret the nature of each component (e.g., understanding which features define each cluster). 
    # However, this is not required for simply clustering and visualizing
    W, H = apply_nmf(dog_X, n_components=5)
    
    # Normalize the W matrix to represent proportions, (so the components for each individual sum to 1)
    W_normalized = normalize_for_visual(W)
    
    # Determine the dominant (largest value) cluster for each dog
    # Sort the dogs by dominant cluster and proportion (highest to lowest)
    W_normalized_sorted = sort_by_dominant_cluster(W_normalized, dog_clades)
    
    # Create a stacked plot where each bar represents a dog, and y-axis represents the proportion to each cluster and visualize
    plot_sorted_dogs(W_normalized_sorted, output_file='NMF_Dogs.png')

    # Uncomment the following line to identify the dominant component in Basenji and Wolf samples
    identify_dominant_component(W_normalized, dog_clades)

