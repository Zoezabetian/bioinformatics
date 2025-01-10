import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.manifold import MDS


'''PART 1: PCA on MNIST dataset'''
mnist_X = np.load("MNIST_X_subset.npy") # each entry is a 784-dimensional vector representing a 28x28 image (flattened)
mnist_y = np.load("MNIST_y_subset.npy") # each entry is the label of the corresponding image

pca = PCA(n_components=2) # Use PCA to reduce the data from 784 to 2 dimensions.
mnist_X_pca = pca.fit_transform(mnist_X) # Fit PCA to the data and transform it to the 2D space

plt.figure(figsize=(10, 8)) # Visualize 2D transformed data with scatter plot, colored by digit label
scatter = plt.scatter(mnist_X_pca[:, 0], mnist_X_pca[:, 1], c=mnist_y, cmap='viridis', alpha=0.7) # Scatter plot colored by digit label [:,0] is the first column of mnist_X_pca, [:,1] is the second column
plt.colorbar(scatter, label='Digit Label')
plt.title("2D PCA of MNIST Subset") 
plt.xlabel("Principal Component 1")
plt.ylabel("Principal Component 2")
plt.savefig("MNIST_PCA_2D.png")

example_2d = mnist_X_pca[0] # Reconstruct the first image from the dataset using only the first 2 principal components
example_reconstructed = pca.inverse_transform(example_2d) # Now convert the example back to 784 dimensions
plt.imsave("MNIST_original.png", mnist_X[0].reshape(28, 28), cmap='gray') # Save original and reconstructed images (28x28 pixels)
plt.imsave("MNIST_reconstructed_2PC.png", example_reconstructed.reshape(28, 28), cmap='gray')

x = -750  # choose coordinates for "1" in the 2D space
y = 200  
chosen_point = [x, y]
chosen_reconstructed = pca.inverse_transform(chosen_point) # Reconstruct the chosen point back to the original 784-dimensional space
plt.imsave("MNIST_reconstructed_1_from_coord.png", chosen_reconstructed.reshape(28, 28), cmap='gray')



'''Part 2: PCA on Dogs SNP Dataset'''
dogs_X = np.load("dogs_X.npy", allow_pickle=True) # represents the feature matrix derived from single nucleotide polymorphism (SNP) data for dogs
# object arrays cannot be loaded without allow_pickle=True
dogs_clades = np.load("dogs_clades.npy", allow_pickle=True) # represents the clade labels for each dog in the dataset

pca = PCA(n_components=2) # Apply PCA to reduce SNP data from 784 to 2 dimensions
dogs_X_pca = pca.fit_transform(dogs_X)

clades_df = pd.DataFrame(dogs_clades, columns=['clade']) # Convert clades to a pandas DataFrame
unique_clades = clades_df['clade'].unique() # Assign a unique color to each clade
colors = plt.cm.viridis(np.linspace(0, 1, len(unique_clades)))  # Generate colors for each unique clade
clade_color_map = dict(zip(unique_clades, colors))  # Map each clade to a color (zip merges the two lists) and create a dictionary of clade-color pairs


plt.figure(figsize=(12, 8)) 
for clade in unique_clades: # Plot each clade separately in the scatter plot
    indices = clades_df['clade'] == clade # Get indices of samples in the current clade if they match
    plt.scatter(dogs_X_pca[indices, 0], dogs_X_pca[indices, 1],  # Scatter plot colored by clade
                color=clade_color_map[clade], label=clade, alpha=0.7)

plt.title("2D PCA of Dogs SNP Dataset")
plt.xlabel("Principal Component 1")
plt.ylabel("Principal Component 2")
plt.legend(title="Clade Labels", loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8, title_fontsize=10) # bbox_to_anchor places the legend outside the plot area
plt.tight_layout()  # Adjust layout to fit everything within the figure area
plt.savefig("Dogs_PCA_2D.png", bbox_inches='tight')



'''Part 3: MDS on Molecule Distance Matrix'''
distance_data = pd.read_csv("molecule_distances.tsv", sep='\t')
distances = distance_data.iloc[:, 2:].values  # Extract distance matrix (from NMR) by excluding the first two columns
atom_types = distance_data['Element'].values  # Extract atom types

mds = MDS(n_components=3, dissimilarity="precomputed", random_state=40) # Apply MDS to reduce the data from 3D to 2D
molecule_coordinates = mds.fit_transform(distances) # Fit MDS to the precomputed distance matrix

coordinates_df = pd.DataFrame(molecule_coordinates, columns=['x', 'y', 'z']) # Convert coordinates to a pandas DataFrame
coordinates_df['element'] = atom_types
coordinates_df.to_csv("molecule_coordinates.csv", index=False) # Save coordinates to CSV