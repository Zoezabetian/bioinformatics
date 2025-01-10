import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

molecule_coordinates = pd.read_csv("molecule_coordinates.csv") # Load MDS coordinates and atom types from coordinates file
distance_matrix = pd.read_csv("molecule_distances.tsv", sep='\t').iloc[:, 2:].values # Load distance matrix from NMR data; '.iloc[:, 2:]' excludes the first two columns

cpk_colors = {
    1: ('white', 'Hydrogen'),     # Hydrogen
    6: ('black', 'Carbon'),       # Carbon
    8: ('red', 'Oxygen'),         # Oxygen
}

atomic_sizes = {    # Define atomic sizes for marker scaling (atomic number as key, size as value)
    1: 50,    # Hydrogen
    6: 100,   # Carbon
    8: 140,   # Oxygen
}

fig = plt.figure(figsize=(10, 8)) # 3D plot setup
ax = fig.add_subplot(111, projection='3d') # Add a 3D subplot to the figure

for i, row in molecule_coordinates.iterrows(): # Plot atoms in 3D with CPK color and adjusted size
    element = row['element']
    x, y, z = row['x'], row['y'], row['z']
    color, name = cpk_colors.get(element, ('gray', 'Unknown'))  # Default to gray and 'Unknown' if element not in dictionary. .get method returns the value for the given key, or a default value if the key is not found.
    size = atomic_sizes.get(element, 80)     # Default size if element not in dictionary
    ax.scatter(x, y, z, color=color, s=size, edgecolor='k', alpha=0.8) # plot the atom with the specified color, size, and edge color

for i in range(len(molecule_coordinates)):
    for j in range(i + 1, len(molecule_coordinates)): # for each pair of atoms (i+1 to avoid self-comparison)
        if distance_matrix[i, j] < 1.6: # Draw bonds between atoms if distance < 1.6
            x_vals = [molecule_coordinates['x'][i], molecule_coordinates['x'][j]] # Get x values for the bond 
            y_vals = [molecule_coordinates['y'][i], molecule_coordinates['y'][j]] # Get y values for the bond
            z_vals = [molecule_coordinates['z'][i], molecule_coordinates['z'][j]] # Get z values for the bond
            ax.plot(x_vals, y_vals, z_vals, color='gray', linewidth=0.8, alpha=0.6) # draw a line between the two atoms

for element, (color, name) in cpk_colors.items(): # Add a legend for atom types with correct names
    ax.scatter([], [], [], color=color, s=100, edgecolor='k', alpha=0.8, label=name)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("3D Visualization of Molecular Structure")
ax.legend(title="Atom Legend")
plt.savefig("molecule_3D_plot.png", bbox_inches='tight')