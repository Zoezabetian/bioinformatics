import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

# load the mnist dataset subset
# x contains the pixel data of 6000 images (each 28x28 pixels, flattened into a 784-element array)
# y contains the true digit labels (0-9) for the images
X = np.load('MNIST_X_subset.npy')
y = np.load('MNIST_y_subset.npy')

# visualize the first example
# reshape the first example (1D array) back into a 28x28 image for visualization
first_example = X[0].reshape(28, 28)
plt.imshow(first_example, cmap='gray')  # display as grayscale image
plt.title('first example from MNIST subset')
plt.axis('off')  
# plt.show()  # display the image (hashing it out so that the script can run without interruption)

# function to perform k-means clustering and return centroids and error
def kmeans_clustering(X, y, k):
    kmeans = KMeans(n_clusters=k, random_state=50)  # set arbitrary random state for reproducibility
    kmeans.fit(X)  # fit the k-means model to the dataset
    
    # get centroids
    centroids = kmeans.cluster_centers_  # centroids are the mean positions of each cluster
    
    # predict clusters
    labels = kmeans.labels_  # labels are the predicted cluster assignments for each data point
    
    # compute clustering error
    total_error = 0
    for cluster in range(k):
        # select the true labels of all points assigned to the current cluster
        cluster_labels = y[labels == cluster] # labels == cluster returns a boolean array, y[...] selects the true labels from that array
        if len(cluster_labels) > 0:
            # find the most common true label within the cluster
            most_common_label = np.bincount(cluster_labels).argmax() # bincount counts the occurrences of each digit (label) in the current cluster. this line finds the digit with the highest count, i.e., the most common digit in the cluster. This represents the best guess for the true label of that cluster. 
            
            # count the number of misclassified points (i.e., points that don't match the most common label)
            misclassified_count = np.sum(cluster_labels != most_common_label)
            total_error += misclassified_count
    
    return centroids, total_error

# perform k-means clustering for k=10 and k=11
for k in [10, 11]:
    centroids, error = kmeans_clustering(X, y, k)  # perform k-means and get centroids and error
    
    # reshape and visualize centroids
    centroids_images = centroids.reshape(k, 28, 28)  # reshape each centroid into a 28x28 image
    plt.figure(figsize=(10, 2))  # set the figure size for better display
    for i in range(k):
        plt.subplot(1, k, i + 1)  #  sets up a grid of k subplots in a single row (1 row, k columns)
        plt.imshow(centroids_images[i], cmap='gray')  # display each centroid image in grayscale
        plt.axis('off')  # hide the axis for clarity
    plt.suptitle(f'Centroids for K={k}')  # add a title for the figure
    plt.savefig(f'centroids_k{k}.png')  # save the centroids image as a png file
    plt.close()  # close the figure
    
    # print the clustering error
    print(f'K={k} Error={error}')  # output the clustering error in the specified format
