# -*- coding: utf-8 -*-


# import libraries
import matplotlib.pyplot as plt
import numpy as np
import random

# define variables for grid dimensions and site occupation probability
x = 50
y = 50
p = 0.59


def build_grid(y_dimension, x_dimension, probability):
    '''
    Generates grid and fills with 1 or 0 dependent on probability.
    
    Additional bounds of zeroes are added to deal with edge cases.
    '''
    
    # define the grid
    grid = np.zeros((y_dimension + 2, x_dimension + 2))
    
    # fill sites with either 1 or 0
    for site in np.nditer(grid[1:-1, 1:-1], op_flags=['readwrite']):
        if random.random() < probability:
            site[...] = 1
        else:
            pass
        
    return grid


def identify_clusters(grid):
    '''
    Find clusters by identifying which sites belong to the same cluster.
    
    Sites neighbouring an occupied site will be relabelled with the label of that site.
    '''
    
    # find all sites which are occupied 
    occupied_sites = np.count_nonzero(grid)
    
    # label all occupied sites with number and finds coords of each
    site_labels = np.arange(occupied_sites)
    site_coords = [list(x) for x in np.argwhere(grid == 1)]
    
    while True:
        
        # array to store site labels which have not yet been changed
        uncorrected_labels = []
         
        for i in np.arange(site_labels.size):
             
             # extract coordinates of ith site
             y, x = site_coords[i]
             
             # relabel site with label of occupied neighbour
             if grid[y-1][x] == 1 and grid[y][x-1] == 0:
                 site_labels[i] = site_labels[site_coords.index([y-1, x])]
             elif grid[y][x-1] == 1 and grid[y-1][x] == 0:
                 site_labels[i] = site_labels[site_coords.index([y, x-1])]
              
             # relabel site with label of either occupied neighbours
             elif grid[y-1][x] == 1 and grid[y][x-1] == 1:
                 neighbour1_label = site_labels[site_coords.index([y-1, x])]
                 neighbour2_label = site_labels[site_coords.index([y, x-1])]
                 site_labels[i] = np.min([neighbour1_label, neighbour2_label])
                 
                 # store labels if neighbour labels are different 
                 if neighbour1_label != neighbour2_label:
                    uncorrected_labels.append([neighbour1_label, neighbour2_label])
         
         # quit loop when all labels are corrected
        if uncorrected_labels == []:
            break
        
         # otherwise correct labels that were stored earlier 
        else:
            for label1, label2 in uncorrected_labels:
                wrong_id = np.max([label1, label2])
                correct_id = np.min([label1, label2])
                site_labels[site_labels == wrong_id] = correct_id
                
    return site_labels, site_coords


def test_percolation(site_labels, site_coords, y_dimension, x_dimension):
    '''
    Finds out if percolation has occured.
    '''
    
    cluster_coords = []
    
    # find coordinates of each site in a cluster 
    for cluster_label in np.unique(site_labels):
        cluster_coords.append([site_coords[k] for k in range(len(site_labels)) if site_labels[k] == cluster_label])
    
    # search for percolated cluster(s)
    percolation = []
    for cluster in cluster_coords:
        cluster = np.array(cluster).T
        if (1 in cluster[0]) and (y_dimension in cluster[0]):
            percolation.append(cluster)
        if (1 in cluster[1]) and (x_dimension in cluster[1]):
            percolation.append(cluster)

    if len(percolation) > 0:
        return "Percolation", percolation
    else:
        return "No Percolation", None


# define populated grid 
grid = build_grid(y, x, p)

# identify clusters and their positions
clusters = identify_clusters(grid)

# run test to see if percolation has occured
status, perc_clusters = test_percolation(clusters[0], clusters[1], y, x)
print(status)

# plot grid
plt.figure(figsize=(6, 6))
plt.imshow(grid[1:-1, 1:-1], cmap='gray_r', origin='upper')

# indentify largest percolating cluster
if perc_clusters:
    largest_cluster = max(perc_clusters, key=lambda c: c.shape[1])
    for (y_coord, x_coord) in largest_cluster.T:
        plt.gca().add_patch(plt.Rectangle((x_coord - 1.5, y_coord - 1.5), 1, 1, color='red', alpha=0.6))

# configure grid 
plt.xticks([])
plt.yticks([])
plt.show()