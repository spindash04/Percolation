# -*- coding: utf-8 -*-


# import libraries
import matplotlib.pyplot as plt
import numpy as np

# define variables for grid dimension and bonding probability
N = 10
p = 0.5

# randomly generate vertical and horizontal bonds
right = np.random.rand(N, N-1) < p  
down  = np.random.rand(N-1, N) < p  

# create parent array
parent = np.arange(N * N) 


# disjoint union algorithm for efficient relabelling
def find(x):
    '''
    Runs up the line of neighbouring bonds to find the site that is its own parent.

    Once 'root' site is found all of the neighbouring bonds are changed to that label.

    '''
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(x, y):
    '''
    Resolves two branches.

    '''
    rx, ry = find(x), find(y)
    if rx != ry:
        parent[ry] = rx  


# loop through each site in the lattice and unify connected sites
for i in range(N):
    for j in range(N):
        idx = i * N + j
        
        # union site with its rightward bond
        if j < N - 1 and right[i, j]:
            union(idx, i * N + (j + 1)) 
        
        # union site with its downward bond 
        if i < N - 1 and down[i, j]:
            union(idx, (i + 1) * N + j)

# create array of cluster labels for each site
clusters = np.zeros((N, N), dtype=int)

# map unique roots to a new cluster label
root_to_cluster = {}
next_cluster = 1

for i in range(N):
    for j in range(N):
        # find root of parent site
        root = find(i * N + j)
        
        # if root not yet assigned to cluster label then establish one
        if root not in root_to_cluster:
            root_to_cluster[root] = next_cluster
            next_cluster += 1
        
        # assign new cluster label to site
        clusters[i, j] = root_to_cluster[root]


# plot grid
fig, ax = plt.subplots(figsize=(6, 6))

for i in range(N):
    for j in range(N):
        
        # draw rightward bonds
        if j < N - 1 and right[i, j]:
            ax.plot([j, j + 1], [i, i], color='black', lw=2)
        
        # draw downward bonds
        if i < N - 1 and down[i, j]:
            ax.plot([j, j], [i, i + 1], color='black', lw=2)

# plot sites as circles
for i in range(N):
    for j in range(N):
        ax.plot(j, i, 'o', color='white', markersize=6, markeredgecolor='black')

# configure grid
ax.set_aspect('equal')
ax.invert_yaxis()
ax.set_xticks([])
ax.set_yticks([])
plt.tight_layout()
plt.show()