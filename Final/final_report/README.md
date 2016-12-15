# Modularity and Louvain Algorithm
Final project for Reed College Bio331: Computational Systems Biology.

## Project Summary

The Louvain algorithm is a greedy, iterative algorithm that tries to optimize the clustering of a graph by finding the clustering that results in the highest modularity (a scalar measure of strengths of clusters within a clustering). The goal of the project was to explore using the Louvain algorithm on two datasets, a small, toy dataset with an intuitive, optimum clustering (according to modularity) and the other a larger social network with pre-determined clusterings.

Since node ordering had an impact on the results of the algorithm, node ordering was randomized in order to, first, perform simulations to find the most frequently outputted clustering for a dataset. These most frequent clusterings were not the optimum clustering and so simulations were then repeated to find the most frequent node pairings. These pairings better reflected the optimum and pre-determined clusterings of the toy graph and social network, respectively. Results deteriorated, however, in the social network (due to its size) when simulations exceeded 100.

##Overview
- `toy_weighted.txt` contains a matrix based on six nodes (A, B, C, D, E, and F) and assigned weights between the edge of each node pair.
- `BadgerMatrix.txt` contains a matrix based on badgers and the contact time in seconds between each pairs of badgers. Each badger was a node in the social network graph and the contact time between each badger pair was the weight of each edge connecting the corresponding node pair.
- `final_project_code.py` contains the Python code for running the Louvain algorithm and simulations to find most frequent node pairings.
- `file_utils.py` is a module for commonly used file functions

## Instructions

1. Open final_project_code.py in a directory that includes file_utils.py and the two textfiles. 
2. Type in the name of the desired textfile in line 13.
3. Uncomment line 23 in order to run the Louvain algorithm on the selected textfile and print the results. Note that these results will appear as a dictionary where each cluster from the final clustering is a key. The entry for each key will be the key's neighbors.
4. Uncomment line 26 in order to print the most frequent node pairings for the textfile.

Note that the number of simulations can be set using variable SIMS in line 15.