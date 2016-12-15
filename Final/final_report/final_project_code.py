##Mina Marden
##BIO 331 Final Project


from file_utils import *
import copy
import numpy as np

##note: BadgerMatrix.txt was edited to work with code ('Badger' was removed from top left)

def main():
	##TEXTFILE is the file that contains the data on the network
	TEXTFILE = 'TEXTFILE'
	##SIMS is the number of desires simulations
	SIMS = 100
	
	##creates variables based on the information within TEXTFILE
	##True is inputted if TEXTFILE is a matrix and False, otherwise
	nodes, edges, matrix, neighbors = convert_file(TEXTFILE, True)
	edge_weights = edge_weight_finder(nodes, matrix)

	##uncomment below to run louvain on TEXTFIILE
	#print(louvain(nodes, matrix, neighbors, edge_weights))

	##uncomment below to find the most frequently paired nodes by the louvain function after SIMS simulations
	#print(d_freq_finder(SIMS, nodes, nodes, matrix, neighbors, edge_weights))
	return

##runs the louvain algorithm for a pre-specified number of simulations and then generates
##a clustering based on the frequency of node pairs appearing in the same cluster
def d_freq_finder(sims, orig_nodes, nodes, matrix, neighbors, edge_weights):
	##creates a dictionary called pairs where each key is a node in the network and
	##the entry for each key is a sub-dictionary where each sub-key is a node in the
	##network that is not the same as the origianl key in the main dictionary
	pairs = {}
	ordered_nodes = orig_nodes[::]
	ordered_nodes.sort()
	tup_in_nodes = False
	for node1 in ordered_nodes:
		pairs[node1] = {}
		if type(node1)==tuple:
			tup_in_nodes = True
		for node2 in ordered_nodes:
			if node1==node2:
				pass
			else:
				pairs[node1][node2] = 0
	keys = []
	checked = []
	for i in range(sims):
		print('i:', i)
		##each simulation computes the louvain algorithm
		mod, d = louvain(nodes,matrix,neighbors, edge_weights)
		##then the clusterings are made into a list of clusters where
		##each element in the cluster is a node in the network
		##example: [(('A','C'),'B'),('D','F','E')] becomes [('A','B','C'),('D','E','F')]
		clusters = []
		for key in d:
			entry = non_tup_finder(key, key, tup_in_nodes)
			entry.sort()
			if entry not in clusters:
				if len(entry)==1 and type(entry[0])==tuple:
					clusters.extend(entry)
				else:
					clusters.append(tuple(entry))
		clusters.sort()
		
		##computes the frequency of each pairing of nodes being in the same cluster
		##based on whether the nodes in the network are tuples or strings
		for c in clusters:
			if tup_in_nodes==False:
				for node1 in c:
					for node2 in c:
						if node1!=node2:
							pairs[node1][node2]+=0.5
							pairs[node2][node1]+=0.5
			else:
				if type(c[0])==tuple:
					for node1 in c:
						for node2 in c:
							if node1!=node2:
								pairs[node1][node2]+=0.5
								pairs[node2][node1]+=0.5
	##creates a dictionary called max_keys where the keys are the nodes in the network
	##and the entry for each key is the list of the node that the key was most frequently
	##paired with during the clustering in the simulations
	max_keys = {}
	for node in pairs:
		max_keys[node] = [node]
	for node in pairs:
		d = pairs[node]
		max_val = max(d.values())
		for key in d:
			if d[key]==max_val:
				max_keys[node].append(key)
				max_keys[key].append(node)
	
	##creates a list called new_clusters based on the entries of max_keys
	new_clusters = []
	for node in max_keys:
		##orders and removes repeats from the lists in the entry of each key
		poss_clus = tuple(sorted(tuple(set(max_keys[node]))))
		if poss_clus not in new_clusters:
			new_clusters.append(poss_clus)

	##creates a list called final_clusters which is identical to new_clusterings, except that
	##clusters that are subsets of other clusters in new_clusters are not present
	##example: if new_clusters = [('A','B','C'),('A','B')], then
	##final_clusters = [('A','B','C')]
	removed = []
	final_clusters = []
	for clust1 in new_clusters:
		remove_truth = False
		for clust2 in new_clusters:
			if clust1!=clust2 and remove_truth==False:
				set1 = set(clust1)
				set2 = set(clust2)
				if set1.issubset(set2) and clust1 not in removed:
					final_clusters.append(clust1)
					removed.append(clust1)
					remove_truth = True

	return final_clusters

##recursive function that takes a cluster c (a tuple of tuples) and creates a list el_list that
##contains all nodes within the tuples of c, based on whether the nodes are tuples or strings
def non_tup_finder(c, prev_c, tup_in_nodes):
	el_list = []
	if type(c)!=tuple:
		return el_list
	elif tup_in_nodes==False:
		for el in c:
			if type(el)==tuple:
				el_list.extend(non_tup_finder(el, c, tup_in_nodes))
			else:
				if el not in el_list:
					el_list.extend([el])
		el_list.sort()
		tuple(el_list)
	elif tup_in_nodes:
		for el in c:
			if type(el)==tuple:
				el_list.extend(non_tup_finder(el, c, tup_in_nodes))
			else:
				if c not in el_list:
					el_list.extend([c])
		el_list.sort()
		tuple(el_list)
	return el_list

##creates a dictionary of node degrees
def node_deg_finder(nodes, edges):
	node_degs = {}
	for node in nodes:
		node_degs[node] = 0
	for edge in edges:
		w = edges[edge]/2.0
		node_degs[edge[0]]+=w
		node_degs[edge[1]]+=w
	return node_degs

##orders the elements of a dictionary low-to-high if low_to_high = True and high-to-low, otherwise
def order(d, low_to_high):
	key_list = []
	val_list = []
	checked = []
	for key in d:
		key_list.append(key)
		val_list.append(d[key])
	val_list.sort()
	new_key_list = []
	for val in val_list:
		for key in key_list:
			if d[key]==val and key not in checked:
				new_key_list.append(key)
				checked.append(key)
	if low_to_high==False:
		val_list = val_list[::-1]
		new_key_list = new_key_list[::-1]
	return new_key_list, val_list

##converts a datafile on a network into a set of nodes, an adjacency matrix, and then a dictionary for
##edges as well as a dictionary for neighbors
def convert_file(textfile_name, textfile_type):
	if textfile_type:
		matrix, nodes = readMatrix(textfile_name)
		edges = edges_maker(nodes, matrix)
		neighbors = edges.copy()
	else:
		nodes, edges, neighbors = readData(textfile_name)
		matrix = []
		for node1 in nodes:
			row = []
			for node2 in nodes:
				if node1!=node2:
					for edge in edges[node1]:
						present = False
						if node2 in edge:
							row.append(float(edge[1]))
							present = True
					if present==False:
						row.append(0)
			matrix.extend([row])

	return nodes, edges, matrix, neighbors

##creates a dictionary of edges based on an adjacency matrix
def edges_maker(nodes, matrix):
	edges = {}
	for i in range(len(matrix)):
		node1 = nodes[i]
		neighs = []
		for j in range(len(matrix[0])):
			node2 = nodes[j]
			if matrix[i][j]!=0:
				neighs.append(node2)
		edges[node1] = neighs
	return edges

##modifies an adjacency matrix after the original nodes in the network are merged
def matrix_changer(new_nodes, old_nodes, matrix):
	new_matrix = []
	index_order = []
	for new_node in new_nodes:
		for i in range(len(old_nodes)):
			old_node = old_nodes[i]
			if new_node==old_node:
				new_matrix.append(matrix[i])
				index_order.append(i)
	for j in range(len(new_matrix)):
		row = new_matrix[j]
		new_row = row[::]
		i = 0
		for index in index_order:
			new_row[i] = row[index]
			i+=1
		new_matrix[j] = new_row
	return new_matrix

##Implementation of the louvain algorithm that greedily finds the optimum clustering of a network
##First, it finds the highest increase in modularity after iteratively moving a node into the cluster of
##its neighbors. It implements the highest increase per node and continues onto the next node
##until each node has been checked. It then checks if there was an overall gain in modularity from the 
##previous clustering. If this is true, then the function repeats until there is not an overall increase
##in modularity and returns the initial clustering that resulted in no further increase.
def louvain(nodes, matrix, neighbors, edge_weights):
	##computes the sum of the weights in the network in order to compute modularity
	total_weight = 0
	for row in matrix:
		total_weight+=sum(row)
	total_weight = float(total_weight)/2

	improve_poss = True
	time = 1
	while improve_poss:
		old_nodes = nodes[::]
		##randomness is introduced for simulations (otherwise every run will be identical)
		nodes = random.sample(nodes[::],len(nodes))
		matrix = matrix_changer(nodes, old_nodes, matrix)
		comms = {}
		for node in nodes:
			comms[node] = [node]		
		old_mod = modularity(nodes, comms, neighbors, matrix, total_weight)
		orig_comms = copy.deepcopy(comms)
		##finds the max change in modularity from moving a node into the clusters of its neighbors
		for node in nodes:
			max_change = 0
			move = []
			for neigh in neighbors[node]:
				orig_comms_node = comms[node]
				orig_comms_neigh = comms[neigh]
				comms[node] = comms[node] + orig_comms_neigh
				comms[neigh] = comms[neigh] + orig_comms_node
				change = mod_change(node, comms, edge_weights, total_weight)
				if change>max_change:
					max_change = change
					move = neigh
				comms[node] = orig_comms_node
				comms[neigh] = orig_comms_neigh
			##the reclustering that led to the highest increase in modularity (if it is greater)
			##than zero) is implemented
			if move:
				#removes node from current cluster
				for x in comms[node]:
					if x!=node:
						comms[x].remove(node)
				#adds node to new cluster
				orig_comms_move = comms[move][::]
				comms[node] = [node] + comms[move]
				for x in orig_comms_move:
					comms[x].append(node)
		##calculates the modularity of the new clustering
		new_mod = modularity(nodes, comms, neighbors, matrix, total_weight)
		##if the new modularity is greater than the old, the new clustering is used; else,
		##the function returns the old clustering as there it decides that it has found the
		##optimum clustering where no higher modularity can be reached (though this can not be true)
		if new_mod>old_mod:
			new_nodes = []
			checked_nodes = []
			for comm_node in comms:
				comm = comms[comm_node]
				do = 'DO_NOT_ADD'
				for node in comm:
					if node in checked_nodes:
						break
					else:
						checked_nodes.append(node)
						do = 'DO_ADD'
				if do=='DO_ADD':
					if len(comm)==1 and type(comm[0])==tuple:
						new_nodes.append(comm[0])
					else:
						new_nodes.append(tuple(comm))
			time+=1
			nodes = new_nodes
			neighbors, edge_weights, matrix = new_neighbors_eweights_and_matrix(nodes, edge_weights)
		else:
			return old_mod, orig_comms

##finds the new neighbor and edge weight dictionaries as well as the adjacency matrix of a network
##after a new clustering
def new_neighbors_eweights_and_matrix(nodes, edge_weights):
	matrix = np.zeros((len(nodes),len(nodes)))
	new_edge_weights = {}
	neighbors = {}
	for tuple1 in nodes:
		neighbors[tuple1] = []
		for tuple2 in nodes:
			new_edge_weights[tuple([tuple1,tuple2])] = 0

	checked = []
	for edge in edge_weights:
		if edge_weights[edge]!=0:
			for i in range(len(nodes)):
				tuple1 = nodes[i]
				for j in range(len(nodes)):
					tuple2 = nodes[j]
					truth = True
					for e1 in edge:
						for e2 in edge:
							if e1!=e2 and truth==True:
								if (e1 in tuple1 and e2 in tuple2) or (e2 in tuple1 and e1 in tuple2):
									w = float(edge_weights[edge])/2
									new_tup = tuple([tuple1,tuple2])
									new_edge_weights[new_tup] = new_edge_weights[new_tup] + w
									matrix[i][j]+=w
									if tuple1!=tuple2 and [tuple1,tuple2] not in checked:
										neighbors[tuple1] = neighbors[tuple1] + [tuple2]
										neighbors[tuple2] = neighbors[tuple2] + [tuple1]
										checked.append([tuple1,tuple2])
										checked.append([tuple2,tuple1])
									truth = False
	return neighbors, new_edge_weights, matrix

##creates a dictionary of edge weights called edge_weights based on an adjacency matrix
def edge_weight_finder(nodes, matrix):
	edge_weights = {}
	for i in range(len(nodes)):
		for j in range(len(nodes)):
			if matrix[i][j]!=0:
				edge_weights[(nodes[i],nodes[j])] = matrix[i][j]
	return edge_weights

##computes the modularity of a clustering of nodes
def modularity(nodes, comms, neighbors, adj_mat, total_weight):
	mod = 0
	i = 0
	for node1 in nodes:
		for node2 in comms[node1]:
			j = nodes.index(node2)
			mod = mod + adj_mat[i][j] - (sum(adj_mat[i])*sum(adj_mat[j]))/(2*total_weight)
		i+=1
	mod = mod/(2*total_weight)
	return mod

##computes shortcut formula for calculating the change in modularity when one node is moved
def mod_change(node, comms, edge_weights, total_weight):
	c = comms[node]
	c_total_sum = 0.0
	c_in_sum = 0.0
	k_i = 0.0
	k_i_in = 0.0
	for edge in edge_weights:
		w = float(edge_weights[edge])/2
		if edge[0] in c or edge[1] in c:
			c_total_sum+=w
		if edge[0] in c and edge[1] in c:
			c_in_sum+=w
		if edge[0]==node or edge[1]==node:
			k_i+=w
		if (edge[0]==node and edge[1] in c) or (edge[1]==node and edge[0] in c):
			k_i_in+=w
	change = (1.0/4.0*(total_weight**2))*((2.0*total_weight*(c_in_sum+k_i_in)) - (c_total_sum+k_i)**2 - (2.0*total_weight*c_in_sum) + c_total_sum**2 + k_i**2)
	return change

if __name__=='__main__':                                                         
	main()