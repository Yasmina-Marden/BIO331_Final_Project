##Opens a text file and then saves the graph info inside as a nodes list and two 
##dictionaries for edges and neighbors. The edge dictionary is created based on edge_type = True 
##or False if the graph in the text file is directed or undirected, respectively. Additional
##information about the edges such as weights is stored in the edges dictionary in a list
##with the target node.
def readData(textFile, edge_type):
	text_file = open(textFile, 'r')
	text = text_file.read()
	text_file.close()
	text = text.split('\n')

	for i in range(len(text)):
		text[i]=text[i].split('\t')

	nodes = []
	edges = {}
	used_nodes = []
	neighbors = {}
	for i in range(len(text)):
		u = text[i][0]
		v = text[i][1]
		new = [v]+text[i][2:]
		if str(u) not in edges:
			edges[u] = [new]
		elif new not in edges[u]:
			edges[u] = edges[u]+[new]

		new = [u]+text[i][2:]
		if edge_type==False:
			if str(v) not in edges:
				edges[v] = [new]
			elif new not in edges[v]:
				edges[v] = edges[v]+[new]

		if str(u) not in neighbors:
			neighbors[u]=[v]
		elif v not in neighbors[u] and v not in neighbors[u]:
			neighbors[u] = neighbors[u]+[v]
		if str(v) not in neighbors:
			neighbors[v]=[u]
		elif u not in neighbors[v] and u not in neighbors[v]:
			neighbors[v] = neighbors[v]+[u]
		if u not in used_nodes:
			used_nodes.append(u)
		if v not in used_nodes:
			used_nodes.append(v)
	for node in nodes:
		if node not in used_nodes:
			neighbors[node] = []

	nodes = []
	for x in neighbors.keys():
		nodes.append(x)

	return nodes, edges, neighbors

def readMatrix(textFile):
	D = []
	headers = []
	with open(textFile) as fin:
		for line in fin:
			row = line.strip().split()
			if len(headers)==0:
				headers = row
			else: 
				row = [int(a) for a in row[1:]]
				D.append(row)
	return D, headers
