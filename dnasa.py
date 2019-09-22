import itertools
import networkx as nx
import sympy
import numpy as np

def partner_indeces(target, half_edge_list):
	indeces = []
	for i, half_edge_tuple in enumerate(half_edge_list):
		if half_edge_tuple[0] == target:
			indeces.append(i)
	return indeces

def graph_construction_generator(half_edge_list):
	if len(half_edge_list) == 2:
		yield [(half_edge_list[0][1], half_edge_list[1][1])]
	else:
		partner_indeces_ = partner_indeces(-1*half_edge_list[0][0],
											half_edge_list)
		for partner_index in partner_indeces_:
			this_half_edge_list = half_edge_list.copy()
			this_half_edge_list.pop(partner_index)
			this_half_edge_list.pop(0)
			for graph_construction in graph_construction_generator(
															this_half_edge_list):
				yield [(half_edge_list[0][1], half_edge_list[partner_index][1])] + graph_construction

def tile_degree_quantities(graph):
	tile_degree_quantities = {}
	adjacency_matrix = nx.to_numpy_matrix(graph)
	for i in range(np.shape(adjacency_matrix)[0]):
		degree = 0
		for j in range(np.shape(adjacency_matrix)[1]):
			if adjacency_matrix[i,j] == 1:
				degree += 1
		if degree in tile_degree_quantities.keys():
			tile_degree_quantities[degree] += 1
		else:
			tile_degree_quantities[degree] = 1
	return tile_degree_quantities

def can_create_complete_graph(tiles):
	if type(tiles[0]) == list:
		half_edges = []
		for tile in tiles:
			half_edges = half_edges + tile
	else:#in the case of tiles being only a single tile
		half_edges = tiles
		
	if (len(half_edges) % 2 != 0):#easy case
		return False
	
	#match up each half-edge with a bond partner
	used_half_edge = [False] * len(half_edges)
	for i, half_edge in enumerate(half_edges):
		#if this half-edge already has a partner, move to the next half-edge
		if used_half_edge[i] == True:
			continue
		else:
			used_half_edge[i] = True
			foundPartner = False
			#search for a partner for this half-edge
			for j in range(i + 1, len(half_edges)):
				if (half_edge + half_edges[j] == 0
						and (used_half_edge[j] == False)):
					used_half_edge[j] = True
					foundPartner = True
					break
			if foundPartner == False:
				return False
	return True

def pot_expansions(pot, graph_order):
	if not (len(pot) < graph_order):
		return pot
	else:
		pot_expansions = []
		partitions = [integer_partition_generator(graph_order, len(pot))]
		for partition in partitions:
			for partition_permutation in sympy.utilities.iterables.multiset_permutations(partition):
				expansion = []
				for i, summand in enumerate(partition_permutation):
					for j in range(summand):
						expansion.append(pot[i])
				pot_expansions.append(expansion)
		return pot_expansions

				
				
		

def can_create_smaller_graph(pot_expansion):
	for i in range(1, len(pot_expansion)):
		for tile_combination in itertools.combinations(pot_expansion, i):
			if can_create_complete_graph(tile_combination):
				return True
	return False

def creates_nonisomorphisms(pot_expansion, graph):
	all_constructions_isomorphic = True
	half_edge_list = []
	for i, tile in enumerate(pot_expansion):
		for half_edge in tile:
			half_edge_list.append((half_edge, i))
	for adjacency_list in graph_construction_generator(half_edge_list):
		g = nx.Graph()
		g.add_edges_from(adjacency_list)
		if not nx.is_isomorphic(g, graph):
			all_constructions_isomorphic = False
			break
	if all_constructions_isomorphic:
		return False
	else:
		return True

def can_create_target(pot_expansion, graph):
	graph_tile_degree_quantities = tile_degree_quantities(graph)
	pot_tile_degree_quantities = {}
	for tile in pot_expansion:
		if len(tile) in pot_tile_degree_quantities.keys():
			pot_tile_degree_quantities[len(tile)] += 1
		else:
			pot_tile_degree_quantities[len(tile)] = 1
	if (graph_tile_degree_quantities == pot_tile_degree_quantities) and can_create_complete_graph(pot_expansion):
		return True
	else:
		return False	

'''This function implements an iterative solution for generating all integer
partitions of some non-negative target integer, n. Not only that, but the number
of summands in the resultant partitions may be arbitrarily limited using the
max_summands argument. The resultant generator generates lists of length
max_summands, wherein each entry refers to a summand in the partition. The
argument max_summands must be a postive integer less than or equal to n.'''
def integer_partition_generator(n,max_summands=None):
	assert type(n) is int
	assert n >= 0
	if max_summands == None:
		max_summands = n
	else:
		assert max_summands > 0
		assert max_summands <= n
	partition = [n] + [0]*(max_summands-1)
	while True:
		yield partition
		evenly_distributed = True
		if str(n/max_summands)[1:] == '0':
			for summand in partition:
				if summand != (n//max_summands):
					evenly_distributed = False
					break
		else:
			for summand in partition:
				if summand != (n//max_summands) and summand != ((n//max_summands)+1):
					evenly_distributed = False
					break
		if evenly_distributed == True:
			break
		try:
			end = partition.index(0)
		except ValueError:
			end = 1
		for i in reversed(range(end)):
			if partition[i] != 1:
				partition[i] -= 1
				anchor = sum(partition[0:i+1])
				for j in range(i+1,len(partition)):
					partition[j]=0
				distribute = n-anchor
				j = i+1
				while distribute > 0:
					if j>=len(partition):
						evenly_distributed = True
						break
					distribution = min(partition[j-1],distribute)
					partition[j] = distribution
					distribute -= distribution
					j += 1
				break
		if evenly_distributed == True:
			break
		
def scenario2(pot, graph, scenario3=False):
	pot_expansions_ = pot_expansions(pot, graph.order())
	if not scenario3:
		one_can_create_graph = False
		for pot_expansion in pot_expansions_:
			if can_create_smaller_graph(pot_expansion):
				return False
			if can_create_target(pot_expansion, graph):
				one_can_create_graph = True
		if one_can_create_graph:
			return True
		else:
			return False
	else:
		pot_expansions_to_be_checked = []
		one_can_create_graph = False
		for pot_expansion in pot_expansions_:
			if can_create_smaller_graph(pot_expansion):
				return (False, None)
			if can_create_complete_graph(pot_expansion):
				pot_expansions_to_be_checked.append(pot_expansion)
			if can_create_target(pot_expansion, graph):
				one_can_create_graph = True
		if one_can_create_graph:
			return (True, pot_expansions_to_be_checked)
		else:
			return (False, None)

def scenario3(pot, graph):
	s2result = scenario2(pot, graph)
	if s2result[0] == False:
		return False
	else:
		pot_expansions = s2result[1]
		for pot_expansion in pot_expansions:
			if creates_nonisomorphisms(pot_expansion):
				return False
		return True