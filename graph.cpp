#include "graph.h"
#include "io.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <set>

// Constructor

Graph::Graph() {}

Graph::Graph(char *filename_graph, char *filename_costs, char *filename_lists) {
	read_graph(filename_graph);
	read_costs(filename_costs);
	read_lists(filename_lists);
}


Graph::~Graph() {

#if BRANCHING_STRATEGY == 0
	free_graph(G);
	free(G);
#endif

	for (int i = 0; i < K.size(); ++i) 
		delete[] V[i].list;

#ifdef STABLE_POOL
	for (int i = 0; i < K.size(); ++i) {
		for (auto &s: global_pool[i]) 
			free(s.list);
		for (auto &s: local_pool[i]) 
			free(s.list);
	}
#endif

}


// read_graph - read a graph in the following format:
//   in the first line, the number of vertices and edges separated by a colon ":"
//   then, for each line, the endpoints of an edge (u,v) where u < v
//   note that vertices starts from 0, i.e. 0 < v < |V|-1
//   example for a diamond graph:
//     4:5
//     0,1
//     0,2
//     0,3
//     1,2
//     1,3
void Graph::read_graph(char *filename) {

	// Open file 
	FILE *stream = fopen(filename, "rt");
	if (!stream) bye("Graph file cannot be opened");

	// Read number of vertices and edges
    int vertices, edges;
	fscanf(stream, "%d:%d\n", &vertices, &edges);

	// Do not accept graph of less than 4 vertices or stable sets
	if (vertices < 4) bye("Number of vertices out range!");
	if (edges < 1 || edges > vertices*(vertices - 1) / 2) bye("Number of edges out of range!");

	// Allocate graph G
	G = (MWSSgraphpnt) malloc(sizeof(MWSSgraph));
	allocate_graph(G, vertices);
	for(int i = 0; i <= G->n_nodes; i++) {
		for(int j = 0; j <= G->n_nodes; j++) {
			G->adj[i][j] = 0;
		}
	}

	// Read edges
	for (int e = 0; e < edges; e++) {
		int u, v;
		fscanf(stream, "%d,%d\n", &u, &v);
		if (u < 0 || u >= v || v >= vertices) bye("Error while reading edges");
		G->adj[u+1][v+1] = 1;
		G->adj[v+1][u+1] = 1;
	}

	fclose(stream);

	// Build G
	build_graph(G);

	// Set some information necessary to solve the MWSSP
	for(int i = 1; i <= G->n_nodes; i++) {
		MWIS_MALLOC(G->node_list[i].adj_last, G->n_nodes+1, nodepnt*);
		G->node_list[i].adj_last[0] = G->adj_last[i];
		G->node_list[i].adj2 = G->adj_last[i];
	}

#if BRANCHING_STRATEGY == 0
	vertex_mapping.resize(get_n_vertices() + 1);
	for (int i = 1; i <= get_n_vertices(); ++i)
		vertex_mapping[i] = i;
#endif

    return;

}

// read_cost - read costs of colors in the following format:
//   in the first line, the number of colors
//   then, the costs of color 0, 1, etc, in succesive order
//   example for C = {0, 1, 2} with costs c_0 = 5, c_1 = 2, c_2 = 8:
//     3
//     5 2 8
void Graph::read_costs(char *filename) {

	// Open file
	FILE *stream = fopen(filename, "rt");
	if (!stream) bye("Cost file cannot be opened");
    int colors;
	fscanf(stream, "%d\n", &colors);

	// Do not accept less than 2 colors
	if (colors < 2) bye("Number of colors out range!");

    w.resize(colors);

	// Read costs
	for (int k = 0; k < colors; k++) {
		int ck;
		fscanf(stream, "%d", &ck);
		if (ck < 0) bye("Color cost must be non negative!");
		w[k] = ck;
	}
	fclose(stream);

    return;

}


// read_list - read list of colors per vertex in the following format:
//   in the first line, the number of vertices and colors separated by a colon ":"
//   then, for each line, the cardinal of L(v) followed by the elements of L(v) in increasing order
//   example for |V| = 3, |C| = 5 and L(0) = {1, 2}, L(1) = {0, 2, 3}, L(2) = {0, 1, 4}:
//     3:5
//     2  1 2
//     3  0 2 3
//     3  0 1 4
void Graph::read_lists(char *filename) {

	// Open file
	FILE *stream = fopen(filename, "rt");
	if (!stream) bye("Cost file cannot be opened");

	// Read number of vertices and colors
	int vertices, colors;
	fscanf(stream, "%d:%d\n", &vertices, &colors);
	if (G->n_nodes != vertices) bye("Number of vertices mismatch!");
	if (w.size() != colors) bye("Number of colors mismatch!");

	// Read lists
	std::vector<std::vector<int>> L (vertices);
	for (int v = 0; v < vertices; v++) {

		int list_size;
		fscanf(stream, "%d", &list_size);
		if (list_size < 1 || list_size > colors) bye("Error reading lists!");
		L[v].resize(list_size);

		int last_read = -1;
		for (int s = 0; s < list_size; s++) {
			int element;
			fscanf(stream, "%d", &element);
			if (element <= last_read || element >= colors) bye("Error reading lists!");
			last_read = element;
			L[v][s] = element;
        }

	}
	fclose(stream);

	// Build subgraphs
	std::vector<std::set<int>> total_V (colors); // V[k] for each 0 <= k < colors
	for (int k = 0; k < colors; k++) {
		// find vertices v such that k in L(v)
		for (int v = 0; v < vertices; v++)
			if (find(L[v].begin(), L[v].end(), k) != L[v].end())
				total_V[k].insert(v+1);
		if (total_V[k].size() == 0) warning("Warning color " + std::to_string(k) + " does not belong to any list.");
	}


    // Build partition
	// The structure of the partition consists of an array and a linked-list. The linked-list is formed with
	// header_part (points to the first) and next_part (points to the next).
	int *header_part = new int[colors]; 	// is the first color (representative) of a given color
	int *next_part = new int[colors];   	// is the next color in the set containing the given color, -1 if the given color is the last
	int *part_set = new int[colors]; 		// array with the first color of each set in the partition
	int *part_card = new int[colors]; 		// array with the cardinal of each set
	int part_size = 0;  		// size of the array

	for (int k = 0; k < colors; k++) {
		bool new_color = true;
		int prev_index; // index to the set of the partition in case color k is not new
		if (k > 0) {
			/* check if the color k is the same as a previous color of the partition */
			for (int i = 0; i < part_size; i++) {
				prev_index = i;
				int prev = part_set[i];
				if (w[k] == w[prev] && total_V[k] == total_V[prev])
					new_color = false;
				if (new_color == false) break;
			}
		}
		if (new_color) {
			// a new color k is added to the partition
			header_part[k] = k;
			next_part[k] = -1;
			part_set[part_size] = k;
			part_card[part_size] = 1;
			part_size++;
		}
		else {
			/* colors k and part_set[prev_index] are indistinguishable from each other  */
			int prev = part_set[prev_index];
			part_card[prev_index]++;
			header_part[k] = prev;
			next_part[k] = -1;
			/* travel the set until reach the last element and append k */
			int r = prev;
			int s;
			do { s = r; r = next_part[s]; } while (r != -1);
			next_part[s] = k;
		}
	}

	// Save partition

	K.resize(part_size);
	for (int i = 0; i < part_size; ++i)
		K[i] = part_set[i];

	C.resize(K.size());
	for (int i = 0; i < K.size(); ++i) {
		C[i].resize(part_card[i]);
		int counter = 0;
		for (int j = part_set[i]; j != -1; j = next_part[j])
			C[i][counter++] = j;
	}

	V.resize(K.size());
	for (int i = 0; i < K.size(); ++i) {
		V[i].n_list = total_V[K[i]].size();
		V[i].list = new nodepnt[V[i].n_list + 1];
		int counter = 1;
		for (int v: total_V[K[i]])
			V[i].list[counter++] = G->node_list + v;
	}

#ifdef STABLE_POOL
	global_pool.resize(K.size());
	local_pool.resize(K.size());
#endif

	delete[] part_card;
	delete[] part_set;
	delete[] next_part;
	delete[] header_part;
	return;

}


int Graph::get_n_vertices() {
    return G->n_nodes;
}


int Graph::get_n_colors() {
    return K.size();
}

int Graph::get_n_V(int i) {
	return V[i].n_list;
}

int Graph::get_C(int i, int j) {
	return C[i][j];
}

int Graph::get_n_C(int i) {
    return C[i].size();
}

int Graph::get_color_cost(int i) {
	return w[K[i]];
}

int Graph::get_m(int v) {
	int max = -1;
	for (int i = 0; i < K.size(); ++i) {
		bool belongs = false;
		for (int j = 1; j <= V[i].n_list; ++j)
			if (v+1 == V[i].list[j]->name) {
				belongs = true;
				break;
			}
		if (belongs) {
			int cost = get_color_cost(i);
			if (cost > max)
				max = cost;
		}
	}
	return max;
}

int Graph::get_n_neighbours(int v, int k) {
	int ret = 0;
	for (int i = 1; i <= V[k].n_list; ++i) {
		int u = V[k].list[i]->name;
		if (u != v+1 && G->adj[u][v+1])
			ret++;
	}
	return ret;
}

void Graph::get_W1(std::vector<int> &W) {
	W.clear();
	W.resize(get_n_vertices(), 0);
	for (int i = 0; i < K.size(); ++i)
		for (int j = 1; j <= V[i].n_list; ++j)
			W[V[i].list[j]->name - 1]++;
}

bool Graph::is_edge(int u, int v) {
	return G->adj[u+1][v+1];
}

void Graph::set_vertex_weight(int v, double y) {
	G->weight[v+1] = (y <= 0 ? 0 : (int) (INTFACTOR * y));
	return;
}

void Graph::print_graph() {
#ifdef VERBOSE
	
	// Print Graph
	std::cout << "Graph:" << std::endl;
	std::cout << "\tVertices: " << G->n_nodes << std::endl;
	std::cout << "\tEdges: " << G->n_edges << std::endl;
	std::cout << "\tAdjacency list: " << std::endl;
	for(int i = 1; i <= G->n_nodes; i++) {
		printf("%5d:",i);
		int cnt = 0;
		nodepnt *pnt = G->node_list[i].adj;
		while(*pnt != NULL) {
			cnt++;
			nodepnt pnt2 = *pnt;
			printf("%5d%s",pnt2->name, (*(++pnt) != NULL) && (cnt % 14) == 0 ? "\n          ":"");
		}
		printf("\n");
	}
	std::cout << std::endl;

	// Print costs
	std::cout << "Costs:" << std::endl;
	for (int i = 0; i < w.size(); ++i)
		std::cout << "\t" << i << ":\t" << w[i] << std::endl;
	std::cout << std::endl;

	// Print colors up to indistinguishability
	std::cout << "Colors:" << std::endl;
	for (int i = 0; i < K.size(); ++i)
		std::cout << "\t" << K[i] << std::endl;
	std::cout << std::endl;

	// Classes of equivalence
	std::cout << "Classes of equivalence:" << std::endl;
	for (int i = 0; i < K.size(); ++i) {
		std::cout << "\t" << K[i] << ":\t";
		for (int k: C[i])
			std::cout << k << "\t";
		std::cout << std::endl;
	}
	std::cout << std::endl;

	// Subgraphs
	std::cout << "Subgraphs:" << std::endl;
	for (int i = 0; i < K.size(); ++i) {
		std::cout << "\t" << K[i] << ":";
		for (int j = 1; j <= V[i].n_list; j++) {
			printf("%5d%s", V[i].list[j]->name, (j % 16) == 0 ? "\n":"");
		}
		printf("\n");
	}
	std::cout << std::endl;

#ifdef STABLE_POOL
	std::cout << "Pool:" << std::endl;
	// Print stable set's global pool
	for (int i = 0; i < K.size(); ++i)
		for (auto &s: global_pool[i]) {
			std::cout << "\t" << K[i] << ":";
			for (int j = 1; j <= s.n_list; j++) {
				printf("%5d%s", s.list[j]->name, (j % 16) == 0 ? "\n":"");
			}
			printf("\n");
		}
#endif

#endif
	return;
}

// preprocess_instance - preprocess the instance:
// For all vertex v such that k(v) = 1, v \in V_j, m(j) >= 2 and j' \in C_j\{j}
// K'= K\cup \{j'\}
// V'_k= V_k, w'(k)=w(k), and m'(k)=m(k), for all k \in K\{j}
// V'_j = V_j\N_{G_j}(v), w'(j)=w(j), and m'(j)=1
// V_{j'}=V_j\{v}, w(j')=w(j), and m(j')=m(j)-1
void Graph::preprocess_instance() {
	
	// Define L[v] = {j \in K: v \in V_j}
	std::vector<std::list<int>> L(get_n_vertices() + 1);
	for (int j = 0; j < K.size(); ++j)
		for (int v = 1; v <= V[j].n_list; ++v)
			L[V[j].list[v]->name].push_back(j);
	
	// Define the list of candidates (v,j) with k(v) = 1, v \in V_j, m(j) >= 2
	std::list<std::pair<int,int>> candidates;
	for (int v = 1; v <= get_n_vertices(); ++v)
		if (L[v].size() == 1 && get_n_C(L[v].front()) >= 2)
			candidates.push_back(std::make_pair(v,L[v].front()));
	
	// Find the candidate (v,j) with the greatest N_{G_j}(v)
	int max_neighbors = -1;
	std::pair<int,int> max_vj (-1,-1);
	for (auto &p: candidates) {
		int n = get_n_neighbours(p.first - 1, p.second);
		if (n > max_neighbors) {
			max_neighbors = n;
			max_vj = p;
		}
	}
	
	// If there is not any candidate, return
	if (max_neighbors == -1)
		return;
	
	// Preprocess the instance
	
	int v = max_vj.first;
	int j = max_vj.second;

	// Build vector of indistinguishable colors
	K.resize(K.size() + 1);
	K[j] = C[j][1];
	K[K.size()-1] = C[j][0];

	// Build the partition into indistinguishable colors
	C.resize(C.size() + 1);
	C[j].assign(C[j].begin() + 1, C[j].end());
	C[C.size()-1].push_back(C[j][0]);

	// Build Vk's
	V.resize(V.size() + 1);
	// Remove N_{G_j}(v) from V[V.size() - 1]
	nodepnt *V1 = new nodepnt[V[j].n_list - max_neighbors + 1];
	int index = 1;
	for (int i = 1; i <= V[j].n_list; ++i) {
		int u = V[j].list[i]->name;
		if (u == v || (u != v && !G->adj[u][v]))
			V1[index++] = V[j].list[i];
	}
	V[V.size() - 1].n_list = V[j].n_list - max_neighbors;
	V[V.size() - 1].list = V1;
	// Remove vertex v from V[j]
	nodepnt *V2 = new nodepnt[V[j].n_list];
	index = 1;
	for (int i = 1; i <= V[j].n_list; ++i)
		if (V[j].list[i]->name != v)
			V2[index++] = V[j].list[i];
	V[j].n_list--; 
	delete[] V[j].list;
	V[j].list = V2;

#ifdef STABLE_POOL
	bye("Stable pool is not implemented yet when preprocessing is on.\n");
#endif	
	
	// Recursion
	preprocess_instance();

	return;

}

bool Graph::solve_MWSSP(int i, double goal, nodepnt **best_stable, int *n_best_stable) {

	// First, we perform a quick check to determine if the answer is negative
	// 	Does the total weight of the vertices of V[i] is lower or equal than the goal?
	int total_weight = 0;
	for (int j = 1; j <= V[i].n_list; ++j)
		total_weight += G->weight[V[i].list[j]->name];
	if (total_weight <= (int) INTFACTOR * goal)
		return false;

#ifdef STABLE_POOL
	// We try to solve the MWSSP heuristically
	if (solve_MWSSP_heuristic(i, goal, best_stable, n_best_stable)) {
		local_pool[i].emplace_back(*best_stable, *n_best_stable);
		return true;
	}
	// Otherwise, proceed with the exact algorithm
	else if (solve_MWSSP_exact(i, goal, best_stable, n_best_stable)) {
		maximize_stable_set(i, best_stable, n_best_stable);	
		local_pool[i].emplace_back(*best_stable, *n_best_stable);
		return true;
	}
	else
		return false;
#else
	// Apply the the exact algorithm
	if (solve_MWSSP_exact(i, goal, best_stable, n_best_stable)) {
		maximize_stable_set(i, best_stable, n_best_stable);
//		if (!check_stable(i, *best_stable, *n_best_stable))
//		 	bye ("Te devolvi algo que no es un estable");
		return true;
	}
	else {
		return false;
	}
#endif

}


bool Graph::solve_MWSSP_heuristic(int i, double goal, nodepnt **best_stable, int *n_best_stable) {

	// Try to find an stable set in the global pool with a total weight greater than goal
	for (auto it = global_pool[i].begin(); it != global_pool[i].end(); ++it) {
		int total_weight = 0;
		for (int j = 1; j <= it->n_list; ++j)
			total_weight += G->weight[it->list[j]->name];
		if (total_weight > (int) INTFACTOR * goal) {
			*best_stable = it->list;
			*n_best_stable = it->n_list;
			it = global_pool[i].erase(it);
			return true;
		}
	}

	return false;

}


bool Graph::solve_MWSSP_exact(int i, double goal, nodepnt **best_stable, int *n_best_stable) {

	int rval = 0;
	int status;
	MWISNW z_best;

	// Initialize info
	wstable_info info;
	info.n_calls = 0;
	info.cpu = 0;
	info.clique_cover_cpu = 0;
	MWIS_MALLOC(info.n_sub_depth, G->n_nodes + 1, int);

	// Initialize data
	MWSSdata data;
	rval = allocate_data(&data, G->n_nodes);

	MWIS_MALLOC(*best_stable, V[i].n_list + 1, nodepnt); // Reserve as much space as |V[i]| + 1

	// Initialize best_stable
	// Initialize parameters
	wstable_parameters parms;
	default_parameters(&parms);
	parms.cpu_limit = MAXTIME_MWIS;

	// Find MWSSP
	rval = max_wstable(G, &data, *best_stable, n_best_stable, &z_best, &info,
						&parms, V[i].list, V[i].n_list, 0, (int) INTFACTOR * (goal + THRESHOLD), &status);

	// Handle time exceeded error
	if (parms.cpu_limit >= 0 && info.cpu > parms.cpu_limit) bye("Error: Time exceeded in MWSSP");	

	// Free info and data
	free_data(&data);
	free_wstable_info(&info);

	// Return
	if (z_best > (int) INTFACTOR * goal) {
		return true;
	}
	else {
		free(*best_stable);
		return false;
	}
}


void Graph::maximize_stable_set(int i, nodepnt **best_stable, int *n_best_stable) {

	for (int j = 1; j <= V[i].n_list; ++j) {
		int v = V[i].list[j]->name;
		bool is_neighbour = false;
		for (int k = 1; k <= *n_best_stable; ++k) {
			int u = (*best_stable)[k]->name;
			if (u == v || G->adj[u][v]) {
				is_neighbour = true;
				break;
			}
		}
		if (!is_neighbour) {
			(*n_best_stable)++;
			(*best_stable)[*n_best_stable] = V[i].list[j];
		}
	}

}


void Graph::coloring_heuristic(std::vector<std::list<nodepntArray>>& stable_sets) {

	stable_sets.resize(get_n_colors());

	// Mark the vertices already colored
    std::vector<bool> colored (get_n_vertices(), false);
	int n_colored = 0;

	// Sort the colors by weight
	std::vector<int> sorted_colors (get_n_colors());
	for (int i = 0; i < get_n_colors(); ++i)
		sorted_colors[i] = i;
	std::sort(sorted_colors.begin(), sorted_colors.end(), [this](int x, int y) {return (get_color_cost(x) < get_color_cost(y));} );

    for (int color_counter = 0; color_counter < get_n_colors(); ++color_counter) {
		int i = sorted_colors[color_counter];
		for (int multiplicity = 0; multiplicity < get_n_C(i); ++multiplicity) {
			nodepnt *stable;
			int n_stable = 0;
			MWIS_MALLOC(stable, V[i].n_list + 1, nodepnt);
			for (int j = 1; j <= V[i].n_list; ++j) {
				int v = V[i].list[j]->name;
                if (colored[v]) continue;

                // Try to add v to the stable
                bool ok = true;
				for (int k = 1; k <= n_stable; ++k){
					int u = stable[k]->name;
					if (u == v || G->adj[u][v]) {
						ok = false;
						break;
					}
				}
                if (ok) {
                    stable[++n_stable] = V[i].list[j];
                    colored[v] = true;
                    n_colored++;
					if (n_colored == get_n_vertices()) break;
                }

			}

            if (n_stable > 0) {
                maximize_stable_set(i, &stable, &n_stable);
                stable_sets[i].emplace_back(stable, n_stable);
            }
			else {
				free(stable);
			}

			if (n_colored == get_n_vertices())
				return;

		}

    }

    return;

}


bool Graph::check_stable(int color, nodepnt *stable, int n_stable) {

	for (int i = 1; i <= n_stable; ++i) {
		bool belongs = false;
		for (int j = 1; j <= V[color].n_list; ++j)
			if (V[color].list[j]->name == stable[i]->name) {
				belongs = true;
				break;
			}
		if (!belongs)
			return false;
	}

	for (int i = 1; i <= n_stable-1; ++i)
		for (int j = i+1; j <= n_stable; ++j) {
			int u = stable[i]->name;
			int v = stable[j]->name;
			if (G->adj[u][v] || G->adj[v][u])
				return false;
		}
	
	return true;
}


bool Graph::check_coloring(std::vector<int> &) {
	// TODO
	return true;
}


void Graph::update_pool() {
	for (int i = 0; i < K.size(); ++i) {
		// Delete the current global_pool
		for (auto &s: global_pool[i]) 
			free(s.list);
		global_pool[i].clear();
		// Copy the local pool to the global pool
		for (auto it = local_pool[i].begin(); it != local_pool[i].end();) {
			global_pool[i].emplace_back(it->list, it->n_list);
			it = local_pool[i].erase(it);
		}
	}
}

void Graph::column_to_stable (bool *column, int size, nodepnt **stable) {
	MWIS_MALLOC(*stable, size + 1, nodepnt);
	int index = 1;
	for (int j = 0; j < get_n_vertices(); ++j)
		if (column[j])
			(*stable)[index++] = G->node_list + (j+1);
	return;
}

void Graph::stable_to_column (nodepnt *stable, int size, bool **column) {
	*column = new bool[get_n_vertices()];
	for (int j = 0; j < get_n_vertices(); ++j)
		(*column)[j] = false;
	for (int j = 1; j <= size; ++j)
		(*column)[stable[j]->name - 1] = true;
	return;
}


void Graph::translate_stable_set (int color, nodepnt *stable_father, int n_stable_father, nodepnt **stable_son, int *n_stable_son) {

#if BRANCHING_STRATEGY == 0

	// Allocate memory for the son's stable set
	MWIS_MALLOC(*stable_son, V[color].n_list + 1, nodepnt);
	*n_stable_son = 0;

	if (st == JOIN) {
		bool is_u = false, is_v = false;
		for (int j = 1; j <= n_stable_father; ++j) {
			int w = stable_father[j]->name;
			if (w == branch_vertex_u) {
				is_u = true;
				if (is_v) continue;
			}
			else if (w == branch_vertex_v) {
				is_v = true;
				if (is_u) continue;
			}
			(*n_stable_son)++;
			(*stable_son)[*n_stable_son] = G->node_list + w;
		}
		if (is_u && is_v)
			maximize_stable_set(color, stable_son, n_stable_son);

		return;	
	}

	else if (st == COLLAPSE) {
		bool is_u = false, is_v = false;
		for (int j = 1; j <= n_stable_father; ++j) {
			int w = stable_father[j]->name;
			if (w < branch_vertex_u) {
				(*n_stable_son)++;
				(*stable_son)[*n_stable_son] = G->node_list + w;				
			}
			else if (w == branch_vertex_u) {
				is_u = true;
				continue;
			}
			else if (w < branch_vertex_v) {
				(*n_stable_son)++;
				(*stable_son)[*n_stable_son] = G->node_list + w;
			}
			else if (w == branch_vertex_v) {
				is_v = true;
				continue;
			}
			else {
				(*n_stable_son)++;
				(*stable_son)[*n_stable_son] = G->node_list + (w-1);
			}
		}
		if (is_u && is_v) {
			(*n_stable_son)++;
			(*stable_son)[*n_stable_son] = G->node_list + branch_vertex_u;
		}
		else if (is_u != is_v)
			maximize_stable_set(color, stable_son, n_stable_son);

		return;	
	}
	
	else
		bye ("Undefined branching status");

#else
	bye ("Undefined translate_stable_set() for the current branching strategy");
#endif

}

Graph *Graph::join_vertices(int u, int v) {

	// Build a new graph
	Graph *H = new Graph();

	// Vertices start at 1
	u = u+1;
	v = v+1;

	// Build Sewell's graph
	H->G = (MWSSgraphpnt) malloc(sizeof(MWSSgraph));
	allocate_graph(H->G, get_n_vertices());
	for(int i = 0; i <= H->G->n_nodes; i++)
		for(int j = 0; j <= H->G->n_nodes; j++)
			H->G->adj[i][j] = G->adj[i][j];
	H->G->adj[u][v] = 1;
	H->G->adj[v][u] = 1;

	build_graph(H->G);

	// Set some information necessary to solve the MWSSP
	for(int i = 1; i <= H->G->n_nodes; i++) {
		MWIS_MALLOC(H->G->node_list[i].adj_last, H->G->n_nodes+1, nodepnt*);
		H->G->node_list[i].adj_last[0] = H->G->adj_last[i];
		H->G->node_list[i].adj2 = H->G->adj_last[i];
	}

	// Build vector of costs
	H->w.assign(w.begin(), w.end());

	// Build vector of indistinguishable colors
	H->K.assign(K.begin(), K.end());

	// Build the partition into indistinguishable colors
	H->C.assign(C.begin(), C.end());

	// Build Vk's
	H->V.resize(H->K.size());
	for (int i = 0; i < H->K.size(); ++i) {
		H->V[i].n_list = V[i].n_list;
		H->V[i].list = new nodepnt[H->V[i].n_list + 1];
		for (int j = 1; j <= H->V[i].n_list; ++j)
			H->V[i].list[j] = H->G->node_list + V[i].list[j]->name;
	}

	// Set branch information
	H->st = JOIN;
	H->branch_vertex_u = u;
	H->branch_vertex_v = v;

	// Build vector of contracted vertices
	H->vertex_mapping.assign(vertex_mapping.begin(), vertex_mapping.end());

#ifdef STABLE_POOL
	// Build the pools
	H->global_pool.resize(H->K.size());
	H->local_pool.resize(H->K.size());
	// Copy into H->global_pool the stables from global_pool
	for (int i = 0; i < H->K.size(); ++i) {
		for (auto &s: global_pool[i]) {
			nodepnt *stable = NULL;
			int n_stable = 0;
			H->translate_stable_set(i, s.list, s.n_list, &stable, &n_stable);
			H->global_pool[i].emplace_back(stable, n_stable);
		}
	}
#endif

	return H;

}

Graph *Graph::collapse_vertices(int u, int v) {

	// Build a new graph
	Graph *H = new Graph();

	// Vertices start at 1
	u = u+1;
	v = v+1;
	if (u > v) { // Swap to ensure u < v
		int t = u;
		u = v;
		v = t;
	}

	// Build Sewell's graph
	H->G = (MWSSgraphpnt) malloc(sizeof(MWSSgraph));
	allocate_graph(H->G, get_n_vertices() - 1); // We have one less vertex
	for(int i = 0; i <= H->G->n_nodes; i++)
		for(int j = 0; j <= H->G->n_nodes; j++)
			H->G->adj[i][j] = G->adj[i >= v ? i+1 : i][j >= v ? j+1 : j]; // Ignore the entry of v
	// Add the neighbours of v to the entry of u
	for (int j = 1; j < v; ++j) {
		if (G->adj[v][j]) {
			H->G->adj[u][j] = 1;
			H->G->adj[j][u] = 1;		
		}
	}
	for (int j = v+1; j <= get_n_vertices(); ++j) {
		if (G->adj[v][j]) {
			H->G->adj[u][j-1] = 1;
			H->G->adj[j-1][u] = 1;
		}
	}

	build_graph(H->G);

	// Set some information necessary to solve the MWSSP
	for(int i = 1; i <= H->G->n_nodes; i++) {
		MWIS_MALLOC(H->G->node_list[i].adj_last, H->G->n_nodes+1, nodepnt*);
		H->G->node_list[i].adj_last[0] = H->G->adj_last[i];
		H->G->node_list[i].adj2 = H->G->adj_last[i];
	}	

	// Build vector of costs
	H->w.assign(w.begin(), w.end());

	// Build vector of indistinguishable colors
	H->K.assign(K.begin(), K.end());

	// Build the partition into indistinguishable colors
	H->C.assign(C.begin(), C.end());

	// Build Vk's
	H->V.resize(H->K.size());
	for (int i = 0; i < H->K.size(); ++i) {
		H->V[i].list = new nodepnt[V[i].n_list + 1]; // At most as many vertices as V[i]
		H->V[i].n_list = 0;
		bool is_u = false;
		bool is_v = false;
		for (int j = 1; j <= V[i].n_list; ++j) {
			int w = V[i].list[j]->name; 
			if (w < u) {
				H->V[i].n_list++;
				H->V[i].list[H->V[i].n_list] = H->G->node_list + w;
			}
			else if (w == u) is_u = true;
			else if (w < v) {
				H->V[i].n_list++;
				H->V[i].list[H->V[i].n_list] = H->G->node_list + w;
			}
			else if (w == v) is_v = true;
			else {
				H->V[i].n_list++;
				H->V[i].list[H->V[i].n_list] = H->G->node_list + (w-1);
			}
		}
		if (is_u && is_v) {
			H->V[i].n_list++;
			H->V[i].list[H->V[i].n_list] = H->G->node_list + u;
		}
	}

	// Set branch information
	H->st = COLLAPSE;
	H->branch_vertex_u = u;
	H->branch_vertex_v = v;

	// Build vector of contracted vertices
	H->vertex_mapping.assign(vertex_mapping.begin(), vertex_mapping.end());
	for (int i = 1; i < H->vertex_mapping.size(); ++i) {
		if (H->vertex_mapping[i] < v)
			continue;
		else if (H->vertex_mapping[i] == v)
			H->vertex_mapping[i] = u;
		else
			H->vertex_mapping[i] = H->vertex_mapping[i] - 1;
	}

#ifdef STABLE_POOL
	// Build the pools
	H->global_pool.resize(H->K.size());
	H->local_pool.resize(H->K.size());
	// Copy into H->global_pool the stables from global_pool
	for (int i = 0; i < H->K.size(); ++i) {
		for (auto &s: global_pool[i]) {
			nodepnt *stable = NULL;
			int n_stable = 0;
			H->translate_stable_set(i, s.list, s.n_list, &stable, &n_stable);
			H->global_pool[i].emplace_back(stable, n_stable);
		}
	}
#endif

	return H;

}

int Graph::get_vertex_u() {
	return branch_vertex_u - 1;
}

int Graph::get_vertex_v() {
	return branch_vertex_v - 1;
}

int Graph::get_n_total_vertices() {
	return vertex_mapping.size() - 1;
}

int Graph::get_current_vertex(int v) {
	return vertex_mapping[v+1]-1;
}

BRANCH_STATUS Graph::get_branch_status() {
	return st;
}

Graph* Graph::choose_color(int v, int k) {

	// Build a new graph
	Graph *H = new Graph();

	// Vertices start at 1
	v = v+1;

	// Reuse Sewell's graph
	H->G = G;

	// Build vector of costs
	H->w.assign(w.begin(), w.end());

	if (get_n_C(k) == 1) {

		// Build vector of indistinguishable colors
		H->K.assign(K.begin(), K.end());

		// Build the partition into indistinguishable colors
		H->C.assign(C.begin(), C.end());

		// Build Vk's
		H->V.resize(H->K.size());
		for (int i = 0; i < H->K.size(); ++i) {
			if (i == k) {
				// Remove N(v)
				int n_neighbours = get_n_neighbours(v-1, k);
				H->V[i].n_list = V[i].n_list - n_neighbours;
				H->V[i].list = new nodepnt[H->V[i].n_list + 1];
				int index = 1;
				for (int j = 1; j <= V[i].n_list; ++j) {
					int u = V[i].list[j]->name;
					if (u == v || (u != v && !G->adj[u][v]))
						H->V[i].list[index++] = H->G->node_list + u;
				}
			}
			else {
				H->V[i].list = new nodepnt[V[i].n_list + 1];
				H->V[i].n_list = 0;
				for (int j = 1; j <= V[i].n_list; ++j) {
					int u = V[i].list[j]->name;
					if (u != v)
						H->V[i].list[++(H->V[i].n_list)] = H->G->node_list + u;
				}
			}
		}
	}
	else {

		// Build vector of indistinguishable colors
		H->K.resize(K.size() + 1);
		for (int i = 0; i < K.size(); ++i) {
			if (i == k)
				H->K[i] = C[k][1]; // Skip C[k][0] 
			else
				H->K[i] = K[i];
		}
		H->K[K.size()] = C[k][0];

		// Build the partition into indistinguishable colors
		H->C.resize(C.size() + 1);
		for (int i = 0; i < K.size(); ++i) {
			if (i == k)
				H->C[i].assign(C[i].begin() + 1, C[i].end());
			else
				H->C[i].assign(C[i].begin(), C[i].end());
		}
		H->C[C.size()].push_back(C[k][0]);

		// Build Vk's
		H->V.resize(V.size() + 1);
		for (int i = 0; i < V.size(); ++i) {
			H->V[i].list = new nodepnt[V[i].n_list + 1];
			H->V[i].n_list = 0;
			for (int j = 1; j <= V[i].n_list; ++j) {
				int u = V[i].list[j]->name;
				if (u != v)
					H->V[i].list[++(H->V[i].n_list)] = H->G->node_list + u;
			}
		}
		int n_neighbours = get_n_neighbours(v-1, k);
		H->V[V.size()].n_list = V[k].n_list - n_neighbours;
		H->V[V.size()].list = new nodepnt[H->V[V.size()].n_list + 1];
		int index = 1;
		for (int j = 1; j <= V[k].n_list; ++j) {
			int u = V[k].list[j]->name;
			if (u == v || (u != v && !G->adj[u][v]))
				H->V[V.size()].list[index++] = H->G->node_list + u;
		}

	}

	// Set branch information
	H->st = CHOOSE;
	H->branch_vertex_v = v;
	H->branch_color_k = k;

#ifdef STABLE_POOL
	// Build the pools
	H->global_pool.resize(H->K.size());
	H->local_pool.resize(H->K.size());
	// Copy into H->global_pool the stables from global_pool
	for (int i = 0; i < H->K.size(); ++i) {
		for (auto &s: global_pool[i]) {
			nodepnt *stable = NULL;
			int n_stable = 0;
			H->translate_stable_set(i, s.list, s.n_list, &stable, &n_stable);
			H->global_pool[i].emplace_back(stable, n_stable);
		}
	}
#endif

	return H;

}

Graph* Graph::remove_color(int v, int k) {

	// Build a new graph
	Graph *H = new Graph();

	// Vertices start at 1
	v = v+1;

	// Reuse Sewell's graph
	H->G = G;

	// Build vector of costs
	H->w.assign(w.begin(), w.end());

	// Build vector of indistinguishable colors
	H->K.assign(K.begin(), K.end());

	// Build the partition into indistinguishable colors
	H->C.assign(C.begin(), C.end());

	// Build Vk's
	H->V.resize(H->K.size());
	for (int i = 0; i < H->K.size(); ++i) {
		if (i == k) {
			// Skip v
			H->V[i].n_list = V[i].n_list - 1;
			H->V[i].list = new nodepnt[H->V[i].n_list + 1];
			int index = 1;
			for (int j = 1; j <= V[i].n_list; ++j) {
				int u = V[i].list[j]->name; 
				if (u != v)
					H->V[i].list[index++] = H->G->node_list + u;
			}

		}
		else {
			H->V[i].n_list = V[i].n_list;
			H->V[i].list = new nodepnt[H->V[i].n_list + 1];
			for (int j = 1; j <= H->V[i].n_list; ++j)
				H->V[i].list[j] = H->G->node_list + V[i].list[j]->name;
		}
	}

	// Set branch information
	H->st = REMOVE;
	H->branch_vertex_v = v;
	H->branch_color_k = k;

#ifdef STABLE_POOL
	// Build the pools
	H->global_pool.resize(H->K.size());
	H->local_pool.resize(H->K.size());
	// Copy into H->global_pool the stables from global_pool
	for (int i = 0; i < H->K.size(); ++i) {
		for (auto &s: global_pool[i]) {
			nodepnt *stable = NULL;
			int n_stable = 0;
			H->translate_stable_set(i, s.list, s.n_list, &stable, &n_stable);
			H->global_pool[i].emplace_back(stable, n_stable);
		}
	}
#endif

	return H;

}


Graph* Graph::choose_indistinguishable_color(int v, int k) {

	// Build a new graph
	Graph *H = new Graph();

	// Vertices start at 1
	v = v+1;

	// Reuse Sewell's graph
	H->G = G;

	// Build vector of costs
	H->w.assign(w.begin(), w.end());

	// Build vector of indistinguishable colors
	H->K.assign(K.begin(), K.end());

	// Build the partition into indistinguishable colors
	H->C.assign(C.begin(), C.end());

	// Build Vk's
	H->V.resize(H->K.size());
	for (int i = 0; i < H->K.size(); ++i) {
		if (i == k) {
			if (get_n_C(k) == 1) {
				// Remove N(v)
				int n_neighbours = get_n_neighbours(v-1, k);
				H->V[i].n_list = V[i].n_list - n_neighbours;
				H->V[i].list = new nodepnt[H->V[i].n_list + 1];
				int index = 1;
				for (int j = 1; j <= V[i].n_list; ++j) {
					int u = V[i].list[j]->name;
					if (u == v || (u != v && !G->adj[u][v]))
						H->V[i].list[index++] = H->G->node_list + u;
				}
			}
			else {
				// Copy Vk
				H->V[i].n_list = V[i].n_list;
				H->V[i].list = new nodepnt[H->V[i].n_list + 1];
				for (int j = 1; j <= V[i].n_list; ++j) {
					int u = V[i].list[j]->name;
					H->V[i].list[j] = H->G->node_list + u;
				}				
			}
		}
		else {
			H->V[i].list = new nodepnt[V[i].n_list + 1];
			H->V[i].n_list = 0;
			for (int j = 1; j <= V[i].n_list; ++j) {
				int u = V[i].list[j]->name;
				if (u != v)
					H->V[i].list[++(H->V[i].n_list)] = H->G->node_list + u;
			}
		}
	}

	// Set branch information
	H->st = CHOOSE;
	H->branch_vertex_v = v;
	H->branch_color_k = k;

#ifdef STABLE_POOL
	// Build the pools
	H->global_pool.resize(H->K.size());
	H->local_pool.resize(H->K.size());
	// Copy into H->global_pool the stables from global_pool
	for (int i = 0; i < H->K.size(); ++i) {
		for (auto &s: global_pool[i]) {
			nodepnt *stable = NULL;
			int n_stable = 0;
			H->translate_stable_set(i, s.list, s.n_list, &stable, &n_stable);
			H->global_pool[i].emplace_back(stable, n_stable);
		}
	}
#endif

	return H;

}
