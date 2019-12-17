#include "graph.h"
#include "io.h"
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <set>
#include <numeric> // iota function
#include <climits>
#include <math.h>

// Construct graph from .graph, .cost and .list
Graph::Graph(char *graph_filename, char *cost_filename, vector<int>& cost_list, char *list_filename) : cost_list(cost_list)
{
    // Read adjacency list
    edges = read_graph(graph_filename, adj);
    // Update number of vertices
    vertices = adj.size();

    // Read cost of colors
    read_cost(cost_filename, cost_list);
    // Update number of colors
    colors = cost_list.size();
    
    // Read list of colors
    read_list(list_filename, vertices, colors, L);

    // If |L[v]| = 1, then that color must be erase from the list of neighbors of v.
    // This is done recursively
    bool repeat = true;
    while (repeat) {
        repeat = false;
        for (int v = 0; v < vertices; ++v) {
            if (L[v].size() == 1) {
                for (int u: adj[v]) {
                    auto it = find(L[u].begin(), L[u].end(), L[v][0]);
                    if (it != L[u].end()) {
                        L[u].erase(it);
                        repeat = true;
                    }
                }
            }
        }
    }
    
    // Construct subgraphs G_k
    V.resize(colors);
    for(int v = 0; v < vertices; v++)
        for(int k: L[v])
            V[k].push_back(v);

    // At first, new_vertex is the identity mapping
    new_vertex.resize(vertices);
    iota(new_vertex.begin(), new_vertex.end(), 0);

#ifdef COLORS_DELETION
    // At first, right_hand_side = -1
    right_hand_side.resize(colors);
    for (int &k: right_hand_side)
        k = -1;
#endif

}

bool Graph::is_edge (int u, int v) {
    return (find(adj[u].begin(), adj[u].end(), v) != adj[u].end());
}

bool Graph::is_admissible(int u, int k) {
    return (find(L[u].begin(), L[u].end(), k) != L[u].end());
}

int Graph:: get_Lv_size (int v) {
    return L[v].size();
}

void Graph::get_Lv (int v, vector<int>& Lv) {
    Lv.resize(L[v].size());
    Lv = L[v];
    return;
}

int Graph::get_Vk_size(int k) {
    return V[k].size();
}

void Graph::get_Vk (int k, vector<int>& Vk) {
    Vk.resize(V[k].size());
    Vk = V[k];
    return;
}

int Graph::get_cost(int k) {
    return cost_list[k];
}

void Graph::get_adjv(int v, vector<int>& adjv) {
    adjv.resize(adj[v].size());
    adjv = adj[v];
    return;    
}

bool Graph::have_common_color(int u, int v) {

    for (int k: L[u])
        if (find(L[v].begin(), L[v].end(), k) != L[v].end())
            return true;

    return false;

}

void Graph::get_new_vertex(vector<int>& ret) {
    ret.resize(new_vertex.size());
    ret = new_vertex;
    return;
}

bool Graph::check_coloring(vector<int>& f) {

	// colors of each vertex belong to the list?
	for (int v = 0; v < vertices; v++) {
		int k_chosen = f[v];
        if (find(L[v].begin(), L[v].end(), k_chosen) == L[v].end()) {
			return false;
        }
	}

	// are there conflicting edges?
	for (int u = 0; u < vertices; u++)
        for (int v: adj[u])
    		if (f[u] == f[v]) {
			    return false;
		    }

	return true;
}

// Maximize the stable set in Gk

void Graph::maximize_stable_set(vector<int> &stable_set, int k) {

    for (unsigned int i = 0; i < V[k].size(); i++) {
        bool add = true;
        for (int s: stable_set) {
            if (s == V[k][i] || is_edge(s,V[k][i])) {
                add = false;
                break;
            }
        }
        if (add) {
            stable_set.push_back(V[k][i]);
        }
    }

}

// Get the weight of the weightest color in L(v) 

int Graph::get_Mv(int v) {

    int max = -1;

    for (int k: L[v]) {
        if (get_cost(k) > max) {
            max = get_cost(k);
        }
    }

    if (max == -1) {
        std::cout << "Error in get_Mv()" << std::endl;
        abort();
    }

    return max;

}

void Graph::show_instance() 
{
	set_color(2);

    cout << "Vertices = " << vertices << endl;
    cout << "Edges = " << edges << endl;

	cout << "Neighborhoods:" << endl;
	int maxdelta = 0;
	for (int v = 0; v < vertices; v++) {
		cout << "N(" << v << ") = {";
		int degree = adj[v].size();
		if (degree > maxdelta) maxdelta = degree;
		for (int d = 0; d < degree; d++) cout << " " << adj[v][d];
		cout << " }, degree = " << degree << endl;
	}
	cout << "Maximum degree = " << maxdelta << endl;

	cout << "Vector of costs: {";
	for (int k = 0; k < colors; k++) {
		cout << " " << k << "->" << cost_list[k];
	}
	cout << " }, colors = " << colors << endl;

	cout << "List of colors:" << endl;
	for (int v = 0; v < vertices; v++) {
		cout << "L(" << v << ") = {";
		for (unsigned int s = 0; s < L[v].size(); s++) cout << " " << L[v][s];
		cout << " }" << endl;
	}

	cout << "Gk subgraphs:" << endl;
	for (int k = 0; k < colors; k++) {
		cout << "V(" << k << ") = {";
		for (unsigned int s = 0; s < V[k].size(); s++) cout << " " << V[k][s];
		cout << " }" << endl;
	}

    return;
}

void Graph::show_statics() {

	/* Show some basic statistics */
	set_color(6);
	cout << "Statistics:" << endl;
	int clique_size = vertices * (vertices - 1) / 2;
	float density = 100.0 * (float)edges / (float)clique_size;
	cout << "  |V| = " << vertices << ", |E| = " << edges << ", |C| = " << colors << "." << endl;

    /* New densities */
    int count1 = 0;
    int count2 = 0;
    for (int v = 0; v < vertices - 1; v++) 
    {
        for (int u = v+1; u < vertices; u++)
        {
            if (is_edge(u,v))
            {
                if (have_common_color(u,v)) 
                {
                    count1++;
                }
                count2++;
            }
            else
            {
                if (!have_common_color(u,v))
                {
                    count2++;
                }
            }
        }
    }
    float dG1 = 100.0 * (float)count1 / (float)clique_size;
    float dG2 = 100.0 * (float)count2 / (float)clique_size;
    cout << "  density(G) = " << density << "%, density(G') = " << dG1 << "%, density(G'') = " << dG2 << "%" << endl;

	/* Average and standard deviation of size of lists */
	float prom = 0.0;
	for (int v = 0; v < vertices; v++) prom += (float)L[v].size();
	prom /= (float)vertices;
	float sigma = 0.0;
	for (int v = 0; v < vertices; v++) {
		float substr = (float)L[v].size() - prom;
		sigma += substr * substr;
	}
	sigma /= (float)(vertices - 1);
	cout << "  Behaviour of |L(v)| ---> prom = " << prom << ", sigma = " << sqrt(sigma) << "." << endl;
	/* Average and standard deviation of vertices of Gk */
	float prom1 = 0.0;
	for (int k = 0; k < colors; k++) prom1 += (float)V[k].size();
	prom1 /= (float)colors;
	float sigma1 = 0.0;
	for (int k = 0; k < colors; k++) {
		float substr = (float)V[k].size() - prom1;
		sigma1 += substr * substr;
	}
	sigma1 /= (float)(colors - 1);
	cout << "  Behaviour of |V(Gk)| ---> prom = " << prom1 << ", sigma = " << sqrt(sigma1) << "." << endl;
	set_color(7);

    return; 
}

// Add an edge between vertices u and v
// CUATION: This function modify the current graph
void Graph::join_vertices (int u, int v) {

    edges++;
    adj[u].push_back(v);
    adj[v].push_back(u);

    return;
}

// Collapse vertices u and v
// CUATION: This function modify the current graph
void Graph::collapse_vertices (int u, int v) {

    if (u > v) {
        int temp = u;
        u = v;
        v = temp;
    }

    // Merge adj[u] and adj[v]
    for (int j: adj[v]) {
        vector<int>::iterator it = find(adj[u].begin(), adj[u].end(), j);
        if (it == adj[u].end()) {
            // j is only adjacent to v
            adj[u].push_back(j);
            *(find(adj[j].begin(), adj[j].end(), v)) = u; // change v to u
        }
        else {
            // j is adjacent to u and v
            edges--;
            adj[j].erase(find(adj[j].begin(), adj[j].end(), v));
        }
    }
    adj.erase(adj.begin()+v);
    vertices--;

    // Rearrange indices of adjacency list
    for(int i = 0; i < vertices; i++)
        for (int& j: adj[i])
            if (j > v)
                j--;

    // Intersection of L(u) and L(v)
    set<int> intersection;
    for (int k: L[u])
        if (find(L[v].begin(), L[v].end(), k) != L[v].end())
            intersection.insert(k);
    if (intersection.size() == 0) {
        bye("Branching error");
    }

    // Delete u from V[k] with k in L(u)/intersection
    for (int k: L[u])
        if (find(intersection.begin(), intersection.end(), k) == intersection.end())
            V[k].erase(find(V[k].begin(), V[k].end(), u));
    // Delete v from V[k] with k in L(v)/intersection
    for (int k: L[v])
        if (find(intersection.begin(), intersection.end(), k) == intersection.end())
            V[k].erase(find(V[k].begin(), V[k].end(), v));
    // Delete v from V[k] with k in intersection
    for (int k: intersection)
        V[k].erase(find(V[k].begin(), V[k].end(), v));

    // Rearrange subgraphs Gk
    for (int k = 0; k < colors; k++)
        for (int& w: V[k]) {
            if (w == v)
                w = u;
            else if (w > v)
                w--;
        }

    // Rearrange L
    L[u].clear();
    L[u].reserve(intersection.size());
    for (int k: intersection)
        L[u].push_back(k);
    L.erase(L.begin()+v);

    // Rearrange new_vertex
    for (int& i : new_vertex)
        if (i == v)
            i = u;
        else if (i > v)
            i--;

    // If |L[v]| = 1, then that color may be erase from the list of neighbors of v
    // Cuando hay colores repetidos este chequeo NO DEBE HACERSE.
    // Suponer que L[v] = {k} y k tiene "multiplicidad" 2.
    // Luego puede pintarse a v con un estable S1 de k y
    // alguno de sus vecinos pueden pintarse con otro estable S2 de k.
    // Luego no se podía borrar a k de las listas de los vecinos de v.
	bool repeat = false;
    if (L[u].size() == 1) {
        if (get_right_hand_side(L[u][0]) == -1) {
            // Erase color L[u][0] from the lists of the neighbors of u
            for (int w: adj[u]) {
                auto it = find(L[w].begin(), L[w].end(), L[u][0]);
                if (it != L[w].end()) {
                    L[w].erase(it);
                    V[L[u][0]].erase(find(V[L[u][0]].begin(), V[L[u][0]].end(), w));
                    repeat = true;
                }
            }
        }
    }
    // The previous is done recursively
    while (repeat) {
        repeat = false;
        for (int w = 0; w < vertices; ++w) {
            if (L[w].size() == 1) {
                if (get_right_hand_side(L[w][0]) == -1) {
                    // Erase color L[w][0] from the lists of the neighbors of w
                    for (int z: adj[w]) {
                        auto it = find(L[z].begin(), L[z].end(), L[w][0]);
                        if (it != L[z].end()) {
                            L[z].erase(it);
                            V[L[w][0]].erase(find(V[L[w][0]].begin(), V[L[w][0]].end(), z));
                            repeat = true;
                        }
                    }
                }
            }
        }
    }

    return;
}

// Color a vertex
// CUATION: This function modify the current graph
void Graph::color_vertex(int v, int k) {

    // Color vertex v with color k
    // (this is equal to set L[v] = k and delete k from the neighbors' lists)

    // Erase v from V[j] forall j in L[v] - {k}
    for (int j: L[v])
        if (j != k)
            V[j].erase(find(V[j].begin(), V[j].end(), v));

    L[v].clear();
    L[v].push_back(k);

    if (get_right_hand_side(L[v][0]) == -1) {
        for (int u: adj[v]) {
            auto it = find(L[u].begin(), L[u].end(), k);
            if (it != L[u].end()) {
                L[u].erase(it);
                V[k].erase(find(V[k].begin(), V[k].end(), u));
            }
        }
    }

    return;

}

void Graph::set_Lv(int v, vector<int>& Lv) {

    for (int k : L[v])
        if (find(Lv.begin(), Lv.end(), k) == Lv.end()) {
            V[k].erase(find(V[k].begin(), V[k].end(), v));
        }

    L[v] = Lv;

    return;
}

void Graph::coloring_heuristic(list<vector<int>>& stable_sets) {

    int n_colored = 0;
    vector<bool> colored (vertices, false);

    for (int k = 0; k < colors; ++k) {

        int multiplicity = 1;
        #ifdef COLORS_DELETION
        multiplicity = -right_hand_side[k];
        #endif

        for (int i = 0; i < multiplicity; ++i) {

            vector<int> stable;
            stable.reserve(V[k].size());

            for (int v: V[k]) {

                if (colored[v]) {continue;}

                // Try to add v to the stable

                bool ok = true;
                for (int s: stable) {
                    if (is_edge(v,s)) {
                        ok = false;
                        break;
                    }
                }
                if (ok) {
                    stable.push_back(v);
                    n_colored++;
                    colored[v] = true;
                    if (n_colored == vertices) {break;}
                }

            }

            if (!stable.empty()) {
                maximize_stable_set(stable, k);
                stable.push_back(k); // Push color
                stable_sets.push_back(stable);
            }

        }

    }

    // Fictional colors

    for (int v = 0; v < vertices; ++v) {

        if (!colored[v]) {
            vector<int> stable;
            stable.push_back(v);
            stable.push_back(colors+1); // Push fictional color
            stable_sets.push_front(stable); // Push front
        } 

    }

    return;

}

#ifdef COLORS_DELETION
void Graph::delete_equal_colors() {

    // CUATION: V[k]'s must be sorted with the same order function

    for (int k1 = 0; k1 < colors-1; ++k1)
        for (int k2 = k1+1; k2 < colors; ++k2) {
            if (get_right_hand_side(k2) == 0) continue;
            if ((get_cost(k1) == get_cost(k2)) && (V[k1] == V[k2])) {

                // Delete k2 from every list
                for (int v: V[k2]) {
                    L[v].erase(find(L[v].begin(), L[v].end(), k2));
                }
                V[k2].clear(); // Clear V[k2]

                // Update right-hand side of k1 and k2
                right_hand_side[k1]--;
                right_hand_side[k2] = 0;

                // Updated eq_colors of k1
                eq_colors[k1].push_back(k2);

            }
        }

}

int Graph::get_right_hand_side(int k) {
    return right_hand_side[k];
}

int Graph::get_eq_colors(int k, int d) {
    return eq_colors[k][d];
}

#endif
