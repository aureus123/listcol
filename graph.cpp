#include "graph.h"
#include "io.h"
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <set>
#include <climits>

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

    // If |L[v]| = 1, then that color must be erase from the list of neighbors of v
    for (int v = 0; v < vertices; ++v)
        if (L[v].size() == 1)
            for (int u: adj[v]) {
                auto it = find(L[u].begin(), L[u].end(), L[v][0]);
                if (it != L[u].end())
                    L[u].erase(it);
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
    right_hand_side.resize(colors,-1);
#endif

}

bool Graph::is_edge (int u, int v) {
    return (find(adj[u].begin(), adj[u].end(), v) != adj[u].end());
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
	cout << "  |V| = " << vertices << ", |E| = " << edges << " (density = " << density << "%), |C| = " << colors << "." << endl;
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
	prom = 0.0;
	for (int k = 0; k < colors; k++) prom += (float)V[k].size();
	prom /= (float)colors;
	sigma = 0.0;
	for (int k = 0; k < colors; k++) {
		float substr = (float)V[k].size() - prom;
		sigma += substr * substr;
	}
	sigma /= (float)(colors - 1);
	cout << "  Behaviour of |V(Gk)| ---> prom = " << prom << ", sigma = " << sqrt(sigma) << "." << endl;
	set_color(7);

    return; 
}

// Add an edge between vertices u and v
// CUATION: This function modify the current graph
void Graph::join_vertices (int u, int v) {

    if (is_edge(u,v)) bye("Branching error");

    edges++;
    adj[u].push_back(v);
    adj[v].push_back(u);

    return;
}

// Collapse vertices u and v
// CUATION: This function modify the current graph
void Graph::collapse_vertices (int u, int v) {

    if (is_edge(u,v)) bye("Branching error");

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
    if (intersection.size() == 0)
        bye("Branching error");

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

    // If |L[v]| = 1, then that color must be erase from the list of neighbors of v
    if (L[u].size() == 1)
        for (int w: adj[u]) {
            auto it = find(L[w].begin(), L[w].end(), L[u][0]);
            if (it != L[w].end()) {
                L[w].erase(it);
                V[L[u][0]].erase(find(V[L[u][0]].begin(), V[L[u][0]].end(), w));
            }
        }

    return;
}

// Color a vertex
// CUATION: This function modify the current graph
void Graph::color_vertex(int v, int k) {

    // Color vertex v with color k
    // (this is equal to set L[v] = k and delete k from the neighbors' lists)

    // Erase v from V[j] forall j in L[v] - k
    for (int j: L[v])
        if (j != k)
            V[j].erase(find(V[j].begin(), V[j].end(), v));

    L[v].clear();
    L[v].push_back(k);

    for (int u: adj[v]) {
        auto it = find(L[u].begin(), L[u].end(), k);
        if (it != L[u].end()) {
            L[u].erase(it);
            V[k].erase(find(V[k].begin(), V[k].end(), u));
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

bool Graph::coloring_heuristic(vector<vector<int>>& stables_set) {

    // Initialize adjacency matrix for quick indexing
    vector<vector<bool>> M (vertices, vector<bool>(vertices,false));
    for (int v = 0; v < vertices; ++v)
        for (int u: adj[v])
            M[v][u] = true;

    // Degrees array
    vector<int> deg (vertices, 0);
    for (int v = 0; v < vertices; ++v)
        deg[v] = adj[v].size();
    
    vector<vector<int>> VV = V; // Copy V[k]'s

    vector<bool> colored (vertices, false); // Currently colored vertices
    int colored_size = 0; // Currently colored vertices size
    vector<int> used (colors,0);

    while (colored_size != vertices) {

        // Sort each VV[k] in non-decreasing degree order
        for (int k = 0; k < colors; ++k)
            sort(VV[k].begin(), VV[k].end(), [deg](int u, int v) {return deg[u] < deg[v];} );

        // Find a big maximal stable with uncolored vertices of some VV[k].
        // Note: a greedy approach is used
        vector<int> max_stable;
        int color;
        for (int k = 0; k < colors; ++k) {

#ifdef COLORS_DELETION
            if (used[k] >= abs(right_hand_side[k])) continue;
#elif
            if (used[k] >= 1) continue;
#endif

            vector<int> stable;
            vector<int> candidates;
            for (int v: VV[k])
                if (!colored[v])
                    candidates.push_back(v);
            while (!candidates.empty()) {
                int v = candidates.back();
                stable.push_back(v);
                candidates.pop_back();
                // Remove N(v) from candidates
                for (auto it = candidates.begin(); it != candidates.end();)
                    if (M[v][*it]) it = candidates.erase(it);
                    else ++it;
            }
            if (stable.size() > max_stable.size()) {
                max_stable = stable;
                color = k;
            }
        }

        if (max_stable.empty())
            return false;

        // Color vertices from max_stable
        for (int v: max_stable) {
            colored[v] = true;
            colored_size++;
        }

        // Rearrange degrees
        for (int v: max_stable)
            for (int u = 0; u < vertices; u++)
                if (M[v][u])
                    deg[u]--;

        max_stable.push_back(color);
        used[color]++;
        stables_set.push_back(max_stable);

    }

    return true;

}

bool Graph::coloring_heuristic2(vector<vector<int>>& stables_set) {

    // For each vertex v, find a maximal stable that covers it
    // The stable is chosen from the larger Gk with k in L(v)
    // The stable is constructed in a greedy fashion:
    //   the vertices are sorted in non-decreasing degree order

    for (int v = 0; v < vertices; ++v) {

        if (L[v].size() == 0) bye("Heuristic error: empty list");

        // Choose unused color
        int r = rand() % L[v].size();
        int k = L[v][r];

        // Construct stable
        vector<int> stable;
        stable.push_back(v);
        vector<bool> used (V[k].size(), false);

        do {

            // Mark N[v] as used
            int u = stable.back();
            for (int i = 0; i < V[k].size(); ++i)
                if (V[k][i] == u)
                    used[i] = true;
            for (int w: adj[u])
                for (int i = 0; i < V[k].size(); ++i)
                    if (V[k][i] == w)
                        used[i] = true;

            // Choose next vertex
            vector<int> candidates;
            for (int i = 0; i < V[k].size(); ++i)
                if (!used[i])
                    candidates.push_back(V[k][i]);

            if (candidates.size() == 0)
                break;

            r = rand() % candidates.size();
            stable.push_back(candidates[r]);

        }
        while (true);

        cout << "Stable: ";
        for (int i: stable)
            cout << " " << i;
        cout << endl;

        stable.push_back(k); // push color at back
        stables_set.push_back(stable);
        
    }

    return true;

}

#ifdef COLORS_DELETION
void Graph::delete_equal_colors() {

    // CUATION: V[k]'s must be sorted with the same order function

    for (int k1 = 0; k1 < colors-1; ++k1)
        for (int k2 = k1+1; k2 < colors; ++k2) {
            if (right_hand_side[k2] == 0) continue;
            if ((cost_list[k1] == cost_list[k2]) && (V[k1] == V[k2])) {

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
