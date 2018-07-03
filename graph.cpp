#include "graph.h"
#include "io.h"
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <set>

//#define RANDOMCOSTS

// Construct graph from .graph, .cost and .list
Graph::Graph(char *graph_filename, char *cost_filename, vector<int>& cost_list, char *list_filename)
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
    // Construct subgraphs G_k
    V.resize(colors);
    for(int v = 0; v < vertices; v++)
        for(int k: L[v])
            V[k].push_back(v);    
}

// urnd - generate numbers with uniform random distribution
//  if flag=false, the interval is [a, b]
//  if flag=true, the interval is [a, b)
static float urnd(float a, float b, bool flag)
{
	return a + rand() * (b - a) / (float)(RAND_MAX + (flag ? 1 : 0));
}

// Construct graph from .graph (.cost and .list are automatically generated)
Graph::Graph(char *graph_filename, vector<int>& cost_list)
{
    // Read adjacency list
    edges = read_graph(graph_filename, adj);
    // Update number of vertices
    vertices = adj.size();
    
    // Update number of colors
    colors = vertices;    
    // Generate costs of colors
    cost_list.resize(colors);
    for (int k = 0; k < colors; k++) {
#ifdef RANDOMCOSTS
        cost_list[k] = (int)urnd(1, 10, false);
#else
        cost_list[k] = 1;
#endif
    }
    
    // Generate list of colors
    L.resize(vertices, vector<int> (colors));
	for (int v = 0; v < vertices; v++)
		for (int k = 0; k < colors; k++)
            L[v][k] = k;
    // Construct subgraphs G_k
    V.resize(colors);
    for(int v = 0; v < vertices; v++)
        for(int k: L[v])
            V[k].push_back(v);     
}

bool Graph::is_edge (int u, int v) {
    return (find(adj[u].begin(), adj[u].end(), v) != adj[u].end());
}

void Graph::get_Vk (int k, vector<int>& Vk) {
    Vk = V[k];
    return;
}

bool Graph::have_common_color(int u, int v) {

    for (int k: L[u])
        if (find(L[v].begin(), L[v].end(), k) != L[v].end())
            return true;

    return false;

}

void Graph::show_instance(vector<int>& costs_list) 
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
		cout << " " << k << "->" << costs_list[k];
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

// Generate a new graph from the current one by collapsing vertices u and v
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

    return;
}
