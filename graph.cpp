#include "graph.h"
#include "io.h"
#include <algorithm>
#include <stdlib.h>

#define RANDOMCOSTS

// Construct graph from .graph, .cost and .list
Graph(char *graph_filename, char *cost_filename, vector<int>& cost_list, char *list_filename)
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
float urnd(float a, float b, bool flag)
{
	return a + rand() * (b - a) / (float)(RAND_MAX + (flag ? 1 : 0));
}

// Construct graph from .graph (.cost and .list are automatically generated)
Graph(char *graph_filename, vector<int>& cost_list)
{
    // Read adjacency list
    read_graph(graph_filename, adj);
    // Update number of vertices
    vertices = adj.size();
    
    // Update number of colors
    colors = vertices;    
    // Generate costs of colors
    costs_list.resize(colors);
    for (int k = 0; k < colors; k++) {
#ifdef RANDOMCOSTS
        costs_list[k] = (int)urnd(1, 10, false);
#else
        costs_list[k] = 1;
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

int Graph::get_Vk_size(int k) {
    return V[k].size();
}


void Graph::get_Nv(int v, vector<int>& Nv) {
    Nv = adj[v];
    return;
}

void Graph::show_instance(vector<int>& costs_list) 
{
	set_color(2);

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

/*
Graph::Graph (vector<pair<int,int> >& edges, vector<vector <int>>& color_list, vector<int>& cost_list) :
    L(color_list), cost(cost_list)
{
    vertices = L.size();
    colors = cost.size();

    adj.resize(vertices); 
    V.resize(colors);

    for (auto p: edges) {
        adj[p.first].push_back(p.second);
        adj[p.second].push_back(p.first);
    }

    for (int v = 0; v < vertices; v++)
        for (int k: L[v])
            V[k].push_back(v);
}
*/
