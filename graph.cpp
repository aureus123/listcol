#include "graph.h"
#include "io.h"
#include <algorithm>
#include <stdlib.h>

#define RANDOMCOSTS

// Construct graph from .graph, .cost and .list
Graph(char *graph_filename, char *cost_filename, vector<int>& cost_list, char *list_filename)
{
    // Read adjacency list
    read_graph(graph_filename, adj);
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
