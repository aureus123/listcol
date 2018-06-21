#include "graph.h"
#include "io.h"
#include <algorithm>


#include <iostream>

// Construct graph from .graph
Graph(char *filename)
{

    // Read adjacency list
    read_graph(graph_filename, adj);
    
    // Update number of vertices
    vertices = adj.size();

}

// Set lists of colors from .list
void Graph::set_L (char *filename) {
    read_list(filename,vertices,colors,L)
    V.resize(colors);
    for(int v = 0; v < vertices; v++)
        for(int k: L[v])
            V[k].push_back(v);
}

// Set lists of colors for classic coloring (the list of each vertex has every colors)
void Graph::set_L () {
    L.resize(vertices, vector<int> (colors));
	for (int v = 0; v < vertices; v++)
		for (int k = 0; k < colors; k++)
            L[v][k] = k;

    V.resize(colors, vector<int> (vertices));
    for(int v = 0; v < vertices; v++)
        for(int k: L[v])
            V[k][v] = v;
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
