#include "graph.h"
#include <algorithm>


#include <iostream>


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

int Graph::get_cost(int k) {
    return cost[k];
}
