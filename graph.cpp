#include "graph.h"

Graph:Graph (vector<pair<int,int> >& edges, vector<vector <int>>& color_list, vector<int>& cost_list)
    vertices(color_list.size()), colors(cost.size()), adj(vertices), L(color_list), V(colors), cost(cost_list)
{
    for (auto p: edges) {
        adj[p.first].push_back(p.second);
        adj[p.second].push_back(p.first);
    }

    for (int v = 0, v < vertices, v++)
        for (int k: L[v])
            V[k].push_back(v);
}

bool Graph:is_edge (int u, int v) {
    return (find(G.adj[u].begin(), G.adj[u].end(), v) != G.adj[u].end());
}

void Graph:get_Vk (int k, vector<int>& Vk) {
    Vk = V[k];
    return;
}

int Graph:get_cost(int k) {
    return cost[k];
}
