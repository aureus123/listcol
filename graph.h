#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <vector>

using namespace std;

class Graph {

    public:
    int vertices;

    Graph(char *filename);
    void set_L(char *filename);
    void set_L();
    //Graph(vector<pair<int,int> >& edges, vector<vector<int> >& L, vector<int>& cost);
    bool is_edge(int u, int v);
    void get_Vk(int k, vector<int>& Vk);
    int get_Vk_size(int k);
    void get_Nv(int v, vector<int>& Nv);

    private:
    vector<vector <int> > adj;    // adj[v]: open neighborhood of vertex v
    vector<vector <int> > L;      // L[v]: allowed colors for vertex v
    vector<vector <int> > V;      // V[k]: vertices coloreable with color k
    
};

#endif
