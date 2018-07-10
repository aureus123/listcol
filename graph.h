#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <vector>

using namespace std;

class Graph {

    public:
    int vertices;
    int colors;
    int edges;

    Graph(char *graph_filename, char *cost_filename, vector<int>& cost_list, char *list_filename);
    Graph(char *graph_filename, vector<int>& cost_list);
    bool is_edge(int u, int v);
    void get_Vk(int k, vector<int>& Vk);
    bool have_common_color(int u, int v);
    void get_new_vertex(vector<int>&);
    void check_coloring(vector<int>& f);

    void show_instance(vector<int>& costs_list);
    void show_statics();
    
    void join_vertices(int u, int v);
    void collapse_vertices(int u, int v);

    private:
    vector<vector <int> > adj;    // adj[v]: open neighborhood of vertex v
    vector<vector <int> > L;      // L[v]: allowed colors for vertex v
    vector<vector <int> > V;      // V[k]: vertices coloreable with color k
    vector<int> new_vertex;       // new_vertex[v]: mapping from original vertices to current vertices
};

#endif
