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
    bool is_edge(int u, int v);
    int get_Lv_size(int v);
    void get_Lv(int v, vector<int>& Lv);
    void get_Vk(int k, vector<int>& Vk);
    int get_cost(int k);
    bool have_common_color(int u, int v);
    void get_new_vertex(vector<int>&);
    bool check_coloring(vector<int>& f);

    void show_instance(vector<int>& costs_list);
    void show_statics();
    
    void join_vertices(int u, int v);
    void collapse_vertices(int u, int v);
    int celim(int v);
    void color_vertex(int v, int k);

    private:
    vector<vector <int> > adj;    // adj[v]: open neighborhood of vertex v
    vector<vector <int> > L;      // L[v]: allowed colors for vertex v
    vector<vector <int> > V;      // V[k]: vertices coloreable with color k
    vector<int>& cost_list;
    vector<int> new_vertex;       // new_vertex[v]: mapping from original vertices to current vertices
};

#endif
