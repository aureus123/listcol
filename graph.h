#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <vector>
#include <map>

#define COLORS_DELETION

using namespace std;

class Graph {

    public:
    int vertices;
    int colors;
    int edges;

    Graph(char *graph_filename, char *cost_filename, vector<int>& cost_list, char *list_filename);  // Constructor

    bool is_edge(int u, int v);
    int get_Lv_size(int v);
    void get_Lv(int v, vector<int>& Lv);
    int get_Vk_size(int k);
    void get_Vk(int k, vector<int>& Vk);
    int get_cost(int k);
    bool have_common_color(int u, int v);
    void get_new_vertex(vector<int>&);
    bool check_coloring(vector<int>& f);

    void show_instance(vector<int>& costs_list);
    void show_statics();
    
    void join_vertices(int u, int v);
    void collapse_vertices(int u, int v);
    void color_vertex(int v, int k);
    void set_Lv(int v, vector<int>& Lv);

    bool coloring_heuristic(vector<vector<int>>& stables_set);

    private:
    vector<vector <int> > adj;    // adj[v]: open neighborhood of vertex v
    vector<vector <int> > L;      // L[v]: allowed colors for vertex v
    vector<vector <int> > V;      // V[k]: vertices coloreable with color k
    vector<int>& cost_list;       // cost_list[k]: cost associated to color k
    vector<int> new_vertex;       // new_vertex[v]: mapping from original vertices to current vertices

#ifdef COLORS_DELETION
    public:
    int get_right_hand_side(int k);
    int get_eq_colors(int k, int d);
    void delete_equal_colors();   // If colors k1 < k2 have got the same associated subgraph and same cost, then
                                  // erase k2 from every list of L and update right_hand_side of k1
    private:
    vector<int> right_hand_side;      // right_hand_side[k]: right hand side of constraint associated to color k
    map<int,vector<int> > eq_colors;  // eq_colors[k]: deleted colors and replaced by k
#endif

};

#endif
