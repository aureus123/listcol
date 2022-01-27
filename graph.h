#ifndef _GRAPH_H_
#define _GRAPH_H_

#ifndef VISUALC
#include "mwis_sewell/mwss.h"
#else
extern "C" {
#include "mwis_sewell/mwss.h"
}
#endif

/* Additional structure for arrays of nodepnt */
typedef struct nodepntArray {
	nodepnt* list;
	int n_list;
	nodepntArray() {};
	nodepntArray(nodepnt *list, int n_list) : list(list), n_list(n_list) {};
	~nodepntArray() {};
} nodepntArray;

#include <vector>
#include <list>
#include <set>

#define INTFACTOR 1000.0
#define THRESHOLD 0.1
#define MAXTIME_MWIS 300.0

//#define STABLE_POOL

#ifndef BRANCHING_STRATEGY
#define BRANCHING_STRATEGY 0  // 0: Branching on edges
                              // 1: Branching on colors
#endif

#ifndef PREPROCESSING
#define PREPROCESSING 0 // 0: never
                        // 1: preprocess and do not change G
                        // 2: preprocess and change G

#endif

enum BRANCH_STATUS {NONE, JOIN, COLLAPSE, CHOOSE, REMOVE};

class Graph {

  public:

  // Build graph
  Graph();
  Graph(char *filename_graph, char *filename_costs, char *filename_lists);

  // Destructor
  ~Graph();

  // Get number of vertices
  int get_n_vertices();

  // Get number of colors
  int get_n_colors();

  // Get size of V[i]
  int get_n_V(int i);

  // Get the j-th color of C[i] 
  int get_C(int i, int j);

  // Get size of C[i]
  int get_n_C(int i);

  // Get cost of color K[i]
  int get_color_cost(int i);

  // Get the max w_k such that vertex v belongs to V_k
  int get_m(int v);

  // Get N(v) \cap Vk
  int get_n_neighbours(int v, int k);

  // Get W1 such that W1[v] = |Lv \cap K| for all v
  void get_W1(std::vector<int> &W1);
  void get_S(std::vector<std::set<int>> &S);

  // Are vertices u and v adjacent?
  bool is_edge(int u, int v);

  // Set the weight of a vertex
  void set_vertex_weight(int v, double y);

  // Print graph
  void print_graph();
	
  // Preprocess the instance
  void preprocess_instance();
  void preprocess_instance(int v, int j);

  // Solve the MWSSP in V[i]: ¿Does it exists a maximum stable set in V[i] with a weight greater than goal?
  //  If the answer is yes, then the result is saved in best_stable and n_best_stable
  bool solve_MWSSP(int i, double goal, nodepnt **best_stable, int *n_best_stable);

  // Check stable set
  bool check_stable(int i, nodepnt *stable, int n_stable);
  
  // Check coloring TODO
  bool check_coloring(std::vector<int> &);

  // LP's stable representation (bool *) to Graph's stable representation (nodepnt *), and viceversa
  void column_to_stable (bool *column, int size, nodepnt **stable);
  void stable_to_column (nodepnt *stable, int size, bool **column);

  // Translate a stable set of the father to a stable set of the current son
  void translate_stable_set (int color, nodepnt *stable_father, int n_stable_father, nodepnt **stable_son, int *n_stable_son);

  // Find a stable set covering
  void coloring_heuristic(std::vector<std::list<nodepntArray>>& stable_sets);

  // Copy the local pool to the global pool
  void update_pool();

  // Build a new graph by joining the vertices u and v
  Graph *join_vertices(int u, int v);

  // Build a new graph by collapsing the vertices u and v
  Graph *collapse_vertices(int u, int v);

  // Build a new graph where vertex v has L'[v] = C[k]
  Graph* choose_color(int v, int k);

  // Build a new graph where vertex v has L'[v] = L[v] \ {C[k] : k \in colors}
  Graph* remove_color(int v, std::set<int> &colors);

  // Get branch status
  BRANCH_STATUS get_branch_status();

  // Get branching's vertices
  int get_vertex_u();
  int get_vertex_v();

  int get_n_total_vertices();     // Get the number of vertices of the original graph
  int get_current_vertex(int i);  // Get the current vertex that stands for the original vertex i
  int get_precoloring(int i);	  // Get the color of the original vertex i
  double get_precoloring_value(); // Get precoloring value

  private:

  MWSSgraphpnt G;                     // Pointer to Sewell's graph
                                      // Recall that vertices start at 1
  std::vector<int> w;                 // Cost vector 
  std::vector<int> K;                 // Vector of indistinguishable colors (K is a subset of [0,...,w.size() - 1])
  std::vector<std::vector<int>> C;    // For each 0 <= i < K.size(), C[i] is the vector with the equivalence class of color K[i]
  std::vector<nodepntArray> V;        // For each 0 <= i < K.size(), V[i] has all the vertices having K[i] in their list

  std::vector<std::list<nodepntArray>> global_pool; // For each 0 <= i < K.size(), global_pool[i] is a list with all the maximal stable sets of V[i] globally found
  std::vector<std::list<nodepntArray>> local_pool;  // For each 0 <= i < K.size(), local_pool[i] is a list with all the maximal stable sets of V[i] locally found

  BRANCH_STATUS st;   // Branching status

  int branch_vertex_u;       // Vertex u
  int branch_vertex_v;       // Vertex v
  
  std::vector<int>  vertex_mapping; // For each i (original vertex), vertex_mapping[i] is the vertex that stands for i in the current graph (or -1 if it was deleted)
  std::vector<int> precoloring; // For each i (original vertex), precoloring[i] is the color of vertex i (or -1 if it is uncolored)
  double precoloring_value; // Precoloring value

  // Read graph
  void read_graph(char *filename);

  // Read costs
  void read_costs(char *filename);

  // Read lists
  void read_lists(char *filename);

  // Heuristic resolution of the MWSSP in V[i]
  //  We only explore the stable sets in stabs[i]
  //  If the answer is yes, then the result is saved in best_stable and n_best_stable
  bool solve_MWSSP_heuristic (int i, double goal, nodepnt **best_stable, int *n_best_stable);

  // Exact resolution of the MWSSP in V[i]
  //  We apply Sewell's algorithm
  //  If the answer is yes, then the result is saved in best_stable and n_best_stable
  bool solve_MWSSP_exact (int i, double goal, nodepnt **best_stable, int *n_best_stable);

  // Maximize best_stable in V[i] and update n_best_stable
  void maximize_stable_set(int i, nodepnt **best_stable, int *n_best_stable);

  // Internal preprocessor functions
  // preprocess_set_color: set a color for a vertex
  // preprocess_remove_vertex: remove a precolored vertex from the graph and save its color
  void preprocess_set_color(int v, int j);
  void preprocess_remove_vertex(int v, int j);

};

#endif
