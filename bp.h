#ifndef _BP_H_
#define _BP_H_

#include "lp.h"
#include "graph.h"
#include <list>
#include <vector>
#include <set>

class Coloring {

    public:

    void get_coloring (vector<int>& f) {
        f.clear();
        f.resize(coloring.size());
        f = coloring;
    }
    void get_active_colors (set<int>& f) {
        f.clear();
        f = active_colors;
    }

    void save(LP& lp) {
        coloring.clear();
        active_colors.clear();
        lp.save_solution(coloring, active_colors);
        value = lp.get_obj_value();
        return;
    };

    void check(Graph& g) {
    if (g.check_coloring(coloring))
        cout << "Valid coloring :)" << endl;
    else
        cout << "Invalid coloring :(" << endl;
    };

    private:

    vector<int> coloring;    // coloring: V -> Nat
    double value;
    set<int> active_colors;  // active colors (subset of C)

};

class Node {

    public:

    Node(LP* lp);
    ~Node();

    double get_obj_value() const;
    LP_STATE solve(double start_t, double root_lower_bound = -1.0);
    void branch(vector<Node*>& sons);

    bool operator< (const Node& n) const;

    template <class Solution>
    void save (Solution&);

    private:

    LP* lp;

};

template <class Solution>
class BP {    

    public:

    BP(Solution&, bool DFS = false, bool EARLY_BRANCHING = false);
    void solve(Node* root);

    // Methods for getting variables' value
    int get_nodes();
    double get_gap();
    double get_primal_bound();
    double get_dual_bound(); 
    double get_time();
    int get_opt_flag();

    private:
    
    list<Node*> L;                   // Priority queue

    Solution& best_integer_solution; // Current best integer solution
    double primal_bound;             // Primal bound (given by the best integer solution)
    double dual_bound;               // Dual bound (given by the worst open relaxation)
    int nodes;                       // Number of processed nodes so far. A node is considered processed if its relaxation has been solved
    double start_t;                  // B&P initial execution time
    int opt_flag;                    // Optimality flag
    double time;                     // Total execution time
    double root_lower_bound;         // Round-up of the LP objective value at the root of the B&P. Used for early branching

    void push (Node* node);
    Node* top();
    void pop();
    void update_primal_bound(Node& node);
    double calculate_dual_bound();
    void show_stats(Node& node);

    // Flags
    bool DFS;
    // **** Early Branching **** 
    //  At any node, if the column generation arrives to a solution with lower value than root_lower_bound,
    //  then is not necessary to keep optimizing the node and it could be inmediatly banched.
    //  CAUTION: Early branching may alter the way nodes are branched (beacause braching depends on the LP solution)
    bool EARLY_BRANCHING;

};

template class BP<Coloring>;

#endif
