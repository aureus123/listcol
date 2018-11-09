#ifndef _BP_H_
#define _BP_H_

#include "lp.h"
#include "graph.h"
#include <list>
#include <vector>


class Coloring {

    public:

    void show() {
        cout << endl << "Optimal coloring information:" << endl;
        //for (unsigned int i = 0; i < f.size(); ++i)
        //    cout << "f(" << i << ") = " << f[i] << endl;
        cout << "cost = " << value << endl;


        // BORRAR
        ofstream out ("a.txt", std::ios::app);
        out << value << "\t";
        out.close();


        return;
    };
    void save(LP& lp) {
        f.clear();
        lp.save_solution(f);
        value = lp.get_obj_value();
        return;
    };
    void check(Graph& g) {
    if (g.check_coloring(f))
        cout << "Valid coloring :)" << endl;
    else
        cout << "Invalid coloring :(" << endl;
    };

    private:

    vector<int> f;
    double value;

};

class Node {

    public:

    Node(LP* lp);
    ~Node();

    double get_obj_value() const;
    LP_STATE solve(double root_lower_bound = -1.0);
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

    int get_processed_nodes();

    private:
    
    list<Node*> L;                   // Priority queue

    Solution& best_integer_solution; // Current best integer solution
    double best_integer_value;       // Value of the current best integer solution
    int processed_nodes;             // Number of processed nodes so far. A node is considered processed if its relaxation has been solved
    double root_lower_bound;         // Round-up of the LP objective value at the root of the B&P. Used for early branching

    void push (Node* node);
    void update_best_integer(Node& node);
    Node* top();
    void pop();
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
