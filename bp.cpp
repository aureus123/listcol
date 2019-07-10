#include "bp.h"
#include "io.h"
#include <iterator>
#include <cfloat>
#include <cmath>

Node::Node (LP* lp) : lp(lp) {
}

Node::~Node () {
    delete lp;
}

double Node::get_obj_value() const {
    return lp->get_obj_value();
}

bool Node::operator< (const Node& n) const {
    return (get_obj_value() < n.get_obj_value());
}

LP_STATE Node::solve(double start_t, double root_lower_bound) {
    return lp->optimize(start_t, root_lower_bound);
}

void Node::branch(vector<Node*>& sons) {

    vector<LP*> lps;
    lp->branch(lps);

    sons.reserve(lps.size());
    for (auto x: lps)
        sons.push_back(new Node(x));

    return;
}

template <class Solution>
void Node::save(Solution& sol) {
    sol.save(*lp);
}

template <class Solution>
BP<Solution>::BP (Solution& sol, bool DFS, bool EARLY_BRANCHING) : 
best_integer_solution(sol), DFS(DFS), EARLY_BRANCHING(EARLY_BRANCHING) {

    primal_bound = DBL_MAX;
    nodes = 0;
    root_lower_bound = -1;
    start_t = ECOclock();

}

template <class Solution>
void BP<Solution>::solve (Node* root) {

    push(root);

    if (!L.empty()) root_lower_bound = ceil(root->get_obj_value());

    while (!L.empty()) {

        // Pop
        Node* node = top();
        show_stats(*node);   // First show_stats, then pop
        pop();

        // Add sons
        vector<Node *> sons;
        node->branch(sons);
        for (auto n: sons) {
            try {
                push(n);
            }
            catch (...) {
                // Time expired
                opt_flag = 0;
                primal_bound = primal_bound == DBL_MAX ? 99999999 : primal_bound;
                dual_bound = calculate_dual_bound();
                dual_bound = dual_bound == -DBL_MAX ? -99999999 : dual_bound;
                nodes = -1;
                time = MAXTIME;
                L.empty();
                delete node;
                return;
            }
        }
    
        delete node;

    }

    if (primal_bound == DBL_MAX) {     // Infeasibility case:
        opt_flag = 2;
        primal_bound = 99999999;
        double db = calculate_dual_bound();
        dual_bound = db == -DBL_MAX ? -99999999 : db; 
        time = ECOclock() - start_t;
    }
    else {    // Optimality case:
        opt_flag = 1;
        dual_bound = primal_bound;
        time = ECOclock() - start_t;
    }


    return;

}

template <class Solution>
void BP<Solution>::push (Node* node) {
    
    // Solve the linear relaxation of the node and prune if possible
    LP_STATE state = EARLY_BRANCHING ? node->solve(start_t, root_lower_bound) : node->solve(start_t); 
    
    nodes++;    
    double obj_value;

    switch (state) {

        case INFEASIBLE:
            // Prune by infeasibility
            delete node;
            return;

        case INTEGER:
            obj_value = node->get_obj_value();
            if (obj_value < primal_bound) update_primal_bound(*node);
            // Prune by optimality
            // The node is deleted after being overwrited by a beter node
            return;

        case FRACTIONAL:
            obj_value = node->get_obj_value();
            if (ceil(obj_value) >= primal_bound) {
                // Prune by bound
                delete node;
                return;
            }
            break;

        case TIME_LIMIT:
            throw std::exception{};
            return;

        case OTHER:
            bye("Unknown LP status");
            return;

    }

    // Place the node in the list according to its priority

    if (DFS) { // DFS strategy
        L.push_back(node);
        return;
    }

    // Otherwise, best-bound strategy

    // list is empty
    if (L.empty()) {
        L.push_back(node);
        return;
    }

    // list is not empty
    for (auto it = L.begin(); it != L.end(); ++it)
        if (*node < **it) {
            L.insert(it, node);
            return;
        }

    // Otherwise, push back
    L.push_back(node);
    return;

}

template <class Solution>
void BP<Solution>::update_primal_bound(Node& node) {

    // Update best integer solution and value
    node.save(best_integer_solution);
    primal_bound = round(node.get_obj_value());

    // Prune if possible
    for (auto it = L.begin(); it != L.end(); ) {
        if ((*it)->get_obj_value() >= primal_bound)
            it = L.erase(it);
        else
            ++it;
    }

}

template <class Solution>
Node* BP<Solution>::top() {
    return L.back();
}

template <class Solution>
void BP<Solution>::pop() {
    L.pop_back();
    return;
}

template <class Solution>
int BP<Solution>::get_nodes() {
    return nodes;
}

template <class Solution>
double BP<Solution>::get_gap() {
    double _dual_bound = calculate_dual_bound();
    double gap = abs(_dual_bound - primal_bound) / (0.0000000001 + abs(primal_bound)) * 100;
    return gap;
}

template <class Solution>
double BP<Solution>:: get_primal_bound() {
    return primal_bound;
}

template <class Solution>
double BP<Solution>:: get_dual_bound() {
    return dual_bound;
};

template <class Solution>
double BP<Solution>:: get_time() {
    return time;
}

template <class Solution>
int BP<Solution>:: get_opt_flag() {
    return opt_flag;
}

template <class Solution>
double BP<Solution>:: calculate_dual_bound() {

    double _dual_bound = DBL_MAX;  // minimum objective value of unpruned nodes
    if (DFS) {
        // Traverse the list and search for the minimum objective value
        for (auto it = L.begin(); it != L.end(); ++it)
            if ((*it)->get_obj_value() < _dual_bound)
                _dual_bound = (*it)->get_obj_value();
    }
    else {
        _dual_bound = L.front()->get_obj_value();
    }
    

    if (_dual_bound == DBL_MAX) {
        return -DBL_MAX;
    }
    else {
        return _dual_bound;
    }

}

template <class Solution>
void BP<Solution>::show_stats (Node& node) {
	static double first_t = start_t;
	static bool first_call = true;

	double now_t = ECOclock();
	if (first_call) first_call = false;
	else {
		if (now_t - first_t < 10.0) return;
		first_t = now_t;
	}
    // Calculate GAP (it is time cosuming when DFS is used)
    //double _dual_bound = calculate_dual_bound();
    //double gap = abs(_dual_bound - primal_bound) / (0.0000000001 + abs(primal_bound)) * 100;

    cout << fixed << setprecision(2);

    cout << "Obj value = " << node.get_obj_value() << "\t Best int = ";
    if (primal_bound == DBL_MAX)
        cout << "inf";
    else
        cout << (int) primal_bound;
    //cout << "\t Gap = " << gap << "%";
    cout << "\t Nodes: processed = " << nodes << ", left = " << L.size() << "\t time = " << now_t - start_t << endl;

}
