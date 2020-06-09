#include "bp.h"
#include "io.h"
#include "lp.h"
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

LP_STATE Node::solve(double start_t) {
    return lp->optimize(start_t);
}

void Node::branch(std::vector<Node*>& sons) {

    std::vector<LP*> lps;
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

int Node::get_n_columns() {
    return lp->get_n_columns();    
}

template <class Solution>
BP<Solution>::BP (Solution& sol, bool DFS, bool EARLY_BRANCHING) : 
best_integer_solution(sol), DFS(DFS), EARLY_BRANCHING(EARLY_BRANCHING) {

    primal_bound = DBL_MAX;
    nodes = 0;
    root_lower_bound = -1;

}

template <class Solution>
void BP<Solution>::solve (Node* root) {

    start_t = ECOclock();

    try {
        push(root);
    }
    catch(...) { // Time or mem expired
        opt_flag = 0;
        primal_bound = 99999999;
        dual_bound = -99999999;
        nodes = -1;
        time = ECOclock() - start_t;
        if (time >= MAXTIME) {
            std::cout << "Time limit reached" << std::endl;
            time = MAXTIME; 
        }
        else {
            std::cout << "Mem limit reached" << std::endl;
        }
        return;
    }

#ifdef ONLY_RELAXATION
    if (!L.empty()) {
        // The initial relaxation is fractional
        opt_flag = 3;
        primal_bound = 99999999;
        double db = calculate_dual_bound();
        dual_bound = db == -DBL_MAX ? -99999999 : db; 
        time = ECOclock() - start_t;
        std::cout << "The initial relaxation is fractional" << std::endl;
        std::cout << "Objective value = " << dual_bound << std::endl;
        return;
    }
    else {
        if (primal_bound == DBL_MAX) {
            // The initial relaxation is infeasible
            opt_flag = 2;
            primal_bound = 99999999;
            dual_bound = -99999999; 
            time = ECOclock() - start_t;
            std::cout << "The initial relaxation is infeasible" << std::endl;
            return;        
        }
        else {
            // The initial relaxation is integer
            opt_flag = 1;
            dual_bound = primal_bound; 
            time = ECOclock() - start_t;
            std::cout << "The initial relaxation is integer" << std::endl;
            std::cout << "Objective value = " << dual_bound << std::endl;
            return;          
        }
    }
#endif

    if (!L.empty()) root_lower_bound = ceil(root->get_obj_value());

    while (!L.empty()) {

        // Pop
        Node* node = top();
        show_stats(*node);   // First show_stats, then pop
        pop();

        // Re-try to prune by bound, since primal_bound could have been improved
        if (ceil(node->get_obj_value()) >= primal_bound) {
            delete node;
            continue;
        }

        // Add sons
        std::vector<Node *> sons;
        node->branch(sons);
        for (auto n: sons) {
            try {
                push(n);
            }
            catch (...) {
                // Time or mem expired
                opt_flag = 0;
                primal_bound = primal_bound == DBL_MAX ? 99999999 : primal_bound;
                dual_bound = calculate_dual_bound();
                dual_bound = dual_bound == -DBL_MAX ? -99999999 : dual_bound;
                nodes = -1;
                time = ECOclock() - start_t;
                if (time >= MAXTIME) {
                    std::cout << "Time limit reached" << std::endl;
                    time = MAXTIME; 
                }
                else {
                    std::cout << "Mem limit reached" << std::endl;
                }
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
        std::cout << "Infeasibility proved" << std::endl;
    }
    else {    // Optimality case:
        opt_flag = 1;
        dual_bound = primal_bound;
        time = ECOclock() - start_t;
        std::cout << "Optimality reached" << std::endl;
        std::cout << "Optimal value = " << primal_bound << std::endl;
    }


    return;

}

template <class Solution>
void BP<Solution>::push (Node* node) {
    
    // Solve the linear relaxation of the node and prune if possible
    LP_STATE state = node->solve(start_t); 


#ifdef ONLY_RELAXATION
    std::cout << "Number of columns = " << node->get_n_columns() << std::endl;
#endif

    nodes++;    
    double obj_value;

    switch (state) {

        case INFEASIBLE:
            // Prune by infeasibility
            delete node;
            return;

        case INTEGER:
            // Prune by optimality
            obj_value = node->get_obj_value();
            if (obj_value < primal_bound) update_primal_bound(*node);
            delete node;
            return;

        case FRACTIONAL:
            obj_value = node->get_obj_value();
            if (ceil(obj_value) >= primal_bound) {
                // Prune by bound
                delete node;
                return;
            }
            break;

        case TIME_OR_MEM_LIMIT:
            delete node;
            throw std::exception{};
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
    primal_bound = node.get_obj_value();

    // Prune if possible
    for (auto it = L.begin(); it != L.end(); ) {
        if ((*it)->get_obj_value() >= primal_bound) {
            delete *it;
            it = L.erase(it);
        }
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

    std::cout << std::fixed << std::setprecision(2);

    std::cout << "Obj value = " << node.get_obj_value() << "\t Best int = ";
    if (primal_bound == DBL_MAX)
        std::cout << "inf";
    else
        std::cout << (int) primal_bound;
    //cout << "\t Gap = " << gap << "%";
    std::cout << "\t Nodes: processed = " << nodes << ", left = " << L.size() << "\t time = " << now_t - start_t << std::endl;

}
