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

LP_STATE Node::solve(double root_lower_bound) {
    return lp->optimize1(root_lower_bound);
}

void Node::branch(vector<Node*>& sons) {

    vector<LP*> lps;
    lp->branch2(lps);

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

    best_integer_value = DBL_MAX;
    processed_nodes = 0;
    root_lower_bound = -1;

}

template <class Solution>
void BP<Solution>::solve (Node* root) {

    push(root);

    if (L.empty()) return; // root was pruned

    // Otherwise, start B&P

    root_lower_bound = ceil(root->get_obj_value());

    while (!L.empty()) {

        // Pop
        Node* node = top();
        show_stats(*node);   // First show_stats, then pop
        pop();

        // Add sons
        vector<Node *> sons;
        node->branch(sons);
        for (auto n: sons)
            push(n);
    
        delete node;

    }

}

template <class Solution>
void BP<Solution>::push (Node* node) {
    
    processed_nodes++;

    // Solve the linear relaxation of the node and prune if possible
    LP_STATE state = EARLY_BRANCHING ? node->solve(root_lower_bound) : node->solve(); 
    double obj_value;
cout << node->get_obj_value() << endl;
    switch (state) {

        case INFEASIBLE:
            // Prune by infeasibility
            delete node;
            return;

        case INTEGER:
            obj_value = node->get_obj_value();
            if (obj_value < best_integer_value) update_best_integer(*node);
            // Prune by optimality
            // The node is deleted after being overwrited by a beter node
            return;

        case FRACTIONAL:
            obj_value = node->get_obj_value();
            if (ceil(obj_value) >= best_integer_value) {
                // Prune by bound
                delete node;
                return;
            }
            break;

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
void BP<Solution>::update_best_integer(Node& node) {

    // Update best integer solution and value
    node.save(best_integer_solution);
    best_integer_value = round(node.get_obj_value());

    // Prune if possible
    for (auto it = L.begin(); it != L.end(); ) {
        if ((*it)->get_obj_value() >= best_integer_value)
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
void BP<Solution>::show_stats (Node& node) {

    // Calculate GAP
    double best_bound = DBL_MAX;  // minimum objective value of unpruned nodes
    if (DFS) {
        // Traverse the list and search for the minimum objective value
        for (auto it = L.begin(); it != L.end(); ++it)
            if ((*it)->get_obj_value() < best_bound)
                best_bound = (*it)->get_obj_value();
    }
    else 
        best_bound = L.front()->get_obj_value();
    double gap = abs(best_bound - best_integer_value) / (0.0000000001 + abs(best_integer_value)) * 100;

    cout << fixed << setprecision(2);

    cout << "Obj value = " << node.get_obj_value() << "\t Best int = ";
    if (best_integer_value == DBL_MAX)
        cout << "inf";
    else
        cout << (int) best_integer_value;
    cout << "\t Gap = " << gap << "%";
    cout << "\t Nodes: processed = " << processed_nodes << ", left = " << L.size() << endl;

}

template <class Solution>
int BP<Solution>::get_processed_nodes() {
    return processed_nodes;
}

/*

#include "io.h"
#include "graph.h"
#include "lp.h"
#include <iostream>
#include <math.h>
#include <queue>
#include <functional>
#include <cfloat>
#include <iomanip> 
#include <ctime>
#include <cmath>  

// for linux users: do not define VISUALC
#ifndef VISUALC
#include <unistd.h>
#include <sys/times.h>
#else
#include <windows.h>
#endif

#define SHOWINSTANCE
#define SHOWSTATICS
//#define DFS

class Node {

    public:
    Graph* graph;               // Graph
    double obj_value;           // Objective value of LP
    int ret;                    // Exit state of the LP (Optimal: fractional or integer, infeasible)
    double prio;                // Priority in the priority queue (equal to 0 in DFS and obj_value in best-bound strategy)

    Node(Graph* G, vector<int>& costs_list) : opt(*G,costs_list) {
        graph = G;
    };
    void solve(double goal) {
        ret = opt.optimize(goal, obj_value);
#ifdef DFS
        prio = 0.0;
#else
        prio = obj_value;
#endif
        return;
    };
    void select_vertices_trick(int& u, int& v) {
        opt.select_vertices_trick(u,v);
        return;
    };
    void select_vertex_sewell(int& v) {
        opt.select_vertex_sewell(v);
        return;
    };
    void save_coloring(vector<int>& f) {
        opt.save_coloring(f);
        return;
    };
    void show_solution() {
        opt.show_solution();
        return;
    };

    private:
    Lopt opt;

};

class ComparePrioQueue {
public:
    bool operator()(Node* n1, Node* n2) {
        return n1->prio <= n2->prio;
    }
};

//
// ECOclock - get a timestamp (in seconds)
//
double ECOclock() {
#ifndef VISUALC
	// measure user-time: use it on single-threaded system in Linux (more accurate)
	struct tms buf;
	times(&buf);
	return ((double)buf.tms_utime) / (double)sysconf(_SC_CLK_TCK);
#else
	// measure standard wall-clock: use it on Windows 
	return ((double)clock()) / (double)CLOCKS_PER_SEC;
#endif
}

void handle_node(Node* node, priority_queue<Node*, vector<Node*>, ComparePrioQueue>& queue, 
double& best_integer, vector<int>& f) {

    if (node->ret == -1) { // LP relaxation is infeasible
        delete node->graph;   // Prune by infeasibility
        delete node;
    }
    else if (node->ret == 1) { // Optimal LP solution is integer
        if (node->obj_value < best_integer) { // Update best integer
            best_integer = round(node->obj_value);
            node->save_coloring(f);
        }
        delete node->graph; // Prune by optimality
        delete node;
    }
    else if (node->ret == 0) { // Optimal LP solution is fractional
        if (node->obj_value >= best_integer) {
            delete node->graph; // Prune by bound
            delete node;
        }
        else if (ceil(node->obj_value) >= best_integer) {
            delete node->graph; // Prune by uper-bound
            delete node;
        }
        else queue.push(node);
    }

    return;
}

void branch_trick(Node *node, priority_queue<Node*, vector<Node*>, ComparePrioQueue>& queue, 
vector<int> costs_list, double& best_integer, vector<int>& f, double root_lower_bound) {

    // Find vertices u and v for branching
    int u, v;
    node->select_vertices_trick(u,v);

    // Create graphs
    Graph* G1 = node->graph;
    Graph* G2 = new Graph(*G1);
    G2->collapse_vertices(u,v);
    G1->join_vertices(u,v);

    // Create nodes
    Node *node1 = new Node(G1, costs_list);
    node1->solve(root_lower_bound);
    Node *node2 = new Node(G2, costs_list);
    node2->solve(root_lower_bound);

    // Handle nodes
    handle_node(node1, queue, best_integer, f);
    handle_node(node2, queue, best_integer, f);

    return;
}

void branch_sewell(Node *node, priority_queue<Node*, vector<Node*>, ComparePrioQueue>& queue, 
vector<int> costs_list, double& best_integer, vector<int>& f, double root_lower_bound) {

    // Select a vertex to branch
    int v;
    node->select_vertex_sewell(v);

    // Create nodes
    vector<int> Lv;
    node->graph->get_Lv(v, Lv);
    for (int k: Lv) {
        Graph *G1 = new Graph(*(node->graph));
        G1->color_vertex(v,k);
        Node *node1 = new Node(G1, costs_list);
        node1->solve(root_lower_bound);
        handle_node(node1, queue, best_integer, f);
    }

    delete node->graph;

    return;
}

int main (int argc, char **argv) {

	set_color(15);
	cout << "LISTCOL - Solves the Minimum Cost List Coloring Problem." << endl;
	set_color(7);
    
	// Read instance and construct graph
	double start_t = ECOclock();
    Graph *G;
    vector<int> costs_list;
	char* filename = argv[1];
    char arg1[300], arg2[300], arg3[300];
    strncpy(arg1,filename,300);
    strcat(arg1,".graph");
    strncpy(arg2,filename,300);
    strcat(arg2,".cost");
    strncpy(arg3,filename,300);
    strcat(arg3,".list");
    G = new Graph(arg1,arg2,costs_list,arg3);
	double stop_t = ECOclock();
	cout << "Time of instance reading = " << stop_t - start_t << " sec." << endl;

#ifndef VERBOSE
	std::cout.clear();
#endif
#ifdef SHOWINSTANCE
    G->show_instance(costs_list);
#endif
#ifdef SHOWSTATICS
    G->show_statics();
#endif

    // BRANCH AND PRICE
    cout << endl << "Branch and Price" << endl;
#ifdef DFS
    cout << "Node selection strategy: DFS" << endl << endl;
#else
    cout << "Node selection strategy: best-bound" << endl << endl;
#endif

    start_t = ECOclock();

    // Priority queue
    priority_queue<Node*, vector<Node*>, ComparePrioQueue> queue;
    double best_integer = DBL_MAX;  // Objective value of the best integer solution
    int explored_nodes = 0;         // Number of explored nodes
    vector<int> f;                  // Best coloring f: V -> C

    // Create initial node and solve it
    Node *node = new Node(G, costs_list);
    node->solve(-1.0);

    // Is integer or infeasible?
    if (node->ret == 1) { // Optimal LP solution is integer
        if (node->obj_value < best_integer) { // Update best integer
            best_integer = round(node->obj_value);
            node->save_coloring(f);
        }
        delete node->graph; // Prune by optimality
        delete node;
    }
    else if (node->ret == -1) { // LP relaxation is infeasible
        delete node->graph;   // Prune by infeasibility
        delete node;
    }

    else {

        queue.push(node);

        // **** Early Branching **** 
        //  root_lower_bound is the round-up of the LP objective value at the root of the B&P
        //  At any node, if the column generation arrives to a solution with lower value than root_lower_bound,
        //  then is not necessary to keep optimizing the node and it could be inmediatly banched.
        //  CAUTION: Early branching may alter the way nodes are branched (beacause braching depends on the LP solution)
        double root_lower_bound = ceil(node->obj_value);
        root_lower_bound = -1; // Comment this line to turn on early branching

        while (!queue.empty()) {

            node = queue.top();
            queue.pop();

            // We already know that the node is fractional, but maybe it could be pruned by bound
            if (node->obj_value >= best_integer) {
                delete node->graph; // Prune by bound
                delete node;
                continue;
            }
            else if (ceil(node->obj_value) >= best_integer) {
                delete node->graph; // Prune by round-up bound
                delete node;
                continue;
            }

            explored_nodes++;

            branch_trick(node, queue, costs_list, best_integer, f, root_lower_bound);
            //branch_sewell(node, queue, costs_list, best_integer, f, root_lower_bound);

            cout << "Obj value = " << fixed << setprecision(2) << node->obj_value << "\t Best int = ";
            if (best_integer == DBL_MAX)
                cout << "inf";
            else
                cout << (int) best_integer;
            //cout << "\t Gap = " << setprecision(2) << fixed << abs(best_bound - best_integer) / (0.0000000001 + abs(best_integer)) * 100 << "%";
            cout << "\t\t Nodes: expored = " << explored_nodes << ", left = " << queue.size() << endl;

            delete node;

        }
    }

    stop_t = ECOclock();

    cout << endl << "Optimal coloring information" << endl;
    cout << "cost = " << best_integer << endl;
    for (unsigned int i = 0; i < f.size(); ++i)
        cout << "f(" << i << ") = " << f[i] << endl;

	G = new Graph(arg1,arg2,costs_list,arg3);
    G->check_coloring(f); 
    delete G;

    cout << "Number of explored nodes = " << explored_nodes << endl;
    cout << "Optimization time = " << stop_t - start_t << "s" << endl;

	return 0;

}

*/
