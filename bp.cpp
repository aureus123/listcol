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

#define SHOWINSTANCE
#define SHOWSTATICS
//#define DFS

class Node {

    public:
    Graph* graph;
    double obj_value;
    int ret;
    double prio;

    Node(Graph* G, vector<int>& costs_list) : opt(*G,costs_list) {
        graph = G;
    };
    void solve() {
        ret = opt.optimize(obj_value);
#ifdef DFS
        prio = 0.0;
#else
        prio = obj_value;
#endif
        return;
    };
    void find_branching_vertices(int& u, int& v) {
        opt.find_branching_vertices(u,v);
        return;
    };
    void save_coloring(vector<int>& f) {
        opt.save_coloring(f);
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

// urnd - generate numbers with uniform random distribution
//  if flag=false, the interval is [a, b]
//  if flag=true, the interval is [a, b)
float urnd(float a, float b, bool flag)
{
	return a + rand() * (b - a) / (float)(RAND_MAX + (flag ? 1 : 0));
}

int main (int argc, char **argv) {

#ifndef VERBOSE
	cout.setstate(ios::failbit);
#endif

	set_color(15);
	cout << "LISTCOL - Solves the Minimum Cost List Coloring Problem." << endl;
	set_color(7);

	if (argc < 2 || argc == 3 || argc > 4) {
		cout << "Usage: listcol file.graph [file.cost file.list]" << endl;
		cout << "  if cost and list file are not given, then it solves the classic coloring problem" << endl;
		bye("Bye!");
	}
    
    // Construct graph
    Graph *G;
    vector<int> costs_list;

	if (argc == 4) { 
        // .cost and .list provided by user
        G = new Graph(argv[1],argv[2],costs_list,argv[3]);
	}
	else {
		// list of costs and list of colors are automatically genereted
        G = new Graph(argv[1],costs_list);
    }

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

    // Priority queue
    priority_queue<Node*, vector<Node*>, ComparePrioQueue> queue;
    Node *node = new Node(G, costs_list);
    node->solve();
    queue.push(node);

    // Best integer solution
    double best_integer = DBL_MAX;

    int state = -1;
    int explored_nodes = 0;
    vector<int> f; // best coloring

    while (!queue.empty()) {

        node = queue.top();
        queue.pop();
        explored_nodes++;

        if (node->ret == 1) { // Optimal LP solution is integer
            if (node->obj_value < best_integer) { // Update best integer
                best_integer = node->obj_value;
                node->save_coloring(f);
            }
            delete node->graph; // Prune by optimality
            state = 0;
        }

        else if (node->ret == 0) { // Optimal LP solution is fractional
            if (node->obj_value >= best_integer) {
                delete node->graph; // Prune by bound
                state = 1;
            }
            else if (ceil(node->obj_value) >= best_integer) {
                delete node->graph; // Prune by uper-bound
                state = 2;
            }
            else {

                // Find vertices u and v for branching
                int u, v;
                node->find_branching_vertices(u,v);

                Graph* G1 = node->graph;
                Graph* G2 = new Graph(*G1);
                G2->collapse_vertices(u,v);
                G1->join_vertices(u,v);

                Node *node1 = new Node(G1, costs_list);
                node1->solve();
                Node *node2 = new Node(G2, costs_list);
                node2->solve();

                queue.push(node2);
                queue.push(node1);

                state = 3;
            }
        }

        else if (node->ret == -1) { // LP relaxation is infeasible
            delete node->graph;   // Prune by infeasibility
            state = 4;
        }

        cout << "Obj value = " << setprecision(2) << fixed << node->obj_value << "\t Best int = ";
        if (best_integer == DBL_MAX)
            cout << "inf";
        else
            cout << best_integer;
        //cout << "\t Gap = " << setprecision(2) << fixed << abs(best_bound - best_integer) / (0.0000000001 + abs(best_integer)) * 100 << "%";
        cout << "\t Nodes: expored = " << explored_nodes << ", left = " << queue.size();
        if (state == 0)
            cout << "\t Pruned (opt)" << endl;
        else if (state == 1)
            cout << "\t Pruned (bound)" << endl;
        else if (state == 2)
            cout << "\t Pruned (round-up bound)" << endl;
        else if (state == 3)
            cout << "\t Branched" << endl;
        else if (state == 4)
            cout << "\t Pruned (infeas)" << endl;

        delete node;

    }

    cout << endl << "Optimal coloring information" << endl;
    cout << "cost = " << best_integer << endl;
    for (unsigned int i = 0; i < f.size(); ++i)
        cout << "f(" << i << ") = " << f[i] << endl;

	if (argc == 4) G = new Graph(argv[1],argv[2],costs_list,argv[3]);
	else G = new Graph(argv[1],costs_list);
    G->check_coloring(f); 
    delete G;

	return 0;

}
