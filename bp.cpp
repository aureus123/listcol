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

/* for linux users: do not define VISUALC */
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


/*
 * ECOclock - get a timestamp (in seconds)
 */
double ECOclock() {
#ifndef VISUALC
	/* measure user-time: use it on single-threaded system in Linux (more accurate) */
	struct tms buf;
	times(&buf);
	return ((double)buf.tms_utime) / (double)sysconf(_SC_CLK_TCK);
#else
	/* measure standard wall-clock: use it on Windows */
	return ((double)clock()) / (double)CLOCKS_PER_SEC;
#endif
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

    double start_t = ECOclock();

    // Priority queue
    priority_queue<Node*, vector<Node*>, ComparePrioQueue> queue;
    Node *node = new Node(G, costs_list);
    node->solve(-1.0);
    queue.push(node);

    
    double best_integer = DBL_MAX;  // Objective value of the best integer solution
    int explored_nodes = 0;         // Number of explored nodes
    vector<int> f;                  // Best coloring f: V -> C

    // **** Early Branching **** 
    //  root_lower_bound is the round-up of the LP objective value at the root of the B&P
    //  At any node, if the column generation arrives to a solution with lower value than root_lower_bound,
    //  then is not necessary to keep optimizing the node and it could be inmediatly banched.
    //  CAUTION: Early branching may alter the way nodes are branched (beacause braching depends on the LP solution)
    double root_lower_bound = -1.0;
    //if (node->ret != -1) root_lower_bound = ceil(node->obj_value); // Comment this line to turn off early branching

    int state = -1;                 // Flag

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
                node1->solve(root_lower_bound);
                Node *node2 = new Node(G2, costs_list);
                node2->solve(root_lower_bound);

                queue.push(node1);
                queue.push(node2);

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

    double stop_t = ECOclock();

    cout << endl << "Optimal coloring information" << endl;
    cout << "cost = " << best_integer << endl;
    for (unsigned int i = 0; i < f.size(); ++i)
        cout << "f(" << i << ") = " << f[i] << endl;

	if (argc == 4) G = new Graph(argv[1],argv[2],costs_list,argv[3]);
	else G = new Graph(argv[1],costs_list);
    G->check_coloring(f); 
    delete G;

    cout << "Optimization time = " << stop_t - start_t << "s" << endl;

	return 0;

}
