#include "io.h"
#include "graph.h"
#include "lp.h"
#include <iostream>
#include <math.h>
#include <queue>
#include <functional>
#include <cfloat>

#define SHOWINSTANCE
#define SHOWSTATICS

class ComparePrioQueue
{
public:
    bool operator()(pair<int,Graph*> n1, pair<int,Graph*> n2) {
        return n1.first<n2.first;
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

    // Priority queue
    priority_queue<pair<int,Graph*>, vector<pair<int,Graph*> >, ComparePrioQueue> queue;
    queue.push(pair<int,Graph*> (1,G));

    // Best integer solution
    double best_integer = DBL_MAX;

    while (!queue.empty()) {

        pair<int,Graph*> p = queue.top();
        queue.pop();

        // Solve LP relaxation by column generation approach
        double obj_value;
        Lopt opt (*p.second, costs_list);
        int ret = opt.optimize(obj_value);
        cout << "Objective value = " << obj_value << "\t Best integer = " << best_integer;

        if (ret == 1) { // Optimal LP solution is integer

            if (obj_value < best_integer) { // Update best integer
                best_integer = obj_value;
                cout << "\t *New incumbent = " << best_integer;
            }
            //delete[] p.second; // Prune by optimality
        }

        else if (ret == 0) { // Optimal LP solution is fractional

            if (obj_value >= best_integer)
                ;//delete[] p.second; // Prune by bound
            else {

                // Find vertices u and v for branching
                int u, v;
                opt.find_branching_vertices(u,v);

                Graph* G1 = p.second;
                Graph* G2 = new Graph(*G1);
                G2->collapse_vertices(u,v); // CAUTION: collapse precedes joint
                G1->join_vertices(u,v);

                queue.push(pair<int,Graph*> (1,G1));
                queue.push(pair<int,Graph*> (1,G2));

            }
        }

        else if (ret == -1) { // LP relaxation is infeasible
            ;
            //delete[] p.second;   // Prune by infeasibility
        }

        cout << "\t Nodes left = " << queue.size() << endl;

    }

	return 0;

}
