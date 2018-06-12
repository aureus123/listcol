#include "io.h"
#include "graph.h"
#include "lp.h"
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <queue>
#include <functional>
#include <cfloat>

#define RANDOMCOSTS
#define SHOWINSTANCE

class ComparePrioQueue
{
public:
    bool operator()(pair<int,reference_wrapper<Graph> > n1,pair<int,reference_wrapper<Graph> > n2) {
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

	// Read graph
    vector<pair<int,int> > edges_list;
	int vertices = read_graph(argv[1], edges_list);

    vector<int> costs_list;
    vector<vector<int> > colors_list;
    int colors;
	if (argc == 4) {
		// Read colors and its costs
		colors = read_cost(argv[2], costs_list);

		// Read list of colors per vertex
		read_list(argv[3], vertices, colors, colors_list);
	}
	else {
		// set C = V, its costs as 1 (or random from {1,...,10}) and L(v)
		colors = vertices;
        costs_list.resize(colors);
		for (int k = 0; k < colors; k++) {
#ifdef RANDOMCOSTS
			costs_list[k] = (int)urnd(1, 10, false);
#else
			costs_list[k] = 1;
#endif
		}
        colors_list.resize(vertices, vector<int> (colors));
		for (int v = 0; v < vertices; v++)
			for (int k = 0; k < colors; k++)
                colors_list[v][k] = k;
	}


    // Construct graph
    Graph G (edges_list, colors_list, costs_list);

#ifdef SHOWINSTANCE

	set_color(2);

	cout << "Neighborhoods:" << endl;
	int maxdelta = 0;
	for (int v = 0; v < vertices; v++) {
		cout << "N(" << v << ") = {";
        vector<int> Nv;
        G.get_Nv(v,Nv);
		int degree = Nv.size();
		if (degree > maxdelta) maxdelta = degree;
		for (int d = 0; d < degree; d++) cout << " " << Nv[d];
		cout << " }, degree = " << degree << endl;
	}
	cout << "Maximum degree = " << maxdelta << endl;

	cout << "Vector of costs: {";
	for (int k = 0; k < colors; k++) {
		cout << " " << k << "->" << costs_list[k];
	}
	cout << " }, colors = " << colors << endl;

	cout << "List of colors:" << endl;
	for (int v = 0; v < vertices; v++) {
		cout << "L(" << v << ") = {";
		for (unsigned int s = 0; s < colors_list[v].size(); s++) cout << " " << colors_list[v][s];
		cout << " }" << endl;
	}

#endif

	/* Show some basic statistics */
	set_color(6);
	cout << "Statistics:" << endl;
    int edges = edges_list.size();
	int clique_size = vertices * (vertices - 1) / 2;
	float density = 100.0 * (float)edges / (float)clique_size;
	cout << "  |V| = " << vertices << ", |E| = " << edges << " (density = " << density << "%), |C| = " << colors << "." << endl;
	/* Average and standard deviation of size of lists */
	float prom = 0.0;
	for (int k = 0; k < vertices; k++) prom += (float)colors_list[k].size();
	prom /= (float)vertices;
	float sigma = 0.0;
	for (int k = 0; k < vertices; k++) {
		float substr = (float)colors_list.size() - prom;
		sigma += substr * substr;
	}
	sigma /= (float)(vertices - 1);
	cout << "  Behaviour of |L(v)| ---> prom = " << prom << ", sigma = " << sqrt(sigma) << "." << endl;
	/* Average and standard deviation of vertices of Gk */
	prom = 0.0;
	for (int k = 0; k < colors; k++) prom += (float)G.get_Vk_size(k);
	prom /= (float)colors;
	sigma = 0.0;
	for (int k = 0; k < colors; k++) {
		float substr = (float)G.get_Vk_size(k) - prom;
		sigma += substr * substr;
	}
	sigma /= (float)(colors - 1);
	cout << "  Behaviour of |V(Gk)| ---> prom = " << prom << ", sigma = " << sqrt(sigma) << "." << endl;
	set_color(7);

#ifndef VERBOSE
	std::cout.clear();
#endif

    // BRANCH AND PRICE
    cout << endl << "Branch and Price" << endl;

    // Priority queue
    priority_queue<pair<int,reference_wrapper<Graph> >, vector<pair<int,reference_wrapper<Graph> > >, ComparePrioQueue> queue;
    queue.push(pair<int,reference_wrapper<Graph> > (1,G));

    // Best integer solution
    double best_integer = DBL_MAX;

    while (!queue.empty()) {

        pair<int,Graph> p = queue.top();
        queue.pop();

        // Solve LP relaxation by column generation approach
        double obj_value;      
        int ret = optimize(p.second, obj_value);

        switch (ret) {
        case 1: // Optimal LP solution is integer
            if (obj_value < best_integer) { // Update best integer
                best_integer = obj_value;
                cout << "New incumbent: " << best_integer << endl;
            }
        case 0: // Optimal LP solution is fractional
            if (obj_value >= best_integer)
                ; // Prune by bound
            else
                ;
        case -1: // LP relaxation is infeasible
            ;   // Prune by infeasibility
        }

    }

	return 0;

}
