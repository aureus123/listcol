#include "graph.h"
#include "bp.h"
#include "lp.h"
#include "io.h"
#include <cstring>

// for linux users: do not define VISUALC
#ifndef VISUALC
#include <unistd.h>
#include <sys/times.h>
#else
#include <windows.h>
#endif

//#define SHOWINSTANCE
#define SHOWSTATICS
#define VERBOSE

int main (int argc, char **argv) {

	set_color(15);
	cout << "LISTCOL - Solves the Minimum Cost List Coloring Problem." << endl;
	set_color(7);
    
	// Read instance and construct graph
    cout << "Reading instance " << argv[1] << endl;
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
#ifdef COLORS_DELETION
    G->delete_equal_colors();
#endif
	double stop_t = ECOclock();
	cout << "Time of instance reading = " << stop_t - start_t << " sec." << endl;

#ifndef VERBOSE
	std::cout.clear();
#endif
#ifdef SHOWINSTANCE
    G->show_instance();
#endif
#ifdef SHOWSTATICS
    G->show_statics();
#endif

    // BRANCH AND PRICE
    cout << endl << "Branch and Price" << endl;
    start_t = ECOclock();
    Node* root = new Node(new LP(G));  // Create root
    Coloring col;                      // Create solution
    BP<Coloring> bp(col, true);        // Initialize B&P
    bp.solve(root);                    // Solve B&P
    stop_t = ECOclock();

    // Check solution (G must be re-initialize because it was deleted)
    if (bp.get_primal_bound() != 99999999) {
        G = new Graph(arg1,arg2,costs_list,arg3);
        col.check(*G);
        delete G;
    }

    // Show statics in std output
    cout << "Optimization time = " << stop_t - start_t << "s" << endl;
    cout << "Number of explored nodes = " << bp.get_nodes() << endl;
    if (bp.get_primal_bound() != 99999999 && bp.get_dual_bound() != -99999999 && stop_t - start_t > MAXTIME)
        cout << "Gap = " << bp.get_gap() << " %" << endl;

    // Write output file
    char out[300];
    strncpy(out,filename,300);
    strcat(out,".out");
    ofstream fout (out, std::ios::trunc);
    double dbound, pbound;
    dbound = bp.get_dual_bound();
    pbound = bp.get_primal_bound();
    fout << bp.get_opt_flag() << ":";
    if (dbound == -99999999) { fout << -99999999 << ":"; }
    else { fout << dbound << ":"; }
    if (pbound == 99999999) { fout << 99999999 << ":"; }
    else { fout << pbound << ":"; }    
    fout << bp.get_nodes() << ":" << bp.get_time() << ":";
    set<int> active_colors;
    col.get_active_colors(active_colors);
    fout << active_colors.size() << endl;
    for (int k: active_colors)
        fout << k << endl;
    vector<int> coloring;
    col.get_coloring(coloring);
    for (int i: coloring)
        fout << i << endl;


    fout.close();

	return 0;

}
