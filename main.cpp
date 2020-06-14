#include "graph.h"
#include "bp.h"
#include "lp.h"
#include "io.h"
#include <cstring>
#include <fstream>

//#define SHOWINSTANCE
//#define SHOWSTATICS
#define VERBOSE

int main (int argc, char **argv) {

    set_color(15);
    std::cout << "LISTCOL - Solves the Minimum Cost List Coloring Problem." << std::endl;
    set_color(5);
    
    // Read instance and construct graph
    std::cout << "Reading instance " << argv[1] << std::endl;
    double start_t = ECOclock();
    std::vector<int> costs_list;
    char* filename = argv[1];
    char arg1[300], arg2[300], arg3[300];
    strncpy(arg1,filename,300);
    strcat(arg1,".graph");
    strncpy(arg2,filename,300);
    strcat(arg2,".cost");
    strncpy(arg3,filename,300);
    strcat(arg3,".list");
    double stop_t = ECOclock();

/*
	std::cout << "Time of instance reading = " << stop_t - start_t << " sec." << endl;
*/

#ifndef VERBOSE
	std::cout.clear();
#endif

/*
#ifdef SHOWINSTANCE
    G->show_instance();
#endif
#ifdef SHOWSTATICS
    G->show_statics();
#endif
*/

    // BRANCH AND PRICE
    set_color(7);
    std::cout << std::endl << "Branch and Price" << std::endl;
    Graph *G = new Graph(arg1,arg2,arg3);   // Create graph
    Node* root = new Node(new LP(G, NULL));  // Create root
    Coloring col;                      // Create solution
    BP<Coloring> bp(col, true);        // Initialize B&P
    bp.solve(root);                    // Solve B&P

/*
    // Check solution (G must be re-initialize because it was deleted)
    if (bp.get_primal_bound() != 99999999) {
        G = new Graph(arg1,arg2,costs_list,arg3);
        col.check(*G);
        delete G;
    }
*/

    // Show statics in std output
    set_color(3);
    std::cout << "Optimization time = " << bp.get_time() << "s" << std::endl;
    std::cout << "Number of explored nodes = " << bp.get_nodes() << std::endl;
    if (bp.get_primal_bound() != 99999999 && bp.get_dual_bound() != -99999999 && stop_t - start_t > MAXTIME)
        std::cout << "Gap = " << bp.get_gap() << " %" << std::endl;

    // Write output file
    char out[300];
    strncpy(out,filename,300);
    strcat(out,".out");
    std::ofstream fout (out, std::ios::trunc);
    double dbound, pbound;
    dbound = bp.get_dual_bound();
    pbound = bp.get_primal_bound();
    fout << bp.get_opt_flag() << ":";
    if (dbound == -99999999) { fout << -99999999 << ":"; }
    else { fout << dbound << ":"; }
    if (pbound == 99999999) { fout << 99999999 << ":"; }
    else { fout << pbound << ":"; }    
    fout << bp.get_nodes() << ":" << bp.get_time() << ":";
    std::set<int> active_colors;
    col.get_active_colors(active_colors);
    fout << active_colors.size() << std::endl;
    for (int k: active_colors)
        fout << k << std::endl;
    std::vector<int> coloring;
    col.get_coloring(coloring);
    for (int i: coloring)
        fout << i << std::endl;

    fout.close();

    return 0;
}
