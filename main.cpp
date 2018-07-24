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

#define SHOWINSTANCE
#define SHOWSTATICS
#define VERBOSE

//
// ECOclock - measure CPU's time
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
#ifdef COLORS_DELETION
    G->delete_equal_colors();
#endif
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
    start_t = ECOclock();
    Node* root = new Node(new LP(G));  // Create root
    Coloring col;                      // Create solution
    BP<Coloring> bp(col,true);         // Initialize B&P
    bp.solve(root);                    // Solve B&P
    stop_t = ECOclock();

    // Show solution
    col.show();

    // Check solution (G must be re-initialize because it was deleted)
    G = new Graph(arg1,arg2,costs_list,arg3);
    col.check(*G);
    delete G;

    cout << "Optimization time = " << stop_t - start_t << "s" << endl;
    cout << "Number of explored nodes = " << bp.get_processed_nodes() << endl;

	return 0;

}
