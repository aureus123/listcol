#ifndef _LP_H_
#define _LP_H_

#define IL_STD 1 // CPLEX lo pide
#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>
#include "graph.h"
#include <set>

#define EPSILON 0.00001
#define MAXTIME 7200.0
//#define SAVELP "form.lp"
#define FICTIONAL_COST 1000
#define THRESHOLD 0.1

#ifndef INITIAL_COLUMN_STRATEGY
#define INITIAL_COLUMN_STRATEGY 0 // 0: dummy columns
                                  // 1: active father's columns
                                  // 2: stable set covering heuristic
#endif

#ifndef BRANCHING_STRATEGY
#define BRANCHING_STRATEGY 0  // 0: Branching on edges
                              // 1: Branching on colors
#endif

#ifndef VARIABLE_SELECTION
#define VARIABLE_SELECTION 0    // 0: (EDG) Trick's strategy
                                //    (CLR) Greatest neighborhood
                                // 1: (EDG) Smallest list of colors
                                //  : (CLR) Smallest list of colors
#endif

#ifndef N_BRANCHS
#define N_BRANCHS 5
#endif

enum LP_STATE {INFEASIBLE, INTEGER, FRACTIONAL, TIME_OR_MEM_LIMIT};

typedef struct Column {
    bool *stable;   // stable[i] = true sii vertex i is in the stable set
    int n_stable;   // size of the stable set
    int color;
    bool fictional;
    Column (bool *stable, int n_stable, int color, bool fictional) : stable(stable), n_stable(n_stable), color(color), fictional(fictional) {} 
} Column;

class LP {

    public: 

    // Constructor
    LP(Graph *G, LP *father, double start_t);

    // Destructor
    ~LP();

    // Optimize LP
    LP_STATE optimize();

    // Save optimal solution
    void save_solution(std::vector<int> &coloring, std::set<int> &active_colors, double &value);

    // Get objective value
    double get_obj_value();

    // Get number of columns
    int get_n_columns();

    // Branch
    void branch(std::vector<LP *> &);

    private:
    
    double start_t;                 // Starting time of b&P
    IloEnv Xenv;                    // CPLEX environment structure
    IloModel Xmodel;                // CPLEX model
    IloObjective Xobj;              // CPLEX objective function
    IloNumVarArray Xvars;           // CPLEX variables
    IloRangeArray Xrestr;           // CPLEX constraints
    IloArray<Column> vars;          // Internal representation of CPLEX's columns 
    IloNumArray values;             // Value of the optimal solution
    IloNum obj_value;               // Objective value of the optimal solution

    std::list<int> pos_vars;        // List with the indexes of the positive variables
    int most_fract_var;             // Index of the most fractional var with at least two vertices in the stable set

    Graph *G;                        // Graph

    // Reset and initialize LP
    void reset();

    // Initialize LP
    void initialize(LP *father);

    // Fill the LP with a feasible set of columns
    void fill_initial_columns(LP *father);

    // Add the column to CPLEX and to our representation
    void add_real_column (bool *column, int size, int color);
    void add_fictional_column (int v);

    // Internal branch functions
    void branch_on_edges(std::vector<LP *> &ret);
    void branch_on_colors(std::vector<LP *> &ret);
    void branch_on_indistinguishable_colors(std::vector<LP *> &ret);

};

#endif
