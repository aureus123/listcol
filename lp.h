#ifndef _LP_H_
#define _LP_H_

#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>
#include "graph.h"
#include <vector>
#include <list>
#include <set>

#define EPSILON 0.00001
#define MAXTIME 7200.0
//#define SAVELP "form.lp"
//#define PUREBYB
#define NOMEMEMPHASIS
#define FICTIONAL_COST 1000
#define THRESHOLD 0.1
//#define INITIAL_HEURISTIC

enum LP_STATE {INFEASIBLE, INTEGER, FRACTIONAL};

struct var {
    vector<int> stable;
    int color;
};

class LP {

    public:

    LP(Graph* G);
    ~LP();

    LP_STATE optimize (double goal);
    void branch (vector<LP*>& lps);     // Trick's branching rutine

    double get_obj_value();
    void save_solution(vector<int>& sol);
    bool check_solution(vector<int>& sol);

    private:
    
    IloEnv Xenv;                    // CPLEX environment structure
    IloModel Xmodel;                // CPLEX model
    IloObjective Xobj;              // CPLEX objective function
    IloNumVarArray Xvars;           // CPLEX variables
    IloRangeArray Xrestr;           // CPLEX constraints
    IloNumArray values;             // Values of Xvars
    list<var> vars;                 // User representation of the LP
    list<var>::iterator it_branch;  // Iterator to the branching variable

    double obj_value;
    bool fictional;

    Graph* G;

    void initialize_LP();
    void set_params(IloCplex&);
    void select_vertices (int& u, int& v);
    void select_vertex (int& v, set<int>& colors);

};

#endif
