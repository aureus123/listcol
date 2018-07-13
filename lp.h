#ifndef _LP_H_
#define _LP_H_

#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>
#include "graph.h"
#include <vector>

enum LP_STATE {INFEASIBLE, INTEGER, FRACTIONAL};

class LP {

    public:

    LP(Graph* G);
    ~LP();

    LP_STATE optimize (double goal);
    void branch1 (vector<LP*>& lps);     // Trick's branching rutine
    void branch2 (vector<LP*>& lps);     // Sewell's branching rutine

    double get_obj_value();
    void save_solution(vector<int>& sol);
    bool check_solution(vector<int>& sol);

    private:
    
    IloEnv Xenv;            // CPLEX environment structure
    IloModel Xmodel;        // CPLEX model
    IloObjective Xobj;      // CPLEX objective function
    IloNumVarArray Xvars;   // CPLEX variables
    IloRangeArray Xrestr;   // CPLEX constraints
    IloNumArray solution;
    double obj_value;

    Graph* G;

    void initialize_LP();
    void set_params(IloCplex&);
    void select_vertices (int& u, int& v);
    void select_vertex (int& v);

};

#endif
