#ifndef _LP_H_
#define _LP_H_

#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>
#include "graph.h"
#include <vector>

class Lopt {

    public:

    Lopt(Graph& G, vector<int>& cost_list);
    ~Lopt();

    int optimize (double& obj_value);
    void find_branching_vertices (int& u, int& v);
    void save_coloring(vector<int>& f);

    private:
    
    IloEnv Xenv;            // CPLEX environment structure
    IloModel Xmodel;        // CPLEX model
    IloObjective Xobj;      // CPLEX objective function
    IloNumVarArray Xvars;   // CPLEX variables
    IloRangeArray Xrestr;   // CPLEX constraints
    IloNumArray solution;

    Graph &G;
    vector<int>& cost_list;

    void initialize_LP(Graph& G);
    void set_params(IloCplex&);

};

#endif
