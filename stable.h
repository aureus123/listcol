#ifndef _STABLE_H_
#define _STABLE_H_

#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>
#include "graph.h"
#include <vector>

#include "mwis_sewell/mwss.h"

template <class SpecificMWSS>
class MWSS {

    public:
    MWSS(Graph& G) : data(G) {};
    void solve(int k, vector<double>& pi, double goal, vector<int>& stable_set, double& weight) {
        data.solve(k,pi,goal,stable_set,weight);
        return;
    };

    private:
    SpecificMWSS data;

};

class Sewell {

    public:
    Sewell(Graph& G);
    ~Sewell();
    void solve(int k, vector<double>& pi, double goal, vector<int>& stable_set, double& weight);

    private:
    vector<MWSSgraph> Mgraph;   // subgraphs of G (one per color) for the MWSS algorithm
    wstable_parameters Mparms;  // parameters for the MWSS algorithm
	MWSSdata Mdata;
	vector<wstable_info> Minfo;

};

class CPLEX {

    public:
    CPLEX(Graph& G);
    ~CPLEX();
    void solve(int k, vector<double>& pi, double goal, vector<int>& stable_set, double& weight);

    private:
    Graph& G;
    vector<IloEnv> Xenv;            // CPLEX environments
    vector<IloModel> Xmodel;        // CPLEX models
    vector<IloObjective> Xobj;      // CPLEX objective functions
    vector<IloNumVarArray> Xvars;   // CPLEX variables

};

#endif
