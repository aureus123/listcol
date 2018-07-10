#include "stable.h"
#include "io.h"
#include <algorithm>
#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>

#define INTFACTOR 1000.0
#define MAXTIME_MWIS 300.0

// Perform some initializations for the MWSS algorithm, including the generation of the subgraphs
Sewell::Sewell (Graph& G) : Mgraph (G.colors), Minfo (G.colors) {

	default_parameters(&Mparms);
	Mparms.cpu_limit = MAXTIME_MWIS;
	// Mparms.reorder = 1; /* reorder vertices of subgraph by degree */
	// if (density < 0.2) Mparms.clique_cover = 2;
	// else Mparms.clique_cover = 1;

	for (int k = 0; k < G.colors; k++) {
		reset_pointers(&Mgraph[k], &Mdata, &Minfo[k]);

        vector<int> Vk;
        G.get_Vk(k, Vk);
		int vertices_Gk = Vk.size();

		/* generate the graph Gk: recall that the first entry (i.e. array[0]) is not used!!! */
		allocate_graph(&Mgraph[k], vertices_Gk);
		for (int v = 1; v <= vertices_Gk; v++) {
			Mgraph[k].node_list[v].name = v;
			Mgraph[k].node_list[v].degree = 0;
		}
		for (int row = 1; row <= vertices_Gk; row++) {
			int u = Vk[row - 1]; /* recall the "0" is not used! */
			Mgraph[k].adj[row][row] = 0;
			for (int col = row + 1; col <= vertices_Gk; col++) {
				int v = Vk[col - 1]; /* recall the "0" is not used! */
				int val = G.is_edge(u,v) ? 1 : 0;
				Mgraph[k].adj[row][col] = val;
				Mgraph[k].adj[col][row] = val;
			}
		}
		build_graph(&Mgraph[k]);
	    if (initialize_max_wstable(&Mgraph[k], &Minfo[k]) > 0) 
            bye("Failed in initialize_max_wstable");
    }

}

Sewell::~Sewell () {
    for(unsigned int k = 0; k < Mgraph.size(); k++) {
        free_graph(&Mgraph[k]);
        free_wstable_info(&Minfo[k]);
    }
}

void Sewell::solve (int k, vector<double>& pi, double goal, vector<int>& stable_set, double& weight) {

	// Perform optimization
	for (int i = 1; i <= Mgraph[k].n_nodes; i++) 
        Mgraph[k].weight[i] = (int)(INTFACTOR*pi[i - 1]); /* recall that "0" is not used! */
	if (call_max_wstable(&Mgraph[k], &Mdata, &Mparms, &Minfo[k], (int)(INTFACTOR*goal), 0) > 0) 
        bye("Failed in call_max_wstable");
	//	if (Mparms.cpu_limit >= 0 && Minfo.cpu > Mparms.cpu_limit) printf("cpu_limit of %f seconds exceeded: %f seconds. Solution may not optimum.\n", Mparms.cpu_limit, Minfo.cpu);
	// printf("Found best stable set of weight %d.\n", Mdata.best_z);

	// Save best stable set (in terms of the vertices in G)
	int stable_set_size = Mdata.n_best;
	for (int i = 1; i <= stable_set_size; i++) 
        stable_set.push_back(Mdata.best_sol[i]->name - 1); /* recall that "0" is not used! */
	weight = ((double)Mdata.best_z) / INTFACTOR;

    free_data(&Mdata);

	return;

}

/*****************************************************************************************/

CPLEX::CPLEX (Graph& G) : G(G), Xenv(G.colors), Xmodel(G.colors), Xobj(G.colors), Xvars(G.colors) {

    for (int k = 0; k < G.colors; k++) {

        vector<int> Vk;
        G.get_Vk(k,Vk);

        Xmodel[k] = IloModel (Xenv[k]);
        Xobj[k] = IloAdd (Xmodel[k], IloMaximize(Xenv[k],0));
        Xvars[k] = IloNumVarArray (Xenv[k], Vk.size(), 0.0, 1.0, ILOBOOL);

        for (unsigned int i = 0; i < Vk.size() - 1; i++)
                for (unsigned int j = i + 1; j < Vk.size(); j++)  {
                        if (G.is_edge(Vk[i],Vk[j])) {
                            IloExpr restr(Xenv[k]);
                            restr += Xvars[k][i] + Xvars[k][j];
                            Xmodel[k].add(restr <= 1);
                            restr.end();
                        }
                }        
    }

}

CPLEX::~CPLEX () {
    for (int k = 0; k < G.colors; k++) {
        Xobj[k].end();
        Xvars[k].end();
        Xmodel[k].end();
        Xenv[k].end();
    }
}

void CPLEX::solve (int k, vector<double>& pi, double goal, vector<int>& stable_set, double& weight) {
     

    IloCplex Xcplex(Xmodel[k]);

    IloNumArray pii(Xenv[k], pi.size());
    for (unsigned int k = 0; k < pi.size(); k++)
        pii[k] = pi[k];
    Xobj[k].setLinearCoefs(Xvars[k], pii);

    Xcplex.setDefaults();
    Xcplex.setOut(Xenv[k].getNullStream());
    Xcplex.setWarning(Xenv[k].getNullStream());
    Xcplex.extract(Xmodel[k]);
    Xcplex.solve();

    vector<int> Vk;
    G.get_Vk(k,Vk);

    for (unsigned int i = 0; i < Vk.size(); i++)
            if (Xcplex.isExtracted(Xvars[k][i]))
                    if (Xcplex.getValue(Xvars[k][i]) > 0.5)
                            stable_set.push_back(i);
    weight = Xcplex.getObjValue();

    Xcplex.end();
}
                   

/*****************************************************************************************/
