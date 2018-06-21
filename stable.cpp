#include "stable.h"
#include "io.h"
#include <algorithm>

#define INTFACTOR 10000.0
#define MAXTIME_MWIS 300.0

// Perform some initializations for the MWSS algorithm, including the generation of the subgraphs
Sewell::Sewell (Graph& G) : Mgraph (G.colors) {

	MWSSdata Mdata;
	wstable_info Minfo;

	default_parameters(&Mparms);
	Mparms.cpu_limit = MAXTIME_MWIS;
	// Mparms.reorder = 1; /* reorder vertices of subgraph by degree */
	// if (density < 0.2) Mparms.clique_cover = 2;
	// else Mparms.clique_cover = 1;

	for (int k = 0; k < G.colors; k++) {
		reset_pointers(&Mgraph[k], &Mdata, &Minfo);

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
    }

}

void Sewell::solve (int k, vector<double>& pi, double goal, vector<int>& stable_set, double& weight) {

	MWSSdata Mdata;
	wstable_info Minfo;

	// Perform optimization
	for (int i = 1; i <= Mgraph[k].n_nodes; i++) 
        Mgraph[k].weight[i] = (int)(INTFACTOR*pi[i - 1]); /* recall that "0" is not used! */
	if (initialize_max_wstable(&Mgraph[k], &Minfo) > 0) 
        bye("Failed in initialize_max_wstable");
	if (call_max_wstable(&Mgraph[k], &Mdata, &Mparms, &Minfo, (int)(INTFACTOR*goal), 0) > 0) 
        bye("Failed in call_max_wstable");
	//	if (Mparms.cpu_limit >= 0 && Minfo.cpu > Mparms.cpu_limit) printf("cpu_limit of %f seconds exceeded: %f seconds. Solution may not optimum.\n", Mparms.cpu_limit, Minfo.cpu);
	// printf("Found best stable set of weight %d.\n", Mdata.best_z);

	// Save best stable set (in terms of the vertices in G)
	int stable_set_size = Mdata.n_best;
	for (int i = 1; i <= stable_set_size; i++) 
        stable_set.push_back(Mdata.best_sol[i]->name - 1); /* recall that "0" is not used! */
	weight = ((double)Mdata.best_z) / INTFACTOR;

	free_data(&Mdata);
	free_wstable_info(&Minfo);

	return;

}

/*****************************************************************************************/
//  RESOLVER MWSS CON CPLEX
/*         
    IloEnv XXenv; // CPLEX environment structure
    IloModel XXmodel(XXenv); // CPLEX model
    IloNumVarArray XXvars(XXenv); // CPLEX variables

    for (int i = 0; i < C_size[k]; i++)
        XXvars.add(IloNumVar(XXenv, 0.0, 1.0, ILOBOOL));

    IloExpr fobj(XXenv, 0);
    for (int i = 0; i < C_size[k]; i++)
            fobj += XXvars[i]*duals[C_set[k][i]];
    XXmodel.add(IloMaximize(XXenv,fobj));
    fobj.end();

    for (int i = 0; i < C_size[k] - 1; i++)
            for (int j = i + 1; j < C_size[k]; j++)  {
                    if (adjacency[C_set[k][i]][C_set[k][j]] > 0) {
                    IloExpr restr(XXenv);
                            restr += XXvars[i] + XXvars[j];
                    XXmodel.add(restr <= 1);
                            restr.end();
                    }
            }

    IloCplex XXcplex(XXmodel);
    XXcplex.setDefaults();
    XXcplex.setOut(XXenv.getNullStream());
    XXcplex.setWarning(XXenv.getNullStream());
    XXcplex.extract(XXmodel);
    XXcplex.solve();

    int counter2 = 0;

    for (int i = 0; i < C_size[k]; i++)
            if (XXcplex.isExtracted(XXvars[i]))
                    if (XXcplex.getValue(XXvars[i]) > 0.5)
                            stable_set[counter2++] = i;
    int stable_set_size = counter2;
    stable_weight = XXcplex.getObjValue();

    XXcplex.end();
    XXvars.end();
    XXmodel.end();
    XXenv.end();                     
*/
/*****************************************************************************************/