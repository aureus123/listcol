#include "lp.h"
#include "stable.h"
#include "io.h"
#include <cmath>
#include <string>
#include <map>
#include <cfloat>

#define EPSILON 0.00001
#define MAXTIME 7200.0
//#define SAVELP "form.lp"
//#define PUREBYB
#define NOMEMEMPHASIS
#define FICTIONAL_COST 1000

Lopt::Lopt(Graph& G, vector<int>& cost_list) : Xmodel(Xenv), Xvars(Xenv), Xrestr(Xenv), solution (Xenv), G(G), cost_list(cost_list) {}

Lopt::~Lopt() {
    Xrestr.end();
    Xvars.end();
    Xmodel.end();
    Xenv.end();
}

int Lopt::optimize (double& obj_value) {

    // COLUMN GENERATION //

	Xobj = IloMinimize(Xenv);

	// We will have "vertices" constraints with r.h.s >= 1 and "colors" constraints with r.h.s >= -1
	for (int v = 0; v < G.vertices; v++) 
        Xrestr.add(IloRange(Xenv, 1.0, IloInfinity));
	for (int k = 0; k < G.colors; k++) 
        Xrestr.add(IloRange(Xenv, -1.0, IloInfinity));

    // Initial LP solution
	initialize_LP(G);

	Xmodel.add(Xobj);
	Xmodel.add(Xrestr);

	IloCplex cplex(Xmodel);

	// Settings
    set_params(cplex);

    // Initialize MWSS solver
    MWSS<Sewell> solver (G);

    // Generate columns
	while (true) {

        // Solve LP
        cplex.solve();

        // Now, find an entering column (if exists)

        // Save dual values
        IloNumArray duals (Xenv,G.vertices+G.colors);
        cplex.getDuals(duals, Xrestr);

        int added_columns = 0;

        for (int k = 0; k < G.colors;k++) {

                // Get Vk
                vector<int> Vk;
                G.get_Vk(k,Vk);

                // Get duals values of Vk
                vector<double> pi (Vk.size());
                for (unsigned int i = 0; i < Vk.size(); i++) {
                        double d = duals[Vk[i]];
                        if (d > -EPSILON && d < 0.0) d = 0.0;
                        pi[i] = d;
                }

                vector<int> stable_set;
                double stable_weight = 0.0; 
                double goal = cost_list[k] + duals[G.vertices + k];
                solver.solve(k, pi, goal, stable_set, stable_weight);

                // Add column if the reduced cost is negative
                if (goal - stable_weight < -EPSILON) {

                        IloNumColumn column = Xobj(cost_list[k]);
                        // fill the column corresponding to ">= 1" constraints
                        for (unsigned int i = 0; i < stable_set.size(); i++) column += Xrestr[Vk[stable_set[i]]](1.0);
                        // and the ">= -1 constraint
                        column += Xrestr[G.vertices + k](-1.0);

                        /* add the column as a non-negative continuos variable */
                        Xvars.add(IloNumVar(column));
                        ++added_columns;
                }
        }

        //cout << added_columns << " columns were added"<< endl;

        if (added_columns == 0)
                break; // optimality reached
    }


    // Cut fictional columns
    IloExpr restr(Xenv);
    for (int n = 0; n < G.vertices; ++n) {
        restr += Xvars[n];
    }
	Xmodel.add(restr <= 0);
    restr.end();
    cplex.solve();

    // Recover LP solution
    cplex.getValues(solution, Xvars);

    /* LP treatment */
    IloCplex::CplexStatus status = cplex.getCplexStatus();

    if (status == IloCplex::Optimal) {

        obj_value = cplex.getObjValue(); 

        // Check for integrality
        bool integral = true;
        for (int i = 0; i < solution.getSize(); i++)
            if (solution[i] > EPSILON && solution[i] < 1 - EPSILON) {
                integral = false;
                break;
            }

        cplex.end();

        if (!integral)
            return 0;
        else
            return 1;
    }

    cplex.end();
    if (status == IloCplex::Infeasible)
        return -1;

    bye("Error solving LP relaxation");
    return 0;

}

void Lopt::find_branching_vertices (int& i, int& j) {

/*
    // ****** Random approach ******

    for (int u = 0; u < G.vertices - 1; ++u)
        for (int v = u+1; v < G.vertices; ++v) {
            
            // u and v must be non-adjacent
            if (G.is_edge(u,v))
                continue;
            else {
                i = u;
                j = v;
                return;
            }
                
        }
*/
  
    // ****** Most fractional approach ******

    // Mapeo de variables a indices
    map<IloNumVarI*,int> map;
    for (int n = 0; n < solution.getSize(); ++n)    
        map[Xvars[n].getImpl()] = n;

    // LP to matrix
    vector<vector<int> > matrix (G.vertices, vector<int> (solution.getSize(),0));
    for(int m = 0; m < G.vertices; ++m)
        for (IloExpr::LinearIterator it = Xrestr[m].getLinearIterator(); it.ok(); ++it)
            if (it.getCoef() == 1)
                matrix[m][map[it.getVar().getImpl()]] = 1;

    double best_rank = DBL_MAX;

    for (int u = 0; u < G.vertices - 1; ++u)
        for (int v = u+1; v < G.vertices; ++v) {
            
            // u and v must be non-adjacent
            if (G.is_edge(u,v))
                continue;

            // L(u) cap L(v) must be non-empty
            if (!G.have_common_color(u,v))
                continue;

            double best_s1 = DBL_MAX; // Stable that contains both u and v
            double best_s2 = DBL_MAX; // Stable that contains only one of {u,v}
            for (int n = 0; n < solution.getSize(); ++n) {
                if ((matrix[u][n] == 1) && (matrix[v][n]) == 1) {
                    if (abs(solution[n] - 0.5) < best_s1)
                        best_s1 = abs(solution[n] - 0.5);
                }
                else if (matrix[u][n] != matrix[v][n]) {
                    if (abs(solution[n] - 0.5) < best_s2)
                        best_s2 = abs(solution[n] - 0.5);
                }
            }
            if (best_s1 + best_s2 < best_rank) {
                best_rank = best_s1 + best_s2;
                i = u;
                j = v;
            }
        }

    if ((i == j) || (G.is_edge(i,j))) bye("Branching error");

    return;

}

void Lopt::save_coloring(vector<int>& f) {

    // Mapeo de variables a indices
    map<IloNumVarI*,int> map;
    for (int n = 0; n < solution.getSize(); ++n)    
        map[Xvars[n].getImpl()] = n;

    // Save optimal integer solution
    vector<int> temp_coloring (G.vertices);
    for (int v = 0; v < G.vertices; v++) {
	    int color_chosen = -1;
	    IloRange expr1 = Xrestr[v];
	    for (IloExpr::LinearIterator it1 = expr1.getLinearIterator(); it1.ok(); ++it1) {
		    IloNumVar var1 = it1.getVar();
		    IloNum coef1 = it1.getCoef();
		    if (solution[map[var1.getImpl()]] > 0.5 && coef1 > 0.5) {
			    for (int k = 0; k < G.colors; k++) {
				    IloRange expr2 = Xrestr[G.vertices + k];
				    for (IloExpr::LinearIterator it2 = expr2.getLinearIterator(); it2.ok(); ++it2) {
					    IloNumVar var2 = it2.getVar();
					    IloNum coef2 = it2.getCoef();
					    if (var2.getId() == var1.getId() && coef2 < -0.5) {
						    color_chosen = k;
						    break;
					    }
				    }
			    }
			    break;
		    }
	    }
	    if (color_chosen == -1) bye("No color chosen!");
	    temp_coloring[v] = color_chosen;
    }
    // Rearrange temp_coloring in terms of the original vertices
    vector<int> new_vertex;
    G.get_new_vertex(new_vertex);
    f.resize(new_vertex.size());
    for (unsigned int i = 0; i < new_vertex.size(); ++i)
        f[i] = temp_coloring[new_vertex[i]];

}

void Lopt::initialize_LP(Graph& G) {

    // Add a fictional column (e_i,-e_i) = (0..010..0,0..0-10..0) for every i = 1,...,vertices
    // This column has got a really high cost
    for (int v = 0; v < G.vertices; v++) {

		IloNumColumn column = Xobj(FICTIONAL_COST);
		// fill the column corresponding to ">= 1" constraint
		column += Xrestr[v](1.0);

		/* add the column as a non-negative continuos variable */
		Xvars.add(IloNumVar(column));

    }

    return;

}

void Lopt::set_params(IloCplex& cplex) {

	cplex.setDefaults();

#ifndef SHOWCPLEX
	cplex.setOut(Xenv.getNullStream());
	cplex.setWarning(Xenv.getNullStream());
#endif
	//cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Algorithm::Primal);
	cplex.setParam(IloCplex::IntParam::MIPDisplay, 3);
	cplex.setParam(IloCplex::NumParam::WorkMem, 2048);
	cplex.setParam(IloCplex::NumParam::TreLim, 2048);
	cplex.setParam(IloCplex::IntParam::NodeFileInd, 0);
	cplex.setParam(IloCplex::NumParam::TiLim, MAXTIME);
	cplex.setParam(IloCplex::NumParam::EpGap, 0.0);
	cplex.setParam(IloCplex::NumParam::EpAGap, 0.0);
	cplex.setParam(IloCplex::NumParam::EpInt, EPSILON);
	cplex.setParam(IloCplex::IntParam::Threads, 1);
	cplex.setParam(IloCplex::IntParam::RandomSeed, 1);
#ifdef NOMEMEMPHASIS
	cplex.setParam(IloCplex::BoolParam::MemoryEmphasis, CPX_OFF);
#else
	cplex.setParam(IloCplex::BoolParam::MemoryEmphasis, CPX_ON);
#endif

#ifdef PUREBYB
	/* set Traditional B&C with pseudo costs branching strategy */
	cplex.setParam(IloCplex::MIPSearch, 1);
	cplex.setParam(IloCplex::VarSel, CPX_VARSEL_PSEUDO);

	/* turn off other features, including presolve */
	cplex.setParam(IloCplex::PreInd, CPX_OFF);
	cplex.setParam(IloCplex::LBHeur, 0);
	cplex.setParam(IloCplex::Probe, -1);
	cplex.setParam(IloCplex::HeurFreq, -1);
	cplex.setParam(IloCplex::RINSHeur, -1);
	cplex.setParam(IloCplex::RepeatPresolve, 0);

	/* turn off cuts */
	cplex.setParam(IloCplex::Cliques, -1);
	cplex.setParam(IloCplex::Covers, -1);
	cplex.setParam(IloCplex::DisjCuts, -1);
	cplex.setParam(IloCplex::FlowCovers, -1);
	cplex.setParam(IloCplex::FlowPaths, -1);
	cplex.setParam(IloCplex::FracCuts, -1);
	cplex.setParam(IloCplex::GUBCovers, -1);
	cplex.setParam(IloCplex::ImplBd, -1);
	cplex.setParam(IloCplex::MIRCuts, -1);
	cplex.setParam(IloCplex::ZeroHalfCuts, -1);
	cplex.setParam(IloCplex::MCFCuts, -1);
	cplex.setParam(IloCplex::LiftProjCuts, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::LocalImplied, -1);
#endif

	cplex.extract(Xmodel);
#ifdef SAVELP
	cplex.exportModel(SAVELP);
	cout << "Integer formulation saved" << endl;
#endif

    return;

}
