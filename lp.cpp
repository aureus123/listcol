#include "graph.h"
#include "MWSS.h"
#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>

#define EPSILON 0.00001

static void initialize_LP();
static void set_params(IloCplex);

bool optimize (Graph& G) {

    // COLUMN GENERATION //

	Xobj = IloMinimize(Xenv);

	// We will have "vertices" constraints with r.h.s >= 1 and "colors" constraints with r.h.s >= -1
	for (int v = 0; v < vertices; v++) 
        Xrestr.add(IloRange(Xenv, 1.0, IloInfinity));
	for (int k = 0; k < colors; k++) 
        Xrestr.add(IloRange(Xenv, -1.0, IloInfinity));

    // Initial LP solution
	initialize_LP();

	Xmodel.add(Xobj);
	Xmodel.add(Xrestr);

	// Settings
	IloCplex cplex(Xmodel);

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
                vector<double> pi;
                for (int i = 0; i < Vk.size(); i++) {
                        double d = duals[Vk[i]];
                        if (d > -EPSILON && d < 0.0) d = 0.0;
                        pi.push_back() = d;
                }

                vector<int> stable_set;
                double stable_weight = 0.0; 
                double goal = G.get_cost(k) + duals[G.vertices + k];
                solver.solve(k, pi, goal, stable_set, &stable_weight);

                // Add column if the reduced cost is negative
                if (goal - stable_weight < -EPSILON) {

                        IloNumColumn column = Xobj(G.get_cost(k));
                        // fill the column corresponding to ">= 1" constraints
                        for (int i = 0; i < stable_set.size(); i++) column += Xrestr[Vk[stable_set[i]]](1.0);
                        // and the ">= -1 constraint
                        column += Xrestr[G.vertices + k](-1.0);

                        /* add the column as a non-negative continuos variable */
                        Xvars.add(IloNumVar(column));
                        ++added_columns;
                }
        }

        //cout << columns_added << " columns were added"<< endl;

        if (added_columns == 0)
                break; // optimality reached
    }

	return true;

}







	cplex.setDefaults();
#ifndef SHOWCPLEX
	cplex.setOut(Xenv.getNullStream());
	cplex.setWarning(Xenv.getNullStream());
#endif
	cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Algorithm::Primal);
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



/*****************************************************************************************/
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
