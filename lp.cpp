#include "lp.h"
#include "stable.h"
#include "io.h"
#include <cmath>
#include <string>

#define EPSILON 0.00001
#define MAXTIME 7200.0
//#define SAVELP "form.lp"
//#define PUREBYB
#define NOMEMEMPHASIS
#define FICTIONAL_COST 10000

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

    // Save solution
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
        if (integral)
            return 1;
        else
            return 0;

    }
    cplex.end();
    if (status == IloCplex::Infeasible)
        return -1;

    bye("Error solving LP relaxation");
    return 0;

}

void Lopt::find_branching_vertices (int& i, int& j) {
    
    // Find most fractional column
    double most_frac = 1;
    int s1 = -1;
    for (int k = 0; k < solution.getSize(); k++)    
        if (abs(solution[k]-0.5) < most_frac) {
            most_frac = abs(solution[k]-0.5);
            s1 = k;
        }

    if (s1 == -1) bye("Error finding vertices for branching");

    // Find first row covered by s1
    for (i = 0; i < G.vertices; i++) {
        bool flag = false; // for double breaking
        for (IloExpr::LinearIterator it = Xrestr[i].getLinearIterator(); it.ok(); ++it)
            if ((Xvars[s1].getImpl() == it.getVar().getImpl()) && (it.getCoef() == 1)) { // it's horrible but it's the best i could do
                flag = true;
                break;  
            }
        if (flag) break;
    }
    if (i == G.vertices) bye("Error finding vertices for branching");

    // Find another column that cover i
    IloNumVar s2;
    IloExpr::LinearIterator it;
    for (it = Xrestr[i].getLinearIterator(); it.ok(); ++it) {
        if (it.getVar().getImpl() == Xvars[s1].getImpl()) continue;
        if (it.getCoef() == 1) {
            s2 = it.getVar().getImpl();
            break; // What a mess!
        }
    }
    if (!it.ok()) bye("Error finding vertices for branching");

    // Find row j such that only one of s1 or s2 cover j
    for (j = 0; j < G.vertices; ++j) {
        if (j == i) continue;
        bool cover_s1 = false;
        bool cover_s2 = false;
        for (IloExpr::LinearIterator it = Xrestr[j].getLinearIterator(); it.ok(); ++it) {
            if (it.getVar().getImpl() == Xvars[s1].getImpl()) cover_s1 = true;
            else if (it.getVar().getImpl() == s2.getImpl()) cover_s2 = true;
        }
        if (cover_s1 != cover_s2) break;
    }
    if (j == G.vertices) bye("Error finding vertices for branching");   

    return;
}

void Lopt::initialize_LP(Graph& G) {

    if (G.colors < G.vertices) {
        cout << "Error initializing LP relaxation: number of vertices is grater than the number of colors" << endl;
        return;
    }

    // Add a fictional column (e_i,-e_i) = (0..010..0,0..0-10..0) for every i = 1,...,vertices
    // This column has got a really high cost
    for (int v = 0; v < G.vertices; v++) {

		IloNumColumn column = Xobj(FICTIONAL_COST);
		// fill the column corresponding to ">= 1" constraint
		column += Xrestr[v](1.0);
		// and the ">= -1 constraint
		column += Xrestr[G.vertices + v](-1.0);

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
