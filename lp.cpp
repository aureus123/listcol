#include "lp.h"
#include "stable.h"
#include "io.h"
#include <cmath>
#include <string>
#include <map>
#include <cfloat>

LP::LP(Graph* G) : Xmodel(Xenv), Xvars(Xenv), Xrestr(Xenv), solution (Xenv), fictional(true), G(G) {}

LP::~LP() {
    Xrestr.end();
    Xvars.end();
    Xmodel.end();
    Xenv.end();
    delete G;
}

double LP::get_obj_value() {
    return obj_value;
}

LP_STATE LP::optimize (double goal) {

    // COLUMN GENERATION //

	Xobj = IloMinimize(Xenv);

	// We will have "vertices" constraints with r.h.s >= 1 and "colors" constraints with r.h.s >= -1
	for (int v = 0; v < G->vertices; v++) 
        Xrestr.add(IloRange(Xenv, 1.0, IloInfinity));
	for (int k = 0; k < G->colors; k++) {
#ifdef COLORS_DELETION
        Xrestr.add(IloRange(Xenv, G->get_right_hand_side(k), IloInfinity));
#elif
        Xrestr.add(IloRange(Xenv, -1.0, IloInfinity));
#endif
    }

    // Initial LP solution
	initialize_LP();

	Xmodel.add(Xobj);
	Xmodel.add(Xrestr);

	IloCplex cplex(Xmodel);

	// Settings
    set_params(cplex);

    // Initialize MWSS solver
    MWSS<Sewell> solver (*G);

    int total_added_columns = 0;

    // Generate columns
	while (true) {

        // Solve LP
        cplex.solve();

        // Early branching
        if (cplex.getObjValue() < goal - EPSILON)
            break;

        // Now, find an entering column (if exists)

        // Save dual values
        IloNumArray duals (Xenv,G->vertices+G->colors);
        cplex.getDuals(duals, Xrestr);

        int added_columns = 0;

        for (int k = 0; k < G->colors;k++) {

            if (G->get_Vk_size(k) == 0) continue;

            // Get Vk
            vector<int> Vk;
            G->get_Vk(k,Vk);

            // Get duals values of Vk
            vector<double> pi (Vk.size());
            for (unsigned int i = 0; i < Vk.size(); i++) {
                    double d = duals[Vk[i]];
                    if (d > -EPSILON && d < 0.0) d = 0.0;
                    pi[i] = d;
            }

            vector<int> stable_set;
            double stable_weight = 0.0; 
            double goal = G->get_cost(k) + duals[G->vertices + k] + THRESHOLD;
            solver.solve(k, pi, goal, stable_set, stable_weight);

            // Add column if the reduced cost is negative
            if (goal - stable_weight < -EPSILON) {

                    IloNumColumn column = Xobj(G->get_cost(k));
                    // fill the column corresponding to ">= 1" constraints
                    for (unsigned int i = 0; i < stable_set.size(); i++) column += Xrestr[Vk[stable_set[i]]](1.0);
                    // and the ">= -1 constraint
                    column += Xrestr[G->vertices + k](-1.0);

                    /* add the column as a non-negative continuos variable */
                    Xvars.add(IloNumVar(column));
                    ++added_columns;
                    ++total_added_columns;

            }
        }

        //cout << added_columns << " columns were added"<< endl;

        if (added_columns == 0)
                break; // optimality reached
    }

    //cout << total_added_columns << " columns were added"<< endl;

    // Cut fictional columns
    if (fictional) {
        IloExpr restr(Xenv);
        for (int n = 0; n < G->vertices; ++n) {
            restr += Xvars[n];
        }
	    Xmodel.add(restr <= 0);
        restr.end();
        cplex.solve();
    }

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
            return FRACTIONAL;
        else
            return INTEGER;
    }

    cplex.end();
    if (status == IloCplex::Infeasible)
        return INFEASIBLE;

    bye("Error solving LP relaxation");
    return INFEASIBLE;

}

void LP::initialize_LP() {

#ifdef INITIAL_HEURISTIC
    vector<vector<int>> stables_set;
    if (G->coloring_heuristic(stables_set)) {

        // Add columns
        for (auto s: stables_set) {

            int k = s.back(); // color
            s.pop_back();

		    IloNumColumn column = Xobj(G->get_cost(k));
            for (int v: s)
		        column += Xrestr[v](1.0); // fill the column corresponding to ">= 1" constraint
            column += Xrestr[G->vertices+k](-1.0); // fill the column corresponding to ">= -1" constraint
		    Xvars.add(IloNumVar(column)); // add the column as a non-negative continuos variable

        }

        fictional = false;
        return;
    }
#endif

    // Add a fictional columns
    for (int v = 0; v < G->vertices; v++) {
		IloNumColumn column = Xobj(FICTIONAL_COST);
		column += Xrestr[v](1.0); // fill the column corresponding to ">= 1" constraint
	    Xvars.add(IloNumVar(column)); // add the column as a non-negative continuos variable
    }

    return;

}

void LP::branch1 (vector<LP*>& lps) {

    // Find vertices u and v to branch
    int u, v;
    select_vertices(u,v);

    // Create LPs
    lps.reserve(2);
    Graph* G2 = new Graph(*G);      // Copy G
    G2->join_vertices(u,v);
#ifdef COLORS_DELETION
    G2->delete_equal_colors();
#endif
    lps.push_back(new LP(G2));
    G->collapse_vertices(u,v);      // Reuse G (delete_equal_colors is not needed!)
    lps.push_back(new LP(G));
    G = NULL;                       // This prevent G being deleted by the father 

    return;

}

void LP::branch2 (vector<LP*>& lps) {

    // Find vertex to branch
    int v;
    set<int> colors;
    select_vertex(v, colors);

    // Create LPs
    vector<int> Lv;
    G->get_Lv(v, Lv);

    // Unused colors
    vector<int> uc;
    for (int k: Lv)
        if (colors.find(k) == colors.end())
            uc.push_back(k);

    lps.reserve(colors.size() + uc.size());

    if (uc.size() > 0) {
        for (int k : colors) {
            Graph* G2 = new Graph(*G);  // copy graph
            G2->color_vertex(v,k);
            lps.push_back(new LP(G2));
        }
        G->set_Lv(v,uc);                // Reuse G
        lps.push_back(new LP(G));
        G = NULL;                       // This prevent G being deleted by the father 
    }
    else {
        for (int k : colors) {
            Graph* G2 = new Graph(*G);  
            G2->color_vertex(v,k);
            lps.push_back(new LP(G2));
        }
    }

    return;

}

void LP::select_vertices (int& i, int& j) {

/*
    // ****** Random approach ******

    for (int u = 0; u < G->vertices - 1; ++u)
        for (int v = u+1; v < G->vertices; ++v) {
            
            // u and v must be non-adjacent
            if (G->is_edge(u,v))
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
    vector<vector<int> > matrix (G->vertices, vector<int> (solution.getSize(),0));
    for(int m = 0; m < G->vertices; ++m)
        for (IloExpr::LinearIterator it = Xrestr[m].getLinearIterator(); it.ok(); ++it)
            if (it.getCoef() == 1)
                matrix[m][map[it.getVar().getImpl()]] = 1;

    double best_rank = DBL_MAX;

    for (int u = 0; u < G->vertices - 1; ++u)
        for (int v = u+1; v < G->vertices; ++v) {
            
            // u and v must be non-adjacent
            if (G->is_edge(u,v))
                continue;

            // L(u) cap L(v) must be non-empty
            if (!G->have_common_color(u,v))
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

    if ((i == j) || (G->is_edge(i,j))) bye("Branching error");

    return;

}

void LP::select_vertex (int& v, set<int>& colors) {

    // Mapeo de variables a indices
    map<IloNumVarI*,int> map;
    for (int n = 0; n < solution.getSize(); ++n)    
        map[Xvars[n].getImpl()] = n;

    // C(v): set of colors associated to some fractional covering that covers v
    // Find v that minimize |C(v)| > 1
    for (int i = 0; i < G->vertices; i++) {
        set<int> C;
        for (IloExpr::LinearIterator it = Xrestr[i].getLinearIterator(); it.ok(); ++it) {
            if (it.getCoef() < 1) continue;
            double val = solution[map[it.getVar().getImpl()]];
            if ((val > EPSILON) && (val < 1 - EPSILON))
                // Find color
                for (int k = 0; k < G->colors; ++k)
		            for (IloExpr::LinearIterator it2 =  Xrestr[G->vertices + k].getLinearIterator(); it2.ok(); ++it2) {
			            if (it.getVar().getImpl() == it2.getVar().getImpl() && it2.getCoef() < -0.5) {
				            if (C.find(k) == C.end())
                                C.insert(k);
                            break;
                        }
                    }

            if (C.size() < 2) continue;

            if (colors.size() == 0) {
                colors = C;
                v = i;
            }
            else if (colors.size() > 0 && C.size() < colors.size()) {
                colors = C;
                v = i;
            }
        }        
    }

    if (colors.size() < 2) bye("Branching error");

    return;

}

void LP::save_solution(vector<int>& f) {

    // Mapping from variable to index in solution array
    map<IloNumVarI*,int> var_index;
    for (int n = 0; n < solution.getSize(); ++n)    
        var_index[Xvars[n].getImpl()] = n;

    map<IloNumVarI*,int> sol_color; // Mapping from used variables to associated color
#ifdef COLORS_DELETION
    vector<int> stables_per_color (G->colors,0);
#endif
    for (int n = 0; n < solution.getSize(); ++n)
        if (solution[var_index[Xvars[n].getImpl()]] > 0.5) {
            // Find associated color
	        for (int k = 0; k < G->colors; k++) {
		        IloRange expr = Xrestr[G->vertices + k];
		        for (IloExpr::LinearIterator it = expr.getLinearIterator(); it.ok(); ++it) {
			        IloNumVar var = it.getVar();
			        IloNum coef = it.getCoef();
			        if (var.getImpl() == Xvars[n].getImpl() && coef < -0.5) {
#ifdef COLORS_DELETION
                        if (abs(G->get_right_hand_side(k)) > 1) {
				            sol_color[Xvars[n].getImpl()] = G->get_eq_colors(k,stables_per_color[k]);
                            stables_per_color[k]++;
                        }
                        else
                            sol_color[Xvars[n].getImpl()] = k;
#elif
                        sol_color[Xvars[n].getImpl()] = k;
#endif
				        break;
			        }
		        }
	        }
        }
   
    // Save optimal integer solution
    vector<int> temp_coloring (G->vertices);
    for (int v = 0; v < G->vertices; v++) {
	    IloRange expr1 = Xrestr[v];
	    for (IloExpr::LinearIterator it1 = expr1.getLinearIterator(); it1.ok(); ++it1) {
		    IloNumVar var1 = it1.getVar();
		    IloNum coef1 = it1.getCoef();
		    if (solution[var_index[var1.getImpl()]] > 0.5 && coef1 > 0.5)
                temp_coloring[v] = sol_color[var1.getImpl()];
	    }
    }

    // Rearrange temp_coloring in terms of the vertices of the original graph
    vector<int> new_vertex;
    G->get_new_vertex(new_vertex);
    f.resize(new_vertex.size());
    for (unsigned int i = 0; i < new_vertex.size(); ++i)
        f[i] = temp_coloring[new_vertex[i]];

    return;
}

bool LP::check_solution(vector<int>& sol) {
    return G->check_coloring(sol);
}

void LP::set_params(IloCplex& cplex) {

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
