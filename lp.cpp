#include "lp.h"
#include "time.h"
#include "io.h"
#include <cmath>

#include <algorithm>


LP::LP(Graph *G, LP *father, double start_t) : 

    start_t(start_t),
    Xmodel(Xenv), 
    Xvars(Xenv), 
    Xrestr(Xenv),
    vars(Xenv),
    G(G) {

    initialize(father);

}

LP::~LP() {
    values.end();
    for (int i = 0; i < vars.getSize(); ++i)
        delete[] vars[i].stable;
    vars.end();
    Xrestr.end();
    Xvars.end();
    Xobj.end();
    Xmodel.end();
    Xenv.end();
    delete G;
}

void LP::reset() {

    // Clear LP
    values.end();
    for (int i = 0; i < vars.getSize(); ++i)
        delete[] vars[i].stable;
    vars.end();
    Xrestr.end();
    Xvars.end();
    Xobj.end();
    Xmodel.end();
    Xenv.end();

    // Reset LP
    Xenv = IloEnv();
    Xmodel = IloModel(Xenv);
    Xvars = IloNumVarArray(Xenv);
    Xrestr = IloRangeArray(Xenv);
    vars = IloArray<Column>(Xenv);
    pos_vars.clear();

    // Initialize LP
    initialize(this);

}

void LP::initialize(LP *father) {

	if (father == NULL) {

#ifdef ONLY_RELAXATION
		std::cout << "Only relaxation is solved\n";
#endif

        // Preprocessor
        if (PREPROCESSING == 0)
            std::cout << "Preprocessing is turned off" << std::endl;
        else if (PREPROCESSING == 1)
            std::cout << "Preprocessing is enabled - Level 1 (the graph is not modified)" << std::endl;
        else if (PREPROCESSING == 2)
            std::cout << "Preprocessing is enabled - Level 2 (the graph is modified)" << std::endl;
        else
            bye("Undefined preprocessing level");

        // LP initialization strategy:
        if (INITIAL_COLUMN_STRATEGY == 0)
            std::cout << "Dummy strategy selected\n";
        else if (INITIAL_COLUMN_STRATEGY == 1) {
            if (PREPROCESSING > 0)
                bye("CNN strategy is not implemented yet when preprocessing is enabled");    
            std::cout << "CCN strategy selected\n";
        }
        else if (INITIAL_COLUMN_STRATEGY == 2)
            std::cout << "PSC strategy selected\n";
        else
            bye("Undefined initial column strategy");

#ifndef ONLY_RELAXATION

        // Branching rule information
        if (BRANCHING_STRATEGY == 0)
            std::cout << "Branching on edges (EDG strategy)\n";
        else if (BRANCHING_STRATEGY == 1) {
            if (INITIAL_COLUMN_STRATEGY == 1)
                bye("Branching on colors is not implemented yet when CNN strategy is selected");
            if (N_BRANCHS < 2)
                bye("Undefined number of branchs");
            std::cout << "Branching on colors (CLR strategy)\n";
            std::cout << "Number of branchs: " << N_BRANCHS << std::endl;        
        }
        else 
            bye("Undefined branching strategy :(");

        if (VARIABLE_SELECTION == 0 && BRANCHING_STRATEGY == 0)
            std::cout << "Variable selection: Sewell's strategy\n";
        else if (VARIABLE_SELECTION == 1 && BRANCHING_STRATEGY == 0)
            std::cout << "Variable selection: Smallest list of colors\n";
        else if (VARIABLE_SELECTION == 0 && BRANCHING_STRATEGY == 1)
            std::cout << "Variable selection: Greatest neighborhood\n";
        else if (VARIABLE_SELECTION == 1 && BRANCHING_STRATEGY == 1)
            std::cout << "Variable selection: Smallest list of colors\n";
        else
            bye("Undefined variable selection strategy");

#ifdef STABLE_POOL
        if (BRANCHING_STRATEGY == 1)
            bye("Pool of columns is not implemented yet when branching on colors is selected");
        if (PREPROCESSING > 0)
            bye("Pool of columns is not implemented yet when preprocessing is enabled");
        std::cout << "Pool of columns enables\n";
#endif

#endif

	}

    G->preprocess_instance();

    // Initialize vertex and color constraints
    // We will have "vertices" constraints with r.h.s >= 1 and "colors" constraints with r.h.s >= -1
    for (int v = 0; v < G->get_n_vertices(); v++) 
        Xrestr.add(IloRange(Xenv, 1.0, IloInfinity));
    for (int k = 0; k < G->get_n_colors(); k++)
        Xrestr.add(IloRange(Xenv, -1 * std::min(G->get_n_C(k), G->get_n_V(k)), IloInfinity));
    Xmodel.add(Xrestr);

    // Initialize objective function
    Xobj = IloMinimize(Xenv, G->get_precoloring_value());
    Xmodel.add(Xobj);

    fill_initial_columns(father);

    return;
}


void LP::fill_initial_columns (LP *father) {

    if (INITIAL_COLUMN_STRATEGY == 0) {

        // Add a fictional column for each vertex 
        for (int v = 0; v < G->get_n_vertices(); ++v)
            add_fictional_column(v);

    }

    else if (INITIAL_COLUMN_STRATEGY == 1) {

        if (father == NULL) {

            // Root node: use PSC strategy
            std::vector<std::list<nodepntArray>> stable_sets;
            G->coloring_heuristic(stable_sets);

            // Mark the colored vertices
            std::vector<bool> colored (G->get_n_vertices(), false);

            for (int i = 0; i < G->get_n_colors(); ++i) {
                for (auto &s: stable_sets[i]) {
                    // Build the column
                    bool *column = NULL;
                    G->stable_to_column(s.list, s.n_list, &column);

                    // Add the column to the son
                    add_real_column(column, s.n_list, i);

                    // Mark the colored vertices
                    for (int j = 0; j < G->get_n_vertices(); ++j)
                        if (column[j]) { colored[j] = true; }

                    // Free the stable set
                    free(s.list);
                }
            }

                // Add dummy columns for uncovered vertices
            for (int j = 0; j < G->get_n_vertices(); ++j) 
                if (!colored[j]) add_fictional_column(j);

        }

        else {

            if (BRANCHING_STRATEGY == 0) {

                switch (G->get_branch_status())
                {
                case JOIN:
                    // Fictional columns to cover u and v
                    add_fictional_column(G->get_vertex_u());
                    add_fictional_column(G->get_vertex_v());
                    break;    
                case COLLAPSE:
                    // Fictional column to cover u
                    add_fictional_column(G->get_vertex_u());
                    break;
                default:
                    bye("Undefined branching status");
                    break;
                }
            }

            else {
                bye("Undefined fill_initial_columns() for the current branching strategy");
            }

            // Add the active columns from the father
            for (int i: father->pos_vars) {

                // Build the stable set of the father
                nodepnt *stable_father = NULL;
                father->G->column_to_stable(father->vars[i].stable, father->vars[i].n_stable, &stable_father);

                // Build the stable set of the son
                nodepnt *stable_son = NULL;
                int n_stable_son = 0;
                G->translate_stable_set(father->vars[i].color, stable_father, father->vars[i].n_stable, &stable_son, &n_stable_son);

                // Build the column of the son
                bool *column_son = NULL;
                G->stable_to_column(stable_son, n_stable_son, &column_son);

                // Add the column to the son
                add_real_column(column_son, n_stable_son, father->vars[i].color);

                free(stable_father);
                free(stable_son);

            }
        }

    }

    else if (INITIAL_COLUMN_STRATEGY == 2) {

        // Apply the stable set covering heuristic
        std::vector<std::list<nodepntArray>> stable_sets;
        G->coloring_heuristic(stable_sets);

        // Mark the colored vertices
        std::vector<bool> colored (G->get_n_vertices(), false);

        for (int i = 0; i < G->get_n_colors(); ++i) {
            for (auto &s: stable_sets[i]) {
                // Build the column
                bool *column = NULL;
                G->stable_to_column(s.list, s.n_list, &column);

                // Add the column to the son
                add_real_column(column, s.n_list, i);

                // Mark the colored vertices
                for (int j = 0; j < G->get_n_vertices(); ++j)
                    if (column[j]) { colored[j] = true; }

                // Free the stable set
                free(s.list);
            }
        }


        // Add dummy columns for uncovered vertices
        for (int j = 0; j < G->get_n_vertices(); ++j)
            if (!colored[j])
                add_fictional_column(j);

    }

}


void LP::add_real_column (bool *stable, int n_stable, int color) {

    // Add column to CPLEX
    IloNumColumn column = Xobj(G->get_color_cost(color));
    // Fill the column corresponding to ">= 1" constraints
    for (int j = 0; j < G->get_n_vertices(); ++j)
        if (stable[j])
            column += Xrestr[j](1.0);
    // and the ">= -1 constraint
    column += Xrestr[G->get_n_vertices() + color](-1.0);
    // add the column as a non-negative continuos variable
    Xvars.add(IloNumVar(column));

    // Add the column to our internal representation
    vars.add(Column(stable, n_stable, color, false)); 

    return;

}

void LP::add_fictional_column (int v) {

    // Add column to CPLEX
    IloNumColumn column = Xobj(FICTIONAL_COST);
    // Fill the column corresponding to ">= 1" constraints
    column += Xrestr[v](1.0);
    // add the column as a non-negative continuos variable
    Xvars.add(IloNumVar(column));

    // Add the column to our internal representation
    bool *myColumn = new bool[G->get_n_vertices()];
    for (int i = 0; i < G->get_n_vertices(); ++i)
        myColumn[i] = false;
    myColumn[v] = true;
    vars.add(Column(myColumn, 1, -1, true)); 

}


LP_STATE LP::optimize() {

    IloCplex cplex(Xmodel);

    // Set CPLEX's parameters

    cplex.setDefaults();

#ifndef SHOWCPLEX
    cplex.setOut(Xenv.getNullStream());
    cplex.setWarning(Xenv.getNullStream());
#endif
#ifdef VISUALC
    cplex.setParam(IloCplex::IntParam::ClockType, 2); /* set wall-clock time */
#else
    cplex.setParam(IloCplex::IntParam::ClockType, 1); /* set user time */
#endif
    cplex.setParam(IloCplex::NumParam::WorkMem, 2048);
    cplex.setParam(IloCplex::IntParam::Threads, 1);
    cplex.setParam(IloCplex::IntParam::RandomSeed, 1);
    //cplex.setParam(IloCplex::IntParam::RootAlg, IloCplex::Algorithm::Barrier); /* Barrier optimizer */
    cplex.extract(Xmodel);
#ifdef SAVELP
    cplex.exportModel(SAVELP);
    show("Linear relaxation formulation saved");
#endif

    // Solve
    double first_t = start_t;
    int akku_added_columns = 0;

    while (true) {

        // Set time limit
        double time_limit = MAXTIME - (ECOclock() - start_t);
        if (time_limit < 0) {
            return TIME_OR_MEM_LIMIT;
        }
        cplex.setParam(IloCplex::NumParam::TiLim, time_limit);

        // Solve LP
        cplex.solve();

        // Handle non-optimality
        IloCplex::CplexStatus status = cplex.getCplexStatus();
        if (status != IloCplex::Optimal) {
            cplex.end();
            if (status == IloCplex::Infeasible) {bye("ERROR: LPext is infeasible");}
            else if (status == IloCplex::AbortTimeLim 
                  || status == IloCplex::MemLimFeas 
                  || status == IloCplex::MemLimInfeas) {return TIME_OR_MEM_LIMIT;}
            else {bye("ERROR: Unknown LPext status");}
        }

        // Save dual values
#ifdef ONLY_RELAXATION
	double now_t = ECOclock();
	if (now_t - first_t >= 10.0) {
		first_t = now_t;
		std::cout << "  Obj value = " << cplex.getObjValue() << ", Akku added columns = " << akku_added_columns << ", time elapsed = " << now_t - start_t << "s.\n";
	}
#endif
        IloNumArray dual_values (Xenv, G->get_n_vertices() + G->get_n_colors());
        cplex.getDuals(dual_values, Xrestr);
        for (int i = 0; i < G->get_n_vertices(); ++i)
            G->set_vertex_weight(i, dual_values[i]);

        // Now, find an entering column (if exists)
        int added_columns = 0;
        for (int i = 0; i < G->get_n_colors(); i++) {

            // Empty Gk
            if (G->get_n_V(i) == 0) continue;

            // Solve MWSSP
            double goal = G->get_color_cost(i) + dual_values[G->get_n_vertices() + i] + EPSILON;
            nodepnt *max_stable = NULL;
            int n_max_stable = 0;

            if (G->solve_MWSSP(i, goal, &max_stable, &n_max_stable)) {
                // Translate the stable set to a column
                bool *column = NULL;
                G->stable_to_column(max_stable, n_max_stable, &column);
                add_real_column(column, n_max_stable, i);
#ifndef STABLE_POOL
                free(max_stable);   // Otherwise, the stable is saved into the pool and must not be free
#endif
                ++added_columns;
            }
        }
	akku_added_columns += added_columns;

        if (added_columns == 0)
            break; // optimality reached
    }

    // Recover primal values and objective value
    values = IloNumArray(Xenv, Xvars.getSize());
    cplex.getValues(values, Xvars);
    obj_value = cplex.getObjValue(); 

    // Check if there are active fictional columns
    for (int i = 0; i < Xvars.getSize(); ++i)
        if (vars[i].fictional && values[i] > EPSILON) {

            // The solution of LPext is not solution of LP. We must decide the feasibility of LP

            // Set the feasibility threshold
            int feas_threshold = 0;
            for (int v = 0; v < G->get_n_vertices(); ++v)
                feas_threshold += G->get_m(v);

            if (obj_value > feas_threshold + EPSILON) {
                // LP is infeasible
                cplex.end();
                return INFEASIBLE;
            }
            else {
                // TODO: 2-phase LP initilization
                cplex.end();
                bye("ERROR: Failure in LPext initialization");
            }

        }

    // Otherwise, the solution of LPext is also solution of LP
    // Check for integrality and save most fractional variable (for the branching procedure)
    // Note: the optimal solution satisfies 0 <= x[S][k] <= 1 for all S, k 
    // (As long as the weights are > 0)

    bool is_integer = true;
    float most_fract_value = 0.5;
    for (int i = 0; i < Xvars.getSize(); ++i) {
        if (values[i] < EPSILON)
            continue;
        if (values[i] < 1 - EPSILON) {
            is_integer = false;
            // Update the most fractional variable
            if (vars[i].n_stable >= 2 && std::abs(values[i] - 0.5) < most_fract_value) {
                most_fract_value = std::abs(values[i] - 0.5);
                most_fract_var = i;
            }
        }
        // Save positive variable
        pos_vars.push_back(i);
    }

    cplex.end();

    if (is_integer) {
        obj_value = round(obj_value);
        return INTEGER;
    }
    else {
        return FRACTIONAL;
    }

}

void LP::save_solution(std::vector<int> &coloring, std::set<int> &active_colors, double &value) {

    value = G->get_precoloring_value();
    std::vector<int> stables_per_color (G->get_n_colors(), 0);
    std::vector<int> temp_coloring (G->get_n_vertices());

    // Build the coloring of the current graph
    for (int i: pos_vars) {
        int color = vars[i].color;
        int true_color = G->get_C(color, stables_per_color[color]);    
        stables_per_color[color]++;
        value += G->get_color_cost(color);
        for (int j = 0; j < G->get_n_vertices(); ++j)
            if (vars[i].stable[j])
                temp_coloring[j] = true_color;
    }

    // Build the coloring of the original graph
    coloring.resize(G->get_n_total_vertices());
    for (int i = 0; i < G->get_n_total_vertices(); ++i) {
        int cv = G->get_current_vertex(i);    
        if (cv == -1)
            coloring[i] = G->get_precoloring(i);
        else
            coloring[i] = temp_coloring[cv];
        active_colors.insert(coloring[i]);
    }
	
    return;
}


double LP::get_obj_value() {
    return obj_value;
}


int LP::get_n_columns() {
    return Xvars.getSize();
}


void LP::branch(std::vector<LP *> &ret) {

#ifdef STABLE_POOL
    // Update the global pool
    G->update_pool();
#endif

    if (BRANCHING_STRATEGY == 0)
        branch_on_edges(ret);
    else if (BRANCHING_STRATEGY == 1)
        branch_on_colors(ret);

    return;

}

void LP::branch_on_edges(std::vector<LP *> &ret) {

    int u = -1, v = -1;

    // ** Variable selection ** 

    if (VARIABLE_SELECTION == 0) {

        // Choose vertices u and v to branch
        // 1) Find the most fractional stable (S1,k) with |S1| > 1
        // 2) Find the first vertex u in S1
        // 3) Find the first active stable set (S2,l) such that u \in S2 [and (S1,k) \neq (S2,l)]
        // 4) If S2 \neq S1, choose v \in S1 (symmetric difference) S2
        // 5) Otherwhise, choose the second vertex v \in S1

        // Find the first vertex u in S1
        if (most_fract_var < 0 || most_fract_var >= Xvars.getSize()) bye("Branching error");
        for (int i = 0; i < G->get_n_vertices(); ++i)
            if (vars[most_fract_var].stable[i]) {
                u = i;
                break;
            }
        if (u < 0) bye("Branching error");

        // Find the first active stable set (S2,l) such that u \in S2 [and (S1,k) \neq (S2,l)]
        for (int i: pos_vars)
            if (vars[i].stable[u] && most_fract_var != i) {
                for (int j = 0; j < G->get_n_vertices(); ++j)
                    // If S2 \neq S1, choose v \in S1 (symmetric difference) S2
                    if (vars[i].stable[j] != vars[most_fract_var].stable[j]) {
                        v = j;
                        break;
                    }
                // Otherwhise, choose the second vertex v \in S1
                if (v < 0) {
                    for (int j = u+1; j < G->get_n_vertices(); ++j)
                        if (vars[most_fract_var].stable[j]) {
                            v = j;
                            break;
                        }
                    if (v < 0) bye("Branching error");
                }
                break;
            }

    }

    else if (VARIABLE_SELECTION == 1) {

        // Choose u and v such that
        //      there is a fractional variable x(S,k) with u,v \in S
        //      u and v maximises k(u) + k(v), where k(u) = |L(u) \¢ap K|
        //      Tie breaking rule: lowest index 

        std::vector<int> W;
        G->get_W1(W);

        int min = INT_MAX;
        std::vector<std::vector<bool>> repeated (G->get_n_vertices(), std::vector<bool> (G->get_n_vertices(), false));
        for (int x: pos_vars)
            if (values[x] > EPSILON && values[x] < 1 - EPSILON)
                for (int i = 0; i < G->get_n_vertices(); ++i)
                    if (vars[x].stable[i]) 
                        for (int j = i+1; j < G->get_n_vertices(); ++j)
                            if (vars[x].stable[j])
                                if (!repeated[i][j]) {
                                    repeated[i][j] = true;
                                    int m = std::min(W[i],W[j]);
                                    if (m < min
                                    || (m == min && i < u)
                                    || (m == min && i == u && j < v)) {
                                        min = m;
                                        u = i;
                                        v = j;
                                    }
                                }

    }

    // ** Branching rule ** 

    // Build the graph with an additional edge between u and v
    Graph *G1 = G->join_vertices(u, v);

    // Build LP from G1
    LP *lp1 = new LP(G1, this, start_t);
  
    // Build the graph by collapsing the vertices u and v
    Graph *G2 = G->collapse_vertices(u, v);

    // Build LP from G2
    LP *lp2 = new LP(G2, this, start_t);

    ret.resize(2);
    ret[0] = lp2;   // First collapse
    ret[1] = lp1;   // Then join

    return;

}


bool is_better (std::tuple<int,int,int> &t1, std::tuple<int,int,int> &t2) {
    int d1, d2, j1, j2, m1, m2;
    std::tie(d1, j1, m1) = t1;
    std::tie(d2, j2, m2) = t2;
    if (d1 > d2)
        return true;
    else if (d1 < d2)
        return false;
    else if (m1 < m2)
        return true;
    else if (m1 > m2)
        return false;
    else if (j1 < j2)
        return true;
    else if (j1 > j2)
        return false;
    return false;
};


void LP::branch_on_colors(std::vector<LP *> &ret) {

    // Choose vertex v and color j to branch with k(v) >= 2

    bool done = false;
    while (!done) {

        int v = -1, j = -1, d = -1, m = -1;

        std::vector<int> W;
        G->get_W1(W); // W[v] = k(v) = |L[v] \cap K|

        if (VARIABLE_SELECTION == 0) {

            // Choose v and j with:
            //      there is some fractional variable (S,j) with v \in S
            //      v maximises |N_{Gj}(v)|
            //      First tie breaking rule: lowest k(v)
            //      Second tie breaking rule: lowest multiplicity
            //      Third tie breaking rule: lowest index of color
            //      Fourth tie breaking rule: lowest index of vertex
            // If k(v) = 1 and m(j) = 1
            //      Report an error
            // Else if k(v) = 1 and m(j) > 1 (this can only happens when the preprocessor is not activated)
            //      Preprocess (v,j) and apply again the variable selection strategy
    
            std::vector<std::vector<bool>> repeated (G->get_n_vertices(), std::vector<bool> (G->get_n_colors(), false));
            int k = INT_MAX;
            for (int i: pos_vars) {
                if (values[i] > EPSILON && values[i] < 1 - EPSILON) { // x(S,color) is fractional
                    int color = vars[i].color;
                    for (int u = 0; u < G->get_n_vertices(); ++u) {
                        if (vars[i].stable[u]) { // u is in S
                            if (!repeated[u][color]) {

                                repeated[u][color] = true;
                                int delta = G->get_n_neighbours(u, color);
                                int mult = G->get_n_C(color);

                                if (W[u] == 1 && delta == 0)
                                    continue;

                                if (delta > d
                                || (delta == d && W[u] < k) // First tie breaking rule
                                || (delta == d && W[u] == k && mult < m) // Second tie breaking rule
                                || (delta == d && W[u] == k && mult == m && color < j) // Third tie breaking rule
                                || (delta == d && W[u] == k && mult == m && color == j && u < v)) { // Fourth tie breaking rule
                                    d = delta;
                                    k = W[u];
                                    j = color;
                                    m = mult;
                                    v = u;
                                }

                            }
                        }
                    }
                }
            }
        }

        else if (VARIABLE_SELECTION == 1) {

            // Choose v such that:
            //      there is some fractional variable x(S,j) with v \in S
            //      minimises k(v)
            //      Tie breaking rule: lowest index
            // Choose j such that
            //      there is some fractional variable x(S,j) with v \in S
            //      maximises |N_{Gj}(v)|
            //      First tie breaking rule: lowest multiplicity
            //      Second tie breaking rule: lowest index of color
            // If k(v) = 1 and m(j) = 1
            //      Report an error
            // Else if k(v) = 1 and m(j) > 1 (this can only happens when the preprocessor is not activated)
            //      Preprocess (v,j) and apply again the variable selection strategy

            std::vector<std::vector<bool>> repeated (G->get_n_vertices(), std::vector<bool> (G->get_n_colors(), false));
            int k = INT_MAX;
            for (int i: pos_vars) {
                if (values[i] > EPSILON && values[i] < 1 - EPSILON) { // x(S,color) is fractional
                    int color = vars[i].color;
                    for (int u = 0; u < G->get_n_vertices(); ++u) {
                        if (vars[i].stable[u]) { // u is in S
                            if (!repeated[u][color]) {

                                repeated[u][color] = true;
                                int delta = G->get_n_neighbours(u, color);
                                int mult = G->get_n_C(color);
                                
                                if (W[u] == 1 && delta == 0) {
                                    continue;
                                }

                                if (W[u] < k
                                || (W[u] == k && u < v) // Vertex's tie breaking rule
                                || (W[u] == k && u == v && delta > d)
                                || (W[u] == k && u == v && delta == d && mult < m) // Color's first tie breaking rule
                                || (W[u] == k && u == v && delta == d && mult == m && color < j)) { // Color's second tie breaking rule
                                    v = u;
                                    j = color;
                                    k = W[u];
                                    d = delta;
                                    m = mult;
                                }

                            }
                        }
                    }
                }
            }

        }


        if (W[v] == 1 && d == 0) { // An isolated vertex was selected
            bye("Branching on color: an isolated vertex was selected");
        }
        else if (W[v] == 1 && d > 0) {

            if (PREPROCESSING > 0)
                bye("Branching on color: internal error");

            G->preprocess_instance(v,j);

            reset();
            optimize();

        }
        else {

            done = true;

            std::list<std::tuple<int,int,int>> candidates;
            candidates.push_back(std::make_tuple(d,j,m)); 

            if (N_BRANCHS > 2) {

                // Find colors k != j such that
                //      v \in S for some x(S,k) fractional
                //      v maximises |N_{G_k}(v)|

                std::vector<bool> repeated (G->get_n_colors(), false);
                repeated[j] = true;
                for (int i: pos_vars) {
                    if (values[i] > EPSILON && values[i] < 1 - EPSILON) { // x(S,k) is fractional
                        int k = vars[i].color;
                        if (!repeated[k] && vars[i].stable[v]) { // v \in S

                            repeated[k] = true;
                            int delta = G->get_n_neighbours(v,k);
                            int mult = G->get_n_C(k);
                            std::tuple<int,int,int> t = std::make_tuple(delta, k, mult);

                            if (candidates.size() == N_BRANCHS - 1 && is_better(candidates.back(), t)) {
                                // If the list is full and the new candidate is worse than the last element, then continue
                                continue;
                            }
                            else {
                                // Traverse the list and find the right position to add the candidate
                                // The list is sorted from best to worst candidates
                                bool added = false;
                                for (auto it = candidates.begin(); it != candidates.end(); ++it)
                                    if (is_better(t,*it)) {
                                        candidates.insert(it,t);
                                        added = true;
                                        break;
                                    }
                                if (!added) {
                                    // Add the element at the back of the list
                                    candidates.push_back(t);
                                }
                                if (candidates.size() == N_BRANCHS) {
                                    // Remove the worst candidate
                                    candidates.pop_back();
                                }
                            }
                        }
                    }
                }
            }

            ret.resize(candidates.size() + 1);
            std::set<int> colors;
            int i = 0;

            for (auto it = candidates.begin(); it != candidates.end(); ++it) {

                int k = std::get<1>(*it);
 
                Graph *G1;
                // Build the graph where v is colored with some color of C[j]
                G1 = G->choose_color(v, k);
                // Build LP from G1
                LP *lp1 = new LP(G1, this, start_t);

                ret[i++] = lp1;
                colors.insert(k);

            }

            Graph *G2;
            // Build the graph where v is not colored with any color of C[k] with k \in colors
            G2 = G->remove_color(v, colors);
            // Build LP from G1
            LP *lp2 = new LP(G2, this, start_t);

            ret[i] = lp2;

        }

    }

    return;

}