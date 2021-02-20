#include "lp.h"
#include "time.h"
#include "io.h"
#include <cmath>

LP::LP(Graph *G, LP *father) : 

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


void LP::initialize(LP *father) {

	if (father == NULL) {
#ifdef ONLY_RELAXATION
		std::cout << "Only relaxation is solved\n";
#elif BRANCHING_STRATEGY == 0
		std::cout << "Branching on edges (EDG strategy)\n";
#elif BRANCHING_STRATEGY == 1
		std::cout << "Branching on colors (CLR strategy)\n";
#elif BRANCHING_STRATEGY == 2
		std::cout << "Branching on indistinguishable (IND strategy)\n";
#else
		bye("Undefined branching strategy :(");
#endif
#ifdef STABLE_POOL
		std::cout << "Pool of columns enabled\n";
#endif
	}

#if BRANCHING_STRATEGY == 0
#ifdef PREPROCESSING
    G->preprocess_instance();
#endif
#endif
	
    // Initialize vertex and color constraints
    // We will have "vertices" constraints with r.h.s >= 1 and "colors" constraints with r.h.s >= -1
    for (int v = 0; v < G->get_n_vertices(); v++) 
        Xrestr.add(IloRange(Xenv, 1.0, IloInfinity));
    for (int k = 0; k < G->get_n_colors(); k++)
        Xrestr.add(IloRange(Xenv, -1 * G->get_n_C(k), IloInfinity));
    Xmodel.add(Xrestr);

    // Initialize objective function
    Xobj = IloMinimize(Xenv);
    Xmodel.add(Xobj);

    fill_initial_columns(father);

    return;
}


void LP::fill_initial_columns (LP *father) {

#if INITIAL_COLUMN_STRATEGY == 0

	if (father == NULL) std::cout << "Dummy strategy selected\n";

    // Add a fictional column for each vertex 
    for (int v = 0; v < G->get_n_vertices(); ++v)
        add_fictional_column(v);

#elif INITIAL_COLUMN_STRATEGY == 1

	if (father == NULL) {
		std::cout << "CCN strategy selected\n";

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
		for (int j = 0; j < G->get_n_vertices(); ++j) if (!colored[j]) add_fictional_column(j);
    	}
	else {

#if BRANCHING_STRATEGY == 0
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
#else
        bye("Undefined fill_initial_columns() for the current branching strategy");
#endif

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

#elif INITIAL_COLUMN_STRATEGY == 2

	if (father == NULL) std::cout << "PSC strategy selected\n";

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

#else
    bye("Undefined initial column strategy");
#endif

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


LP_STATE LP::optimize(double start_t) {

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

#if BRANCHING_STRATEGY == 0

 
    value = 0.0;
    std::vector<int> stables_per_color (G->get_n_colors(), 0);
    std::vector<int> temp_coloring (G->get_n_vertices());
 
    // Build the coloring of the current graph
    for (int i: pos_vars) {
        int color = vars[i].color;
        int true_color = G->get_C(color, stables_per_color[color]);
        stables_per_color[color]++;
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
	
    // Compute the cost
    value = G->get_precoloring_value();
    for (int k: active_colors)
        value += G->get_color_cost(k);

#else

    coloring.resize(G->get_n_vertices());
    std::vector<int> stables_per_color (G->get_n_colors(), 0);
    for (int i: pos_vars) {
        int color = vars[i].color;
        int true_color = G->get_C(color, stables_per_color[color]);
        stables_per_color[color]++;
        for (int j = 0; j < G->get_n_vertices(); ++j) {
            if (vars[i].stable[j])
                coloring[j] = true_color;
        }
        active_colors.insert(true_color);
    }
    value = get_obj_value();

#endif

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

#if BRANCHING_STRATEGY == 0
    branch_on_edges(ret);
#elif BRANCHING_STRATEGY == 1
    branch_on_colors(ret);
#elif BRANCHING_STRATEGY == 2
    branch_on_colors(ret);
#endif

    return;

}

void LP::branch_on_edges(std::vector<LP *> &ret) {

    // ** Variable selection ** 

    // Choose vertices u and v to branch
    // 1) Find the most fractional stable (S1,k) with |S1| > 1
    // 2) Find the first vertex u in S1
    // 3) Find the first active stable set (S2,l) such that u \in S2 [and (S1,k) \neq (S2,l)]
    // 4) If S2 \neq S1, choose v \in S1 (symmetric difference) S2
    // 5) Otherwhise, choose the second vertex v \in S1

    int u = -1, v = -1;

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

    // ** Branching rule ** 

    // Build the graph with an additional edge between u and v
    Graph *G1 = G->join_vertices(u, v);

    // Build LP from G1
    LP *lp1 = new LP(G1, this);
  
    // Build the graph by collapsing the vertices u and v
    Graph *G2 = G->collapse_vertices(u, v);

    // Build LP from G2
    LP *lp2 = new LP(G2, this);

    ret.resize(2);
    ret[0] = lp2;   // First collapse
    ret[1] = lp1;   // Then join

    return;

}


void LP::branch_on_colors(std::vector<LP *> &ret) {

    // Choose vertex v and color k to branch
    
    // Define W1 = {v : v is in a single color class}
    // Choose (v,k) such that
    //  v is in W1
    //  v is in some Sk with x(S,k) fractional
    //  v maximises |N(v) \cap Vk| + rand(0,1)  (random tie break) 
    // If such (v,k) exists, then reduce the instance
    // Otheriwse, choose (v,k) such that
    //  v is in some Sk with x(S,k) fractional
    //  v maximises |N(v) \cap Vk| + rand(0,1)  (random tie break) 
    // Such (v,k) always exists, then branch according to CLR or IND

    std::vector<int> W1;
    G->get_W1(W1); // W[v] = |Lv \cap K|

    bool empty_W1 = true;
    double max_value = -1.0;
    int max_v = -1;
    int max_k = -1;
    for (int i: pos_vars)
        if (values[i] > EPSILON && values[i] < 1 - EPSILON) // x(S,k) is fractional
            for (int v = 0; v < G->get_n_vertices(); ++v)
                if (vars[i].stable[v]) { // v is in S

                    double value = G->get_n_neighbours(v, vars[i].color) + (double) rand() / RAND_MAX; 

                    if (W1[v] == 1 && empty_W1) {
                        empty_W1 = false;
                        max_value = value;
                        max_v = v;
                        max_k = vars[i].color;
                    }
                    else if (W1[v] == 1 && !empty_W1 && value > max_value) {
                        max_value = value;
                        max_v = v;
                        max_k = vars[i].color;
                    }
                    else if (W1[v] > 1 && empty_W1 && value > max_value) {
                        max_value = value;
                        max_v = v;
                        max_k = vars[i].color;
                    } 
                    
                }

    if (max_v == -1 || max_k == -1)
        bye("Branching error");
    
    if (!empty_W1) {

        // Build the graph where v is colored with k and reduce the instance 
        Graph *G1 = G->choose_color(max_v, max_k);

        // Build LP from G1
        LP *lp1 = new LP(G1, this);

        ret.resize(1);
        ret[0] = lp1;   // First choose

    }
    else {

        Graph *G1, *G2;

#if BRANCHING_STRATEGY == 1

        // Build the graph where v is colored with k and reduce the instance
        G1 = G->choose_color(max_v, max_k);

#elif BRANCHING_STRATEGY == 2

        // Build the graph where v is colored with k, do not change the partition
        G1 = G->choose_indistinguishable_color(max_v, max_k);

#endif

        // Build LP from G1
        LP *lp1 = new LP(G1, this);

        // Build the graph where v cannot be colored with k
        G2 = G->remove_color(max_v, max_k);

        // Build LP from G2
        LP *lp2 = new LP(G2, this);

        ret.resize(2);
        ret[0] = lp1;   // First choose
        ret[1] = lp2;   // Then remove
        
        return;

    }

}

