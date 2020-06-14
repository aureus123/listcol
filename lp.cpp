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

    // Add a fictional column for each vertex 
    for (int v = 0; v < G->get_n_vertices(); ++v)
        add_fictional_column(v);

#elif INITIAL_COLUMN_STRATEGY == 1

    if (father == NULL)
        // Root node: add a fictional column for each vertex 
        for (int v = 0; v < G->get_n_vertices(); ++v)
            add_fictional_column(v);
    
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

    // Apply the stable set covering heuristic
    std::vector<std::list<nodepntArray>> stable_sets;
    G->coloring_heuristic(stable_sets);

    // Mark the colored vertices
    std::vector<bool> colored (G->get_n_vertices(), false);

    for (int i = 0; i < G->get_n_colors(); ++i)
        for (auto &s: stable_sets[i]) {
            // Build the column
            bool *column = NULL;
            G->stable_to_column(s.list, s.n_list, &column);

            // Add the column to the son
            add_real_column(column, s.n_list, i);

            // Mark the colored vertices
            for (int j = 0; j < G->get_n_vertices(); ++j)
                if (column[j]) {
                    colored[j] = true;
                }

            // Free the stable set
            free(s.list);
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
#ifndef VISUALC
	cplex.setParam(IloCplex::IntParam::ClockType, 1); /* set user time */
#endif
	cplex.setParam(IloCplex::NumParam::WorkMem, 2048);
	cplex.setParam(IloCplex::IntParam::Threads, 1);
	cplex.setParam(IloCplex::IntParam::RandomSeed, 1);
	cplex.extract(Xmodel);
#ifdef SAVELP
	cplex.exportModel(SAVELP);
    show("Linear relaxation formulation saved");
#endif

    // Solve

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

void LP::save_solution(std::vector<int> &coloring, std::set<int> &active_colors) {

#if BRANCHING_STRATEGY == 0

    // Build the coloring of the current graph
    std::vector<int> stables_per_color (G->get_n_colors(), 0);
    std::vector<int> temp_coloring (G->get_n_vertices());
    for (int i: pos_vars) {
        int color = vars[i].color;
        int true_color = G->get_C(color, stables_per_color[color]);
        stables_per_color[color]++;
        for (int j = 0; j < G->get_n_vertices(); ++j)
            if (vars[i].stable[j])
                temp_coloring[j] = true_color;
        active_colors.insert(true_color);
    }

    // Build the coloring of the original graph
    coloring.resize(G->get_n_total_vertices());
    for (int i = 0; i < G->get_n_total_vertices(); ++i)
        coloring[i] = temp_coloring[G->get_current_vertex(i)];

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
    branch_on_indistinguishable_colors(ret);
#else
    bye("Undefined branching strategy");
#endif

    return;

}

void LP::branch_on_edges(std::vector<LP *> &ret) {

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
    // 1) Find the pair (v,k) that maximizes the cardinality of N(v) \cap V[k] such that
    //  -) |L(v)| > 1
    //  -) Exists some stable set S with v \in S and x(S,k) is fractional
    // 2) In case of tie, choose the one with the largest value of st(v,k) where:
    //  -) st(v,k) = st1(v,k) + st2(v,k) + st3(v,k)     if |C[k]| = 1
    //  -) st(v,k) = st1(v,k) + st2(v,k)                if |C[k]| > 1 

    // n_L[v] represents |L(v)| for every vertex v
    std::vector<int> n_L;
    G->get_n_L(n_L);

    // Set of candidates (a set is needed to avoid repeated elements)
    std::set<std::pair<int,int>> candidates; // (v,k) candidates
    int max_neighbours = -1;

    for (int i: pos_vars)
        if (values[i] > EPSILON && values[i] < 1 - EPSILON)
            for (int j = 0; j < G->get_n_vertices(); ++j)
                if (vars[i].stable[j]) {
                    if (n_L[j] <= 1) continue;  // Ignore vertex v with |L(v)| = 1
                    int v = j;
                    int k = vars[i].color;
                    int neighbours = G->get_n_neighbours(v, k);
                    if (neighbours > max_neighbours) {
                        candidates.clear();
                        candidates.emplace(v,k);
                        max_neighbours = neighbours;
                    }
                    else if (neighbours == max_neighbours)
                        candidates.emplace(v,k);
                }

    if (candidates.empty())
        bye("Branching error");
    
    int max_v = -1;
    int max_k = -1;

    if (candidates.size() == 1) {
        max_v = candidates.begin()->first;
        max_k = candidates.begin()->second;

    }

    else {

        // Calculate st(v,k) of every candidate
        std::vector<int> st (candidates.size(), 0);
        for (int i: pos_vars)
            if (values[i] > EPSILON && values[i] < 1 - EPSILON) {            
                int c = 0; // st index
                for (auto &t: candidates) {
                    int v = t.first;
                    int k = t.second;
                    if (vars[i].stable[v]) st[c]++; // Increase st1 or st2 by 1
                    else
                        // Check if there is a neighbor of v
                        if (k == vars[i].color && G->get_n_C(k) == 1) 
                            for (int j = 0; j < G->get_n_vertices(); ++j)
                                if (vars[i].stable[j] && G->is_edge(v,j)) {
                                    st[c]++; // Increase st3 by 1
                                    break;
                                }
                    ++c;
                }
            }

        // Choose max_v and max_k
        int max = -1;
        int c = 0; // st index
        for (auto &t: candidates) {
            int v = t.first;
            int k = t.second;
            if (st[c] > max) {
                max = st[c];
                max_v = v;
                max_k = k;
            }
            c++;
        }
    }

    // Build the graph where v is colored with k
    Graph *G1 = G->choose_color(max_v, max_k);

    // Build LP from G1
    LP *lp1 = new LP(G1, this);

    // Build the graph where v cannot be colored with k
    Graph *G2 = G->remove_color(max_v, max_k);
  
    // Build LP from G2
    LP *lp2 = new LP(G2, this);

    ret.resize(2);
    ret[0] = lp1;   // First choose
    ret[1] = lp2;   // Then remove
    return;

}

void LP::branch_on_indistinguishable_colors(std::vector<LP *> &ret) {

    bye("Not implemented yet");

/*

    // Choose vertex v and color k to branch
    // Let W = {v : v \in V_{k_1} \cap V_{k_2}, for some k_1 != k_2 } 
    // 0) If W = {}, then branch on colors
    // 1) Otherwise, find the most fractional (S,k) with |S| > 2
    // 2) If W \cap S = {}, find the next most fractional one.
    // 3) This process is repeated until W \cap S != {}.
    // 4) After that, choose from W \cap S the vertex v with largest value of $|{ x^{*j}_{S'} \notin \ZZ : v \in S', j \in K}|

    std::vector<std::list<int>> branch_vars (G->get_n_vertices());
    // branch_vars[v]: list of indexes i such that
    //  *) values[i] is fractional
    //  *) v \in vars[i].stable
    // and the index of the most fractional variable is in the front of the list

    for (int i: pos_vars)
        if (values[i] > EPSILON && 
            values[i] < 1 - EPSILON &&
            vars[i].n_stable >= 2) {
            for (int v = 0; v < G->get_n_vertices(); ++v)
                if (vars[i].stable[v]) {
                    if (branch_vars.empty())
                        branch_vars[v].push_front(i);
                    else {
                        int front = branch_vars[v].front();
                        if (std::abs(values[i] - 0.5) < std::abs(values[front] - 0.5))
                            branch_vars[v].push_front(i);
                        else
                            branch_vars[v].push_back(i);
                    }
                }

    // Choose the most fractional branching variable
    double min_value = 2.0;
    int branch_vertex;
    int branch_color;
    for (int v = 0; v < G->get_n_vertices(); ++v) {
        if (branch_vars[v].size() > 1) { // Y si son todos estables asociados al mismo color?
            int i = branch_vars[v].front();
            if (std::abs(values[i] - 0.5) < min_value ||
               (std::abs(values[i] - 0.5) == min_value && branch_vars[branch_vertex].size() < branch_vars[v].size())) {   
                min_value = std::abs(values[i] - 0.5);
                branch_vertex = v;
                branch_color = vars[i].color;
            }
        }
    }

    if (min_value > 1)
        branch_on_colors(ret);
    else {

        // Build the graph where v is colored with k
        Graph *G1 = G->choose_indistinguishable_color(branch_vertex, branch_color);

        // Build LP from G1
        LP *lp1 = new LP(G1, this);

        // Build the graph where v cannot be colored with k
        Graph *G2 = G->remove_indistinguishable_color(branch_vertex, branch_color);

        // Build LP from G2
        LP *lp2 = new LP(G2, this);

        ret.resize(2);
        ret[0] = lp1;   // First choose
        ret[1] = lp2;   // Then remove

    }

*/

    return;

}