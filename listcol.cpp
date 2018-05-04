/*
 * LISTCOL - Solves the Minimum Cost List Coloring Problem
 * Made in 2018 by Daniel Severin and Mauro Lucci
 *
 * Requires IBM ILOG CPLEX 12.7
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>

#include<vector>
#include<list>
#include<numeric> // iota function
#include<algorithm> // sort function

/* for linux users: do not define VISUALC */
#ifndef VISUALC
#include <unistd.h>
#include <sys/times.h>
#include "mwis_sewell/mwss.h"
#else
#include <windows.h>
extern "C" {
#include "mwis_sewell/mwss.h"
}
#endif

ILOSTLBEGIN

using namespace std;

/* CONSTANTS */

#define EPSILON 0.00001
#define INTFACTOR 10000.0
#define INFDIST 9999999
#define MAXTIME 7200.0
#define MAXTIME_MWIS 300.0
#define VERBOSE
#define RANDOMCOSTS
//#define SHOWINSTANCE
//#define SHOWALLSTABLES
//#define SHOWSOLUTION
//#define SHOWCPLEX
//#define SAVELP "form.lp"

/* FLAGS OF THE OPTIMIZATION */

//#define PUREBYB
#define USEBARRIER
#define NOMEMEMPHASIS
#define ONLYRELAXATION

/* GLOBAL VARIABLES */

int vertices, edges; /* number of vertices and edges */
int *edge_u, *edge_v; /* array of endpoints of edges */
int *degrees; /* degree of each vertex */
int **neigh_vertices; /* neighbors of each vertex */
int **antineigh_vertices; /* anti-neighbors of each vertex */
int **adjacency; /* adjacency matrix: 0 means no adjacency; >0 gives the index to the edge + 1 */
int **dist; /* distance matrix */

int colors; /* number of colors */
int *cost; /* cost of each color */
int *L_size; /* size of each set L(v), per vertex */
int **L_set; /* elements of L(v), per vertex */
int *C_size; /* size of each set C(k) = L^-1(k), per color */
int **C_set; /* elements of C(k), per color */

int *optimal_coloring; /* optimal solution given by CPLEX */

int sizeP, sizeR, sizeX, sizeV; /* size of the sets for the enumeration algorithm (Bron–Kerbosch) */
bool *setP, *setR, *setX, *setV; /* elements of the sets for the enumeration algorithm */
int current_color; /* color representing stable sets being generated */

IloEnv Xenv; /* CPLEX environment structure */
IloModel Xmodel(Xenv); /* CPLEX model */
IloObjective Xobj; /* CPLEX objective function */
IloNumVarArray Xvars(Xenv); /* CPLEX variables */
IloRangeArray Xrestr(Xenv); /* CPLEX constraints */

MWSSgraph *Mgraph; /* subgraphs of G (one per color) for the MWSS algorithm */
wstable_parameters Mparms;  /* parameters for the MWSS algorithm */
int **V_to_Vk; /* V_to_Vk[k] : V --> V' given the vertex in G, it returns the one in the subgraph k (i.e. C_set[k][V_to_Vk[k][v]] = v) */

/* FUNCTIONS */

/*
 * ECOclock - get a timestamp (in seconds)
 */
double ECOclock() {
#ifndef VISUALC
	/* measure user-time: use it on single-threaded system in Linux (more accurate) */
	struct tms buf;
	times(&buf);
	return ((double)buf.tms_utime) / (double)sysconf(_SC_CLK_TCK);
#else
	/* measure standard wall-clock: use it on Windows */
	return ((double)clock()) / (double)CLOCKS_PER_SEC;
#endif
}

/*
 * set_color - change color of text
 */
void set_color(int color)
{
#ifdef VERBOSE
#ifdef VISUALC
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), color);
#endif
#endif
}

/*
 * bye - finish executing and show a message
 */
void bye(char *string)
{
#ifdef VERBOSE
	set_color(12);
	cout << string << endl;
	set_color(7);
#endif
	exit(1);
}

/*
* urnd - generate numbers with uniform random distribution
*  if flag=false, the interval is [a, b]
*  if flag=true, the interval is [a, b)
*/
float urnd(float a, float b, bool flag)
{
	return a + rand() * (b - a) / (float)(RAND_MAX + (flag ? 1 : 0));
}

/*
 * read_graph - read a graph in the following format:
 *   in the first line, the number of vertices and edges separated by a colon ":"
 *   then, for each line, the endpoints of an edge (u,v) where u < v
 *   note that vertices starts from 0, i.e. 0 < v < |V|-1
 *   example for a diamond graph:
 *     4:5
 *     0,1
 *     0,2
 *     0,3
 *     1,2
 *     1,3
 */
void read_graph(char *filename)
{
	/* open file */
	FILE *stream = fopen(filename, "rt");
	if (!stream) bye("Graph file cannot be opened");
	fscanf(stream, "%d:%d\n", &vertices, &edges);

	/* do not accept graph of less than 4 vertices or stable sets */
	if (vertices < 4) bye("Number of vertices out range!");
	if (edges < 1 || edges > vertices*(vertices - 1) / 2) bye("Number of edges out of range!");

	/* ask for memory */
	degrees = new int[vertices];
	adjacency = new int*[vertices];
	for (int u = 0; u < vertices; u++) {
		degrees[u] = 0;
		adjacency[u] = new int[vertices];
		for (int v = 0; v < vertices; v++) adjacency[u][v] = 0;
	}
	edge_u = new int[edges];
	edge_v = new int[edges];

	/* read edges */
	for (int e = 0; e < edges; e++) {
		int u, v;
		fscanf(stream, "%d,%d\n", &u, &v);
		if (u < 0 || u >= v || v >= vertices) {
			cout << "Error reading edge " << e + 1 << "!" << endl;
			bye("Bye!");
		}
		if (adjacency[u][v] != 0) {
			cout << "A repeated edge was found: (" << u << ", " << v << ")" << endl;
			bye("Bye!");
		}
		else {
			degrees[u]++;
			degrees[v]++;
			edge_u[e] = u;
			edge_v[e] = v;
			adjacency[u][v] = e + 1;
			adjacency[v][u] = e + 1;
		}
	}
	fclose(stream);

	/* also neighborhoods and anti-neighborhoods are computed */
	neigh_vertices = new int*[vertices];
	antineigh_vertices = new int*[vertices];
	for (int v = 0; v < vertices; v++) {
		int degree = degrees[v];

		/* ask for more memory and fill it */
		neigh_vertices[v] = new int[degree];
		antineigh_vertices[v] = new int[vertices - 1 - degree];
		int d = 0, ad = 0;
		for (int w = 0; w < vertices; w++) {
			if (w != v) {
				if (adjacency[v][w] > 0) {
					neigh_vertices[v][d] = w;
					d++;
				}
				else {
					antineigh_vertices[v][ad] = w;
					ad++;
				}
			}
		}
		if (d != degree || ad != vertices - 1 - degree) bye("Internal error!");
	}
}

/*
 * read_cost - read costs of colors in the following format:
 *   in the first line, the number of colors
 *   then, the costs of color 0, 1, etc, in succesive order
 *   example for C = {0, 1, 2} with costs c_0 = 5, c_1 = 2, c_2 = 8:
 *     3
 *     5 2 8
 */
void read_cost(char *filename)
{
	/* open file */
	FILE *stream = fopen(filename, "rt");
	if (!stream) bye("Cost file cannot be opened");
	fscanf(stream, "%d\n", &colors);

	/* do not accept less than 2 colors */
	if (colors < 2) bye("Number of colors out range!");

	/* ask for memory */
	cost = new int[colors];

	/* read costs */
	for (int k = 0; k < colors; k++) {
		int ck;
		fscanf(stream, "%d", &ck);
		if (ck < 0) bye("Color cost must be non negative!");
		cost[k] = ck;
	}
	fclose(stream);
}

/*
 * read_list - read list of colors per vertex in the following format:
 *   in the first line, the number of vertices and colors separated by a colon ":"
 *   then, for each line, the cardinal of L(v) followed by the elements of L(v) in increasing order
 *   example for |V| = 3, |C| = 5 and L(0) = {1, 2}, L(1) = {0, 2, 3}, L(2) = {0, 1, 4}:
 *     3:5
 *     2  1 2
 *     3  0 2 3
 *     3  0 1 4
 */
void read_list(char *filename)
{
	/* open file */
	FILE *stream = fopen(filename, "rt");
	if (!stream) bye("Cost file cannot be opened");

	/* check vv and cc */
	int vv, cc;
	fscanf(stream, "%d:%d\n", &vv, &cc);
	if (vv != vertices) bye("Number of vertices mismatch!");
	if (cc != colors) bye("Number of colors mismatch!");

	/* ask for memory */
	L_size = new int[vertices];
	C_size = new int[colors];
	L_set = new int*[vertices];
	C_set = new int*[colors];

	for (int k = 0; k < colors; k++) C_size[k] = 0;

	/* read lists */
	for (int v = 0; v < vertices; v++) {
		int list_size;
		fscanf(stream, "%d", &list_size);
		if (list_size < 1 || list_size > colors) bye("Error reading lists!");
		L_size[v] = list_size;
		L_set[v] = new int[list_size];

		/* read a list */
		int last_read = -1;
		for (int s = 0; s < list_size; s++) {
			int element;
			fscanf(stream, "%d", &element);
			if (element <= last_read || element >= colors) bye("Error reading lists!");
			last_read = element;
			L_set[v][s] = element;
			C_size[element]++;
		}
	}
	fclose(stream);

	/* now, fill the inverse of L, i.e. C_set */
	for (int k = 0; k < colors; k++) {
		int set_size = C_size[k];
		if (set_size == 0) {
			set_color(13);
			cout << "Warning: Color " << k << " does not belong to any list." << endl;
			set_color(7);
			continue;
		}
		C_set[k] = new int[set_size];
		int p = 0;
		/* find vertices v such that k in L(v), and add v to C_set */
		for (int v = 0; v < vertices; v++) {
			for (int s = 0; s < L_size[v]; s++) {
				if (L_set[v][s] == k) C_set[k][p++] = v;
			}
		}
		if (p != set_size) bye("Internal error!");
	}
}

/*
 * connected() - check if G is a connected graph
 * in addition, we compute a matrix of distances between vertices
 */
bool connected()
{
	/* we ask for memory and fill the distance matrix with 0 in the diagonal, 1 for neighbors and +inf for the remaining entries */
	dist = new int*[vertices];
	for (int u = 0; u < vertices; u++) {
		dist[u] = new int[vertices];
		for (int v = 0; v < vertices; v++) {
			int d = INFDIST;
			if (u == v)	d = 0;
			else { if (adjacency[u][v] > 0) d = 1; }
			dist[u][v] = d;
		}
	}

	/* we use a simple implementation of Floyd algorithm (note: it is not the best way of knowing if G is connected!) */
	for (int v = 0; v < vertices; v++) {
		for (int u = 0; u < vertices; u++) {
			for (int w = 0; w < vertices; w++) {
				int sum = dist[u][v] + dist[v][w];
				if (sum < dist[u][w]) dist[u][w] = sum;
			}
		}
	}

	/* compute diameter */
	int diameter = 0;
	for (int u = 0; u < vertices - 1; u++) {
		for (int v = u + 1; v < vertices; v++) {
			int d = dist[u][v];
			if (d >= INFDIST) {
				cout << "There is no path between " << u << " and " << v << "." << endl;
				return false;
			}
			if (diameter < d) diameter = d;
		}
	}
	cout << "Diameter of G: " << diameter << endl;

	return true;
}

/*
 * enumerate_stable_sets - Apply the Bron-Kerbosch algorithm
 * Here, we use the Tomita variant as in:
 *   http://www.dcs.gla.ac.uk/~pat/jchoco/clique/enumeration/report.pdf
 */
void enumerate_stable_sets()
{
	// cout << "R = {";
	// for (int v = 0; v < vertices; v++) if (setR[v]) cout << " " << v;
	// cout << " } (" << sizeR << "),  P = {";
	// for (int v = 0; v < vertices; v++) if (setP[v]) cout << " " << v;
	// cout << " } (" << sizeP << "),  X = {";
	// for (int v = 0; v < vertices; v++) if (setX[v]) cout << " " << v;
	// cout << " } (" << sizeX << ")" << endl;

	/* propose R if |P|=|X|=0 */
	if (sizeP == 0) {
		if (sizeX == 0) {
			/* R is a maximal stable set of G[C(k)] */
			IloNumColumn stable_set = Xobj(cost[current_color]);
			/* fill the column corresponding to ">= 1" constraints (insert "1" in constraint indexed by v) */
			for (int v = 0; v < vertices; v++) if (setR[v]) stable_set += Xrestr[v](1.0);
			/* and the ">= -1 constraint (insert "-1" in constraint indexed by color) */
			stable_set += Xrestr[vertices + current_color](-1.0);
#ifdef ONLYRELAXATION
			/* add the column as a non-negative continuos variable */
			Xvars.add(IloNumVar(stable_set));
#else
			/* add the column as a non-negative integer variable */
			Xvars.add(IloIntVar(stable_set));
#endif

#ifdef SHOWALLSTABLES
			set_color(4);
			cout << "Set found: {";
			for (int v = 0; v < vertices; v++) if (setR[v]) cout << " " << v;
			cout << " }, size = " << sizeR << " (color = " << current_color << ")" << endl;
			set_color(7);
#endif
		}
		return;
	}
	/* choose u from P union X such that u has the highest number of anti-neighbors in P */
	int best_u = -1, max_antineighbors = -1;
	for (int u = 0; u < vertices; u++) {
		if (setP[u] || setX[u]) {
			/* count the number of anti-neighbors in P */
			int antidegree = vertices - 1 - degrees[u];
			int antineighbors = 0;
			for (int d = 0; d < antidegree; d++) {
				int w = antineigh_vertices[u][d];
				if (setP[w]) antineighbors++;
			}
			if (antineighbors > max_antineighbors) {
				best_u = u;
				max_antineighbors = antineighbors;
			}
		}
	}

	/* backup P and X */
	bool *setPP = new bool[vertices];
	bool *setXX = new bool[vertices];
	for (int w = 0; w < vertices; w++) {
		setPP[w] = setP[w];
		setXX[w] = setX[w];
	}
	int sizePP = sizeP;
	int sizeXX = sizeX;

	/* check every v in P - AN(u), or equivalently v in P cap N[u] (since P is in C(k), v is a vertex of G[C(k)]) */
	int *cap = new int[sizeP];
	int cap_size = 0;
	if (setP[best_u]) cap[cap_size++] = best_u;
	for (int d = 0; d < degrees[best_u]; d++) {
		int v = neigh_vertices[best_u][d];
		if (setP[v]) cap[cap_size++] = v;
	}

	for (int d = 0; d < cap_size; d++) {
		int v = cap[d];
		/* Tomita call: R <- R union {v},  P <- P cap AN(v),   X <- X cap AN(v) */
		/* or equivalently P <- P - N[v],   X <- X - N[v] */
		if (setR[v]) bye("Internal error!");
		setR[v] = true;
		sizeR++;
		for (int w = 0; w < vertices; w++) {
			if (setV[w] && (w == v || adjacency[v][w] > 0)) {
				if (setP[w]) {
					setP[w] = false;
					sizeP--;
				}
				if (setX[w]) {
					setX[w] = false;
					sizeX--;
				}
			}
		}
		enumerate_stable_sets();

		/* restore R, P and X */
		setR[v] = false;
		sizeR--;
		for (int w = 0; w < vertices; w++) {
			setP[w] = setPP[w];
			setX[w] = setXX[w];
		}
		sizeP = sizePP;
		sizeX = sizeXX;

		/* pass processed vertices from P to X (except the last iteration which is not used) */
		if (d < cap_size - 1) {
			for (int s = 0; s <= d; s++) {
				int w = cap[s];
				if (!setP[w] || setX[w]) bye("Internal error!");
				setP[w] = false;
				setX[w] = true;
			}
			sizeP -= d + 1;
			sizeX += d + 1;
		}
	}

	/* restore P and X */
	for (int w = 0; w < vertices; w++) {
		setP[w] = setPP[w];
		setX[w] = setXX[w];
	}
	sizeP = sizePP;
	sizeX = sizeXX;

	delete[] cap;
	delete[] setXX;
	delete[] setPP;
}

/*
 * print_solution - show the obtained solution on screen
 */
void print_solution(IloCplex cplex)
{
#ifdef SHOWSOLUTION
	cout << "Solution (non-zero variables): " << endl;
	int count = 0;
	for (int x = 0; x < Xvars.getSize(); x++) {
		IloNumVar var1 = Xvars[x];
		double xval = cplex.getValue(var1);
		if (xval > EPSILON) {
			cout << "  Stable {";
			for (int v = 0; v < vertices; v++) {
				IloRange expr = Xrestr[v];
				for (IloExpr::LinearIterator it = expr.getLinearIterator(); it.ok(); ++it) {
					IloNumVar var2 = it.getVar();
					IloNum coef2 = it.getCoef();
					if (var2.getId() == var1.getId() && coef2 > 0.5) {
						cout << " " << v;
						break;
					}
				}
			}
			cout << " }, ";
			/* find the color associated to that stable */
			int color_chosen = -1;
			for (int k = 0; k < colors; k++) {
				IloRange expr = Xrestr[vertices + k];
				for (IloExpr::LinearIterator it = expr.getLinearIterator(); it.ok(); ++it) {
					IloNumVar var2 = it.getVar();
					IloNum coef2 = it.getCoef();
					if (var2.getId() == var1.getId() && coef2 < -0.5) {
						color_chosen = k;
						goto colors_was_chosen;
					}
				}
			}
			if (color_chosen == -1) bye("No color chosen for that stable!");
colors_was_chosen:;
			cout << "k = " << color_chosen << ", x*(" << x << ") = " << xval << endl;
			count++;
		}
	}
	cout << "  count = " << count << endl;
#endif
}

/*
 * optimize1 - make an exhaustive enmeration of stable sets and solve the set-cover formulation
 */
bool optimize1()
{
	setP = new bool[vertices];
	setR = new bool[vertices];
	setX = new bool[vertices];
	setV = new bool[vertices];

	Xobj = IloMinimize(Xenv);

	/* we will have "vertices" constraints with r.h.s >= 1 and "colors" constraints with r.h.s >= -1 */
	for (int v = 0; v < vertices; v++) Xrestr.add(IloRange(Xenv, 1.0, IloInfinity));
	for (int k = 0; k < colors; k++) Xrestr.add(IloRange(Xenv, -1.0, IloInfinity));

	/* enumerate all the stable sets of the graph G[C(k)] with Bron–Kerbosch algorithm */
	for (int k = 0; k < colors; k++) {
		sizeP = 0;
		sizeR = 0;
		sizeX = 0;
		sizeV = 0;
		for (int v = 0; v < vertices; v++) {
			setP[v] = false;
			setR[v] = false;
			setX[v] = false;
			setV[v] = false;
		}
		for (int s = 0; s < C_size[k]; s++) {
			int v = C_set[k][s];
			setP[v] = true; sizeP++;
			setV[v] = true; sizeV++;
		}
		current_color = k;
		enumerate_stable_sets();
	}

	Xmodel.add(Xobj);
	Xmodel.add(Xrestr);

	IloCplex cplex(Xmodel);
	cplex.setDefaults();
#ifndef SHOWCPLEX
	cplex.setOut(Xenv.getNullStream());
	cplex.setWarning(Xenv.getNullStream());
#endif
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
#ifdef USEBARRIER
	cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Algorithm::Barrier);
#else
	cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Algorithm::Primal);
#endif
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

	/* solve it! */
	cplex.solve();
	IloCplex::CplexStatus status = cplex.getCplexStatus();

#ifdef ONLYRELAXATION
	/* LP treatment */
	if (status != IloCplex::Optimal) {
		switch (status) {
		case IloCplex::InfOrUnbd:
		case IloCplex::Infeasible: bye("Infeasible :(");
		case IloCplex::AbortTimeLim: cout << "Time limit reached!" << endl; break;
		default: bye("Unexpected error :(");
		}
		return false;
	}
#else
	/* MIP treatment */
	IloInt nodes = cplex.getNnodes();
	cout << "Number of nodes evaluated: " << nodes << endl;
	if (status != IloCplex::Optimal) {
		switch (status) {
		case IloCplex::InfOrUnbd:
		case IloCplex::Infeasible: bye("Infeasible :(");
		case IloCplex::AbortTimeLim: cout << "Time limit reached!" << endl; break;
		default: bye("Unexpected error :(");
		}
		/* read floating point bounds and convert them to integer values */
		double lower = cplex.getObjValue();
		double upper = cplex.getBestObjValue();
		cout << "Best bounds are " << (int)(lower + (1.0 - EPSILON)) << " <= objective <= " << (int)(upper + EPSILON) << "." << endl;
		cout << "Relative gap = " << 100.0 * (upper - lower) / lower << "." << endl; /* Rel Gap = |bestbound - bestinteger|/|bestinteger| */
		return false;
	}

	/* save optimal solution */
	for (int v = 0; v < vertices; v++) {
		int color_chosen = -1;
		IloRange expr1 = Xrestr[v];
		for (IloExpr::LinearIterator it1 = expr1.getLinearIterator(); it1.ok(); ++it1) {
			IloNumVar var1 = it1.getVar();
			IloNum coef1 = it1.getCoef();
			if (cplex.getValue(var1) > 0.5 && coef1 > 0.5) {
				for (int k = 0; k < colors; k++) {
					IloRange expr2 = Xrestr[vertices + k];
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
		optimal_coloring[v] = color_chosen;
	}
#endif

	/* optimality reached, show stables */
	set_color(10);
	cout << "Optimality reached! :)" << endl;
	set_color(3);
	print_solution(cplex);
	set_color(7);
	cout << "objective = " << cplex.getObjValue() << endl;
	cout << "variables = " << cplex.getNcols() << endl;

	delete[] setV;
	delete[] setX;
	delete[] setR;
	delete[] setP;
	return true;
}

/*
 * initialize_MWSS - perform some initializations for the MWSS algorithm, including the generation of the subgraphs
 */
void initialize_MWSS()
{
	MWSSdata Mdata;
	wstable_info Minfo;

	default_parameters(&Mparms);
	Mparms.cpu_limit = MAXTIME_MWIS;
	// Mparms.reorder = 1; /* reorder vertices of subgraph by degree */
	// if (density < 0.2) Mparms.clique_cover = 2;
	// else Mparms.clique_cover = 1;

	Mgraph = new MWSSgraph[colors];
	V_to_Vk = new int*[colors];
	for (int k = 0; k < colors; k++) {
		reset_pointers(&Mgraph[k], &Mdata, &Minfo);

		/* set the lookup function V -> Vk */
		int vertices_Gk = C_size[k];
		V_to_Vk[k] = new int[vertices];
		for (int v = 0; v < vertices; v++) V_to_Vk[k][v] = -1;
		for (int s = 0; s < vertices_Gk; s++) {
			int v = C_set[k][s];
			V_to_Vk[k][v] = s;
		}

		/* generate the graph Gk: recall that the first entry (i.e. array[0]) is not used!!! */
		allocate_graph(&Mgraph[k], vertices_Gk);
		for (int v = 1; v <= vertices_Gk; v++) {
			Mgraph[k].node_list[v].name = v;
			Mgraph[k].node_list[v].degree = 0;
		}
		for (int row = 1; row <= vertices_Gk; row++) {
			int u = C_set[k][row - 1]; /* recall the "0" is not used! */
			Mgraph[k].adj[row][row] = 0;
			for (int col = row + 1; col <= vertices_Gk; col++) {
				int v = C_set[k][col - 1]; /* recall the "0" is not used! */
				int val = adjacency[u][v] > 0 ? 1 : 0;
				Mgraph[k].adj[row][col] = val;
				Mgraph[k].adj[col][row] = val;
			}
		}
		build_graph(&Mgraph[k]);
	}
}

/*
 * solve_MWSS - solves the Maximum Weighted Stable Set Problem of a subgraph of G with vector cost pi >= 0
 *              and returns the size of the stable set such that: 1) has maximum weight or 2) is greater than "goal" in its cost
 */
int solve_MWSS(int k, double *pi, double goal, int *stable_set, double *weight)
{
	MWSSdata Mdata;
	wstable_info Minfo;

	/* perform optimization */
	for (int i = 1; i <= C_size[k]; i++) Mgraph[k].weight[i] = (int)(INTFACTOR*pi[i - 1]); /* recall that "0" is not used! */
	if (initialize_max_wstable(&Mgraph[k], &Minfo) > 0) bye("Failed in initialize_max_wstable");
	if (call_max_wstable(&Mgraph[k], &Mdata, &Mparms, &Minfo, (int)(INTFACTOR*goal), 0) > 0) bye("Failed in call_max_wstable");
	//	if (Mparms.cpu_limit >= 0 && Minfo.cpu > Mparms.cpu_limit) printf("cpu_limit of %f seconds exceeded: %f seconds. Solution may not optimum.\n", Mparms.cpu_limit, Minfo.cpu);
	// printf("Found best stable set of weight %d.\n", Mdata.best_z);

	/* save best stable set (in terms of the vertices in G) */
	int stable_set_size = Mdata.n_best;
	for (int i = 1; i <= stable_set_size; i++) stable_set[i - 1] = Mdata.best_sol[i]->name - 1; /* recall that "0" is not used! */
	*weight = ((double)Mdata.best_z) / INTFACTOR;

	free_data(&Mdata);
	free_wstable_info(&Minfo);
	return stable_set_size;
}

/*
 * cost_comparison_function - compare a pair of colors by their cost
*/
bool cost_comparison_function(int i, int j) {
	return (cost[i] < cost[j]);
}

/*
 * stable_covering_heuristic - construct a maximal stable covering
 * 		             - if succeeds, adds the corresponding columns to the master problem
*/
void stable_covering_heuristic() {

	// Sort colors in non-decreasing order
	vector<int> ordered_colors(colors);
	iota(ordered_colors.begin(), ordered_colors.end(), 0); // initialization
	sort(ordered_colors.begin(), ordered_colors.end(), cost_comparison_function);

	int covered = 0; // number of covered vertices
	vector<bool> covered_array(vertices, false); // at first every vertex is uncovered

	// Foreach color k try to find a maximal stable set in G[C(k)] prioritizing non-covered verteces
        int counter = 0;
	for (int k : ordered_colors) {

		// Termination criteria checking
		if (covered == vertices)
			break;

		// Sort C[k] in such a way that non-covered vertices appear first
		list<int> ordered_vertices;
		bool flag = false; // flag = true if G[C(k)] has at least one non-covered vertex
		for (int i = 0; i < C_size[k]; i++) {
			if (covered_array[C_set[k][i]])
				ordered_vertices.push_back(C_set[k][i]); // Push it at back
			else {
				ordered_vertices.push_front(C_set[k][i]); // Push it at front
				flag = true;
			}
		}

		if (!flag)
			continue;

		// Maximal stable construction
		vector<int> stable;
		stable.reserve(C_size[k]);
		while (ordered_vertices.size() != 0) {

			int v = ordered_vertices.front();
			ordered_vertices.pop_front(); // Add v to the stable set
			stable.push_back(v);
			if (!covered_array[v]) {
                                covered_array[v] = true;
                                covered++;
			}

			// Delete neighbors of v
			// TODO: rewrite in case G[C(k)] is stored
			for (int i = 0; i < degrees[v]; i++) {
				ordered_vertices.remove(neigh_vertices[v][i]); // Time consuming implementation
			}

		}

		// Add column
		IloNumColumn column = Xobj(cost[k]);
		// fill the column corresponding to ">= 1" constraints (insert "1" in constraint indexed by v)
		for (int v : stable)
			column += Xrestr[v](1.0);
		// and the ">= -1 constraint (insert "-1" in constraint indexed by color)
		column += Xrestr[vertices + k](-1.0);

		/* add the column as a non-negative continuos variable */
		Xvars.add(IloNumVar(column));
                counter++;

	}

	// Check if there is a remaning non-covered vertex
	if (covered != vertices)
		bye("Heuristic failed to find an initial maximal stable covering");
        else
                cout << counter << " columns were added"<< endl;

}

/*
 * optimize2 - solve the set-cover formulation via column generation
 */
bool optimize2()
{
#ifndef ONLYRELAXATION
	bye("Branch-and-price not implemented yet!");
#endif
	/* generate the subgraphs Gk */
	initialize_MWSS();

	Xobj = IloMinimize(Xenv);

	/* we will have "vertices" constraints with r.h.s >= 1 and "colors" constraints with r.h.s >= -1 */
	for (int v = 0; v < vertices; v++) Xrestr.add(IloRange(Xenv, 1.0, IloInfinity));
	for (int k = 0; k < colors; k++) Xrestr.add(IloRange(Xenv, -1.0, IloInfinity));

	stable_covering_heuristic();

	Xmodel.add(Xobj);
	Xmodel.add(Xrestr);

	// Settings

	IloCplex cplex(Xmodel);
	cplex.setDefaults();
#ifndef SHOWCPLEX
	cplex.setOut(Xenv.getNullStream());
	cplex.setWarning(Xenv.getNullStream());
#endif
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
#ifdef USEBARRIER
	cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Algorithm::Barrier);
#else
	cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Algorithm::Primal);
#endif
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

	// Column generation approach
  	int *stable_set = new int[vertices];
        double stable_weight = 0.0; 
        double *pi = new double[vertices];
	while (true) {

                cplex.solve();

                // Find an entering column
                bool exists = false;
                IloNumArray duals (Xenv,vertices+colors);
                cplex.getDuals(duals, Xrestr); // Save dual values

                int counter = 0;

                for (int k = 0; k < colors; k++) {
                        for (int i = 0; i < C_size[k]; i++) {
                                double d = duals[C_set[k][i]];
                                if (d > -EPSILON && d < 0.0) d = 0.0;
                                pi[i] = d;
                        }
                        //for (int i = 0; i < C_size[k]; i++) cout << "pi[" << C_set[k][i] << "] = " << pi[i] << endl;
                        double goal = cost[k] + duals[vertices + k];
                        //cout << "goal = " << goal << endl;

                        int stable_set_size = solve_MWSS(k, pi, goal, stable_set, &stable_weight);

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

                        if (stable_weight - goal > EPSILON) {
                                // Add column
                                exists = true;
                                //cout << "Stable set = {";
                                //for (int i = 0; i < stable_set_size; i++) cout << " " << C_set[k][stable_set[i]];
                                //cout << " }" << endl;
                                IloNumColumn column = Xobj(cost[k]);
                                // fill the column corresponding to ">= 1" constraints
                                for (int i = 0; i < stable_set_size; i++) column += Xrestr[C_set[k][stable_set[i]]](1.0);
                                // and the ">= -1 constraint
                                column += Xrestr[vertices + k](-1.0);

                                /* add the column as a non-negative continuos variable */
                                Xvars.add(IloNumVar(column));
                                counter++;

                                //break; // comment this line for multiple column addition per iteration
                        }
                }

                cout << counter << " columns were added"<< endl;

                if (!exists)
                        break; // optimality reached
        }

        IloCplex::CplexStatus status = cplex.getCplexStatus();
#ifdef ONLYRELAXATION
        /* LP treatment */
        if (status != IloCplex::Optimal) {
	        switch (status) {
	        case IloCplex::InfOrUnbd:
	        case IloCplex::Infeasible: bye("Infeasible :(");
	        case IloCplex::AbortTimeLim: cout << "Time limit reached!" << endl; break;
	        default: bye("Unexpected error :(");
	        }
	        return false;
        }
        else {
			/* optimality reached, show stables */
			set_color(10);
			cout << "Optimality reached! :)" << endl;
			set_color(3);
			print_solution(cplex);
			set_color(7);
			cout << "objective = " << cplex.getObjValue() << endl;
			cout << "variables = " << cplex.getNcols() << endl;
		}
#endif

	/* free memory */
	delete[] pi;
	delete[] stable_set;
	for (int k = 0; k < colors; k++) {
		free_graph(&Mgraph[k]);
		delete[] V_to_Vk[k];
	}
	delete[] V_to_Vk;
	delete[] Mgraph;
	return true;
}

/*
 * check_coloring - check if the given coloring is valid
 */
bool check_coloring(int *f) {
	/* colors of each vertex belong to the list ??? */
	for (int v = 0; v < vertices; v++) {
		int k_chosen = f[v];
		bool happy = false;
		for (int s = 0; s < L_size[v]; s++) {
			if (L_set[v][s] == k_chosen) {
				happy = true;
				break;
			}
		}
		if (!happy) {
			set_color(12);
			cout << "Ouch! Vertex " << v << " is colored with " << k_chosen << " which does not belong to its list :S" << endl;
			set_color(7);
			return false;
		}
	}

	/* are there conflicting edges ??? */
	for (int e = 0; e < edges; e++) {
		int u = edge_u[e];
		int v = edge_v[e];
		if (f[u] == f[v]) {
			set_color(12);
			cout << "Ouch! Vertices " << u << " and " << v << " have the same color :S" << endl;
			set_color(7);
			return false;
		}
	}
	return true;
}

/*
 * main - Main program
 */
int main(int argc, char **argv)
{
#ifndef VERBOSE
	std::cout.setstate(std::ios::failbit);
#endif

	set_color(15);
	cout << "LISTCOL - Solves the Minimum Cost List Coloring Problem." << endl;
	set_color(7);

	if (argc < 2 || argc == 3 || argc > 4) {
		cout << "Usage: listcol file.graph [file.cost file.list]" << endl;
		cout << "  if cost and list file are not given, then it solves the classic coloring problem" << endl;
		bye("Bye!");
	}

	/* read graph */
	read_graph(argv[1]);

	if (argc == 4) {
		/* read colors and its costs */
		read_cost(argv[2]);

		/* read list of colors per vertex */
		read_list(argv[3]);
	}
	else {
		/* set C = V, its costs as 1 (or random from {1,...,10}) and L(v) = C */
		colors = vertices;
		cost = new int[vertices];
		for (int k = 0; k < colors; k++) {
#ifdef RANDOMCOSTS
			cost[k] = (int)urnd(1, 10, false);
#else
			cost[k] = 1;
#endif
		}
		L_size = new int[vertices];
		for (int v = 0; v < vertices; v++) L_size[v] = colors;
		C_size = new int[colors];
		for (int k = 0; k < colors; k++) C_size[k] = vertices;
		L_set = new int*[vertices];
		for (int v = 0; v < vertices; v++) {
			L_set[v] = new int[colors];
			for (int k = 0; k < colors; k++) L_set[v][k] = k;
		}
		C_set = new int*[colors];
		for (int k = 0; k < colors; k++) {
			C_set[k] = new int[vertices];
			for (int v = 0; v < vertices; v++) C_set[k][v] = v;
		}
	}

#ifdef SHOWINSTANCE
	set_color(2);
	cout << "Neighborhoods:" << endl;
	int maxdelta = 0;
	for (int v = 0; v < vertices; v++) {
		cout << "N(" << v << ") = {";
		int degree = degrees[v];
		if (degree > maxdelta) maxdelta = degree;
		for (int d = 0; d < degree; d++) cout << " " << neigh_vertices[v][d];
		cout << " }, degree = " << degree << endl;
	}
	cout << "Maximum degree = " << maxdelta << endl;
	cout << "Vector of costs: {";
	for (int k = 0; k < colors; k++) {
		cout << " " << k << "->" << cost[k];
	}
	cout << " }, colors = " << colors << endl;
	cout << "List of colors:" << endl;
	for (int v = 0; v < vertices; v++) {
		cout << "L(" << v << ") = {";
		int list_size = L_size[v];
		for (int s = 0; s < list_size; s++) cout << " " << L_set[v][s];
		cout << " }" << endl;
	}
#endif

	/* Show some basic statistics */
	set_color(6);
	cout << "Statistics:" << endl;
	int clique_size = vertices * (vertices - 1) / 2;
	float density = 100.0 * (float)edges / (float)clique_size;
	cout << "  |V| = " << vertices << ", |E| = " << edges << " (density = " << density << "%), |C| = " << colors << "." << endl;
	/* Average and standard deviation of size of lists */
	float prom = 0.0;
	for (int k = 0; k < vertices; k++) prom += (float)L_size[k];
	prom /= (float)vertices;
	float sigma = 0.0;
	for (int k = 0; k < vertices; k++) {
		float substr = (float)L_size[k] - prom;
		sigma += substr * substr;
	}
	sigma /= (float)(vertices - 1);
	cout << "  Behaviour of |L(v)| ---> prom = " << prom << ", sigma = " << sqrt(sigma) << "." << endl;
	/* Average and standard deviation of vertices of Gk */
	prom = 0.0;
	for (int k = 0; k < colors; k++) prom += (float)C_size[k];
	prom /= (float)colors;
	sigma = 0.0;
	for (int k = 0; k < colors; k++) {
		float substr = (float)C_size[k] - prom;
		sigma += substr * substr;
	}
	sigma /= (float)(colors - 1);
	cout << "  Behaviour of |V(Gk)| ---> prom = " << prom << ", sigma = " << sqrt(sigma) << "." << endl;
	set_color(5);
	if (connected() == false) cout << "Note: G is not connected!" << endl;
	set_color(7);

	/* optimize using an exhaustive enmeration of stable sets */
	optimal_coloring = new int[vertices];
	double start_t = ECOclock();
	bool status = optimize2();
	double stop_t = ECOclock();
	if (status) {
#ifndef ONLYRELAXATION
		/* show the answer on screen */
		cout << "Coloring:" << endl;
		for (int v = 0; v < vertices; v++) cout << "  f(" << v << ") = " << optimal_coloring[v] << endl;
		check_coloring(optimal_coloring);
#endif
	}
	cout << "Time of optimization = " << stop_t - start_t << " sec." << endl;

	/* free memory */
	delete[] optimal_coloring;
	for (int v = 0; v < vertices; v++) delete[] dist[v];
	delete[] dist;
	for (int k = 0; k < colors; k++) delete[] C_set[k];
	delete[] C_set;
	for (int v = 0; v < vertices; v++) delete[] L_set[v];
	delete[] L_set;
	delete[] C_size;
	delete[] L_size;
	delete[] cost;
	for (int v = 0; v < vertices; v++) {
		delete[] antineigh_vertices[v];
		delete[] neigh_vertices[v];
	}
	delete[] antineigh_vertices;
	delete[] neigh_vertices;
	delete[] edge_v;
	delete[] edge_u;
	for (int v = 0; v < vertices; v++) delete[] adjacency[v];
	delete[] adjacency;
	delete[] degrees;

#ifndef VERBOSE
	std::cout.clear();
#endif
	return 0;
}
