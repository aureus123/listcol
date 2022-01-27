/*
 * LISTCOLA - Solves the Minimum Cost List Coloring Problem w/Alternative Models
 * Made in 2018-2019 by Daniel Severin
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
#else
#include <windows.h>
#endif

ILOSTLBEGIN

using namespace std;

/* CONSTANTS */

#define EPSILON 0.00001
#define BIGNUMBER 99999999
#define MAXTIME 7200.0
#define VERBOSE
//#define SHOWSTATS
//#define SHOWINSTANCE
//#define SHOWALLSTABLES
#define SHOWCPLEX
//#define SAVELP "form.lp"

/* FLAGS OF THE OPTIMIZATION */

/* Both models: */
#define SYMMETRYCOL
//#define ONLYRELAXATION
//#define TUNEDPARAMS
//#define STABLEMODEL

/* Set Covering model: */
#define MAXCOLUMNS 2000000

/* Compact Formulation: */
//#define EDGEINEQ
//#define SYMMETRYRESTR1
//#define SYMMETRYRESTR2
//#define SYMMETRYRESTR3

/* GLOBAL VARIABLES */

int vertices, edges; /* number of vertices and edges */
int *edge_u, *edge_v; /* array of endpoints of edges */
int *degrees; /* degree of each vertex */
int **neigh_vertices; /* neighbors of each vertex */
int **antineigh_vertices; /* anti-neighbors of each vertex */
int **adjacency; /* adjacency matrix: 0 means no adjacency; >0 gives the index to the edge + 1 */

int colors; /* number of colors */
int *cost; /* cost of each color */
int *L_size; /* size of each set L(v), per vertex */
int **L_set; /* elements of L(v), per vertex */
int *C_size; /* size of each set C(k) = L^-1(k), per color */
int **C_set; /* elements of C(k), per color */

/* The structure of the partition consists of an array and a linked-list. The linked-list is formed with
   header_part (points to the first) and next_part (points to the next). */
int *header_part; /* is the first color (representative) of a given color */
int *next_part; /* is the next color in the set containing the given color, -1 if the given color is the last */
int *part_set; /* array with the first color of each set in the partition */
int *part_card; /* array with the cardinal of each set */
int part_size; /* size of the array */

int *optimal_coloring; /* optimal solution given by CPLEX */
bool *colors_used; /* colors used in the optimal solution */
IloInt nodes_explored; /* number of explored nodes */
double lowerbound, upperbound; /* LB and UB of optimization */

int sizeP, sizeR, sizeX, sizeV; /* size of the sets for the enumeration algorithm (Bron–Kerbosch) */
bool *setP, *setR, *setX, *setV; /* elements of the sets for the enumeration algorithm */
int variables; /* number of variables generated */
int current_part; /* partition (of colors) representing stable sets being generated */

bool **adj_mat, *S_set, *Q_set; /* graph (as adjacency matrix), input and output sets for the find_max_clique procedure */

IloEnv Xenv; /* CPLEX environment structure */
IloModel Xmodel(Xenv); /* CPLEX model */
IloObjective Xobj; /* CPLEX objective function */
IloNumVarArray Xvars(Xenv); /* CPLEX variables */
IloRangeArray Xrestr(Xenv); /* CPLEX constraints */

/* FUNCTIONS */

/*
 * cmpfunc- for using with qsort (sort "val" in descending order)
 * usage: qsort(data, cardinal, sizeof(sortstr), cmpfunc);
 */
struct sortstr { int i, val; };
int cmpfunc(const void *a, const void *b) {
	return (((struct sortstr*)b)->val - ((struct sortstr*)a)->val);
}

/*
 * ECOclock - get a timestamp (in seconds)
 */
double ECOclock() {
#ifdef VISUALC
	/* measure standard wall-clock: use it on Windows */
	return ((double)clock()) / (double)CLOCKS_PER_SEC;
#else
	/* measure user-time: use it on single-threaded system in Linux (more accurate) */
	struct tms buf;
	times(&buf);
	return ((double)(buf.tms_utime)) / (double)sysconf(_SC_CLK_TCK);
#endif
}

/*
 * set_color - change color of text
 *     0 - Black, 1 - Blue, 2 - Green, 3 - Cyan, 4 - Red, 5 - Purple, 6 - Yellow, 7 - Gray
 *     8-15 - Brighter colors
 */
void set_color(int color)
{
#ifdef VERBOSE
#ifdef VISUALC
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), color);
#else
	const char *codes[] = {
	  "30", "34", "32", "36", "31", "35", "33", "37",
	  "90", "94", "92", "96", "91", "95", "93", "97"
	};
	printf("\e[%sm", codes[color]);	
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
	exit(-1);
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
	if (edges < 0 || edges > vertices*(vertices - 1) / 2) bye("Number of edges out of range!");

	/* ask for memory */	degrees = new int[vertices];
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
		if (list_size < 1 || list_size > colors) {
			cout << "list_size = " << list_size << ", colors = " << colors << endl;
			bye("Error reading lists!");
		}
		L_size[v] = list_size;
		L_set[v] = new int[list_size];

		/* read a list */
		int last_read = -1;
		for (int s = 0; s < list_size; s++) {
			int element;
			fscanf(stream, "%d", &element);
			if (element <= last_read || element >= colors) {
				cout << "last_read = " << last_read << ", element = " << element << ", colors = " << colors << endl;
				bye("Error reading lists!");
			}
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
 * find_color_part - Find a partition of indistinguishable colors
 */
void find_color_part()
{
#ifdef SYMMETRYCOL
	/* Symmetry-col flag turned on -> find those colors that are undistinguishable */
	for (int k = 0; k < colors; k++) {
		bool new_color = true;
		int prev_index; /* <-- index to the set of the partition in case color k is not new */
		if (k > 0) {
			/* check if the color k is the same as a previous color of the partition */
			for (int i = 0; i < part_size; i++) {
				prev_index = i;
				int prev = part_set[i];
				if (cost[k] == cost[prev] && C_size[k] == C_size[prev]) {
					/* if costs are equal and sizes of Vk too, check vertex by vertex
					   (requires that elements of C_set be ordered in ascending order) */
					bool equal = true;
					for (int c = 0; c < C_size[k]; c++) {
						if (C_set[k][c] != C_set[prev][c]) {
							equal = false;
							break;
						}
					}
					if (equal) new_color = false;
				}
				if (new_color == false) break;
			}
		}
		if (new_color) {
			/* a new color k is added to the partition */
			header_part[k] = k;
			next_part[k] = -1;
			part_set[part_size] = k;
			part_card[part_size] = 1;
			part_size++;
		}
		else {
			/* colors k and part_set[prev_index] are indistinguishable from each other  */
			int prev = part_set[prev_index];
			part_card[prev_index]++;
			header_part[k] = prev;
			next_part[k] = -1;
			/* travel the set until reach the last element and append k */
			int r = prev;
			int s;
			do { s = r; r = next_part[s]; } while (r != -1);
			next_part[s] = k;
		}
	}
#else
	/* Symmetry-col flag turned off -> say that every color is indistinguishable */
	for (int k = 0; k < colors; k++) {
		header_part[k] = k;
		next_part[k] = -1;
		part_set[k] = k;
		part_card[k] = 1;
	}
	part_size = colors;
#endif

//#ifdef SHOWINSTANCE
	/* show the partition obtained */
	cout << "Partition of colors:";
	for (int i = 0; i < part_size; i++) {
		int k = part_set[i];
		cout << " {";
		do {
			cout << " " << k;
			k = next_part[k];
		} while (k != -1);
		cout << " }";
	}
	cout << ", size = " << part_size << endl;
//#endif SHOWINSTANCE
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
			IloNumColumn stable_set = Xobj(cost[part_set[current_part]]);
			/* fill the column corresponding to ">= 1" constraints (insert "1" in constraint indexed by v) */
			for (int v = 0; v < vertices; v++) if (setR[v]) stable_set += Xrestr[v](1.0);
			/* and the ">= -|C| constraint (insert "-1" in constraint indexed by color) */
			stable_set += Xrestr[vertices + current_part](-1.0);
#ifdef ONLYRELAXATION
			/* add the column as a non-negative continuos variable */
			Xvars.add(IloNumVar(stable_set));
#else
			/* add the column as a non-negative integer variable */
			Xvars.add(IloIntVar(stable_set));
#endif
			variables++;
			if (variables % 10000 == 0) {
				cout << " " << variables << " columns generated     \r";
				if (variables > MAXCOLUMNS) bye("Max. limit of columns reached!");
			}

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
 * optimize1 - make an exhaustive enmeration of stable sets and solve the set-cover formulation
 */
int optimize1()
{
	cout << "Using the set-cover formulation with an exhaustive enumeration of stable sets." << endl;
	setP = new bool[vertices];
	setR = new bool[vertices];
	setX = new bool[vertices];
	setV = new bool[vertices];

	Xobj = IloMinimize(Xenv);

	/* we will have "vertices" constraints with r.h.s >= 1 and "colors" constraints with r.h.s >= -|indistiguish.| */
	for (int v = 0; v < vertices; v++) Xrestr.add(IloRange(Xenv, 1.0, IloInfinity));
	for (int p = 0; p < part_size; p++) Xrestr.add(IloRange(Xenv, -part_card[p], IloInfinity));

	/* enumerate all the stable sets of the graph G[C(k)] with Bron–Kerbosch algorithm */
	variables = 0;
	for (int p = 0; p < part_size; p++) {
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
		int k = part_set[p];
		for (int s = 0; s < C_size[k]; s++) {
			int v = C_set[k][s];
			setP[v] = true; sizeP++;
			setV[v] = true; sizeV++;
		}
		current_part = p;
		enumerate_stable_sets();
	}
	cout << "Number of variables: " << variables << endl;

	Xmodel.add(Xobj);
	Xmodel.add(Xrestr);

	IloCplex cplex(Xmodel);
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
	cplex.setParam(IloCplex::BoolParam::MemoryEmphasis, CPX_ON); /* for memory saving */

#ifdef TUNEDPARAMS
	//cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Algorithm::Primal);
	//cplex.setParam(IloCplex::Param::NodeAlgorithm, IloCplex::Algorithm::Primal);
	//cplex.setParam(IloCplex::Param::Preprocessing::Presolve, CPX_OFF);
	//cplex.setParam(IloCplex::Param::Preprocessing::RepeatPresolve, 0);
	//cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 0);
	cplex.setParam(IloCplex::Param::MIP::Cuts::Cliques, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::ZeroHalfCut, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::BQP, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::Covers, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::Disjunctive, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::FlowCovers, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::Gomory, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::GUBCovers, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::Implied, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::LiftProj, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::LocalImplied, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::MCFCut, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::MIRCut, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::PathCut, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::RLT, -1);
	//cplex.setParam(IloCplex::Param::MIP::Strategy::Search, CPX_MIPSEARCH_TRADITIONAL);
	cplex.setParam(IloCplex::Param::MIP::Strategy::VariableSelect, CPX_VARSEL_MAXINFEAS);
	cplex.setParam(IloCplex::Param::MIP::Strategy::NodeSelect, CPX_NODESEL_DFS);
	cplex.setParam(IloCplex::Param::MIP::Strategy::Branch, CPX_BRDIR_UP);
	//cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);
	//cplex.setParam(IloCplex::Param::MIP::Strategy::Probe, -1);
	//cplex.setParam(IloCplex::Param::MIP::Strategy::FPHeur, -1);
	//cplex.setParam(IloCplex::Param::MIP::Strategy::RINSHeur, -1);
	//cplex.setParam(IloCplex::Param::MIP::Strategy::LBHeur, CPX_OFF);
#endif

	cplex.extract(Xmodel);
#ifdef SAVELP
	cplex.exportModel(SAVELP);
	cout << "Integer formulation saved" << endl;
#endif

	/* solve it! */
	nodes_explored = -1;
	lowerbound = -BIGNUMBER;
	upperbound = BIGNUMBER;
	cplex.solve();
	IloCplex::CplexStatus status = cplex.getCplexStatus();

	int optim_flag = 1;
#ifdef ONLYRELAXATION
	/* LP treatment */
	if (status != IloCplex::Optimal) {
		switch (status) {
		case IloCplex::InfOrUnbd:
		case IloCplex::Infeasible:  cout << "Infeasible :(" << endl; return 2;
		case IloCplex::AbortTimeLim: cout << "Time limit reached!" << endl; break;
		default: bye("Unexpected error :(");
		}
		optim_flag = 0;
	}
	else {
		/* optimality reached */
		set_color(10);
		cout << "LP relaxation solved! :)" << endl;
		lowerbound = cplex.getObjValue();
		cout << "  objective = " << lowerbound << endl;
		optim_flag = 3;
	}
#else
	/* MIP treatment */
	if (status != IloCplex::Optimal) {
		switch (status) {
		case IloCplex::InfOrUnbd:
		case IloCplex::Infeasible:  cout << "Infeasible :(" << endl; return 2;
		case IloCplex::AbortTimeLim: cout << "Time limit reached!" << endl; break;
		default: bye("Unexpected error :(");
		}
		optim_flag = 0;
	}

	/* read bounds */
	lowerbound = cplex.getBestObjValue();
	upperbound = cplex.getObjValue();

	if (lowerbound <= 0.0) lowerbound = -BIGNUMBER;
	if (upperbound >= BIGNUMBER) upperbound = BIGNUMBER;
	else {
		/* save optimal solution */
		for (int k = 0; k < colors; k++) colors_used[k] = false;
		int *available_color = new int[part_size];
		for (int p = 0; p < part_size; p++) available_color[p] = part_set[p];
		IloInt *Id_assigned = new IloInt[colors];
		int *color_assigned = new int[colors];
		int assigned = 0;
		for (int v = 0; v < vertices; v++) {
			int part_chosen = -1;
			IloInt Id_chosen = -1;
			IloRange expr1 = Xrestr[v];
			for (IloExpr::LinearIterator it1 = expr1.getLinearIterator(); it1.ok(); ++it1) {
				if (it1.getCoef() > 0.5) {
					IloNumVar var1 = it1.getVar();
					if (cplex.getValue(var1) > 0.5) {
						IloInt Id1 = var1.getId();
						for (int p = 0; p < part_size; p++) {
							IloRange expr2 = Xrestr[vertices + p];
							for (IloExpr::LinearIterator it2 = expr2.getLinearIterator(); it2.ok(); ++it2) {
								if (it2.getCoef() < -0.5) {
									IloNumVar var2 = it2.getVar();
									IloInt Id2 = var2.getId();
									if (Id1 == Id2) {
										part_chosen = p;
										Id_chosen = Id1;
										goto color_ch;
									}
								}
							}
						}
						bye("No color chosen!");
					}
				}
			}
			bye("No variable in 1 for current vertex!");
color_ch:;
			/* search for the color assigned to that variable, and assign it to v */
			int color_chosen = -1;
			for (int i = 0; i < assigned; i++) {
				if (Id_chosen == Id_assigned[i]) {
					color_chosen = color_assigned[i];
					break;
				}
			}
			/* in case the variable has no color previously assigned, assign one from the partition */
			if (color_chosen == -1) {
				color_chosen = available_color[part_chosen];
				if (color_chosen == -1) bye("Elements of partition exhausted!");
				available_color[part_chosen] = next_part[color_chosen];
				Id_assigned[assigned] = Id_chosen;
				color_assigned[assigned] = color_chosen;
				assigned++;
			}
			optimal_coloring[v] = color_chosen;
			colors_used[color_chosen] = true;
		}
		delete[] available_color;
	}

	if (lowerbound < upperbound) {
		set_color(12);
		cout << "Bounds: LB = " << lowerbound << ", UB = " << upperbound << "." << endl;
		if (0.0 < lowerbound && upperbound < BIGNUMBER) {
			/* Rel Gap = |best-integer - bestbound|/|bestinteger| */
			cout << "Relative gap = " << 100.0 * (upperbound - lowerbound) / upperbound << "." << endl;
		}
	}

	if (optim_flag == 1) {
		/* optimality reached */
		set_color(10);
		cout << "Optimality reached! :)" << endl;
		cout << "  objective = " << lowerbound << endl;
		nodes_explored = cplex.getNnodes();
		cout << "  nodes evaluated = " << nodes_explored << endl;
	}
#endif
	set_color(7);

	delete[] setV;
	delete[] setX;
	delete[] setR;
	delete[] setP;
	return optim_flag;
}

/*
 * find_max_clique - for a given adjacency matrix and a given S, find a maximal clique Q on G[S]
 * warning: S_set is modified
 */
void find_max_clique()
{
	/* initialize Q */
	for (int v = 0; v < vertices; v++) Q_set[v] = false;
	bool flag;
	do {
		/* find the vertex from S with largest degree */
		int bestv = -1;
		int bestdegree = -1;
		for (int v = 0; v < vertices; v++) {
			if (S_set[v]) {
				int degree = 0;
				for (int w = 0; w < vertices; w++) {
					if (v != w && S_set[w] && adj_mat[v][w]) degree++;
				}
				if (degree > bestdegree) {
					bestv = v;
					bestdegree = degree;
				}
			}
		}
		if (bestv == -1) bye("Internal error!");
		/* add the vertex to the clique and delete it and its non-neighbors from S */
		Q_set[bestv] = true;
		S_set[bestv] = false;
		for (int w = 0; w < vertices; w++) {
			if (w != bestv && S_set[w] && adj_mat[bestv][w] == false) S_set[w] = false;
		}
		/* if there are still vertices in S, repeat the procedure */
		flag = false;
		for (int v = 0; v < vertices; v++) {
			if (S_set[v]) {
				flag = true;
				break;
			}
		}
	} while (flag);
	/* show the clique found */
//	cout << "Clique found: {";
//	for (int v = 0; v < vertices; v++) if (Q_set[v]) cout << " " << v;
//	cout << " }" << endl;
}

/*
* optimize2 - use the standard vertex-color formulation
*/
int optimize2()
{
	cout << "Using the standard vertex-color formulation." << endl;

	variables = 0;

	/* create a binary variable per color */
	int *w_index = new int[colors];
	for (int k = 0; k < colors; k++) {
#ifdef ONLYRELAXATION
		Xvars.add(IloNumVar(Xenv, 0.0, 1.0, ILOFLOAT));
#else
		Xvars.add(IloNumVar(Xenv, 0.0, 1.0, ILOBOOL));
#endif
		w_index[k] = variables;
		variables++;
	}

	/* create a binary variable per vertex and color */
	int **x_index = new int*[vertices];
	for (int v = 0; v < vertices; v++) {
		x_index[v] = new int[colors];
		for (int k = 0; k < colors; k++) x_index[v][k] = -1;
		for (int c = 0; c < L_size[v]; c++) {
#ifdef ONLYRELAXATION
			Xvars.add(IloNumVar(Xenv, 0.0, 1.0, ILOFLOAT));
#else
			Xvars.add(IloNumVar(Xenv, 0.0, 1.0, ILOBOOL));
#endif
			int k = L_set[v][c];
			x_index[v][k] = variables;
			variables++;
		}
	}
	cout << "Number of variables: " << variables << endl;

	/* generate objective function */
	IloExpr fobj(Xenv, 0);
	for (int k = 0; k < colors; k++) fobj += cost[k] * Xvars[w_index[k]];
	Xmodel.add(IloMinimize(Xenv, fobj));

	/* generate assignment constraints */
	int count = 0;
	for (int v = 0; v < vertices; v++) {
		IloExpr restr(Xenv);
		for (int c = 0; c < L_size[v]; c++) {
			int k = L_set[v][c];
			restr += Xvars[x_index[v][k]];
		}
		Xmodel.add(restr == 1);
		count++;
	}
	cout << " assingment constraints: " << count << endl;

	/* we need to know which colors are not covered by the adjacency constraints */
	bool **C_covered = new bool*[colors];
	for (int k = 0; k < colors; k++) {
		C_covered[k] = new bool[C_size[k]];
		for (int s = 0; s < C_size[k]; s++) C_covered[k][s] = false;
	}

#ifdef EDGEINEQ
	/* generate adjacency constraints */
	count = 0;
	for (int e = 0; e < edges; e++) {
		int u = edge_u[e];
		int v = edge_v[e];
		for (int k = 0; k < colors; k++) {
			int iuk = x_index[u][k];
			int ivk = x_index[v][k];
			if (iuk != -1 && ivk != -1) {
				Xmodel.add(Xvars[iuk] + Xvars[ivk] - Xvars[w_index[k]] <= 0);
				count++;
				for (int s = 0; s < C_size[k]; s++) { /* also take into account that u and v cover k */
					if (C_set[k][s] == u || C_set[k][s] == v) C_covered[k][s] = true;
				}
			}
		}
	}
	cout << " edge constraints: " << count << endl;
#else
	/* generate a cover of edges by cliques */
	adj_mat = new bool*[vertices];
	for (int v = 0; v < vertices; v++) adj_mat[v] = new bool[vertices];
	S_set = new bool[vertices];
	Q_set = new bool[vertices];
	count = 0;
	for (int k = 0; k < colors; k++) {
		/* let k be a color used by at least two vertices */
		if (C_size[k] < 2) continue;
		/* construct the adjacency matrix of Gk */
		for (int u = 0; u < vertices; u++) {
			for (int v = 0; v < vertices; v++) adj_mat[u][v] = false;
		}
		for (int s1 = 0; s1 < C_size[k] - 1; s1++) {
			int u = C_set[k][s1];
			for (int s2 = s1 + 1; s2 < C_size[k]; s2++) {
				int v = C_set[k][s2];
				if (adjacency[u][v] > 0) {
					adj_mat[u][v] = true;
					adj_mat[v][u] = true;
				}
			}
		}
		bool flag;
		do {
			/* construct the set such that each vertex has at least one non-covered edge */
			for (int v = 0; v < vertices; v++) S_set[v] = false;
			flag = false;
			for (int s = 0; s < C_size[k]; s++) {
				int v = C_set[k][s];
				for (int w = 0; w < vertices; w++) {
					if (v != w && adj_mat[v][w]) {
						S_set[v] = true;
						flag = true;
						break;
					}
				}
			}
			if (flag) {
				/* find a maximal clique */
				//cout << "Color " << k << "." << endl;
				find_max_clique();
				/* try to maximalize it (if it is not maximal yet) */
				for (int s = 0; s < C_size[k]; s++) {
					int v = C_set[k][s];
					if (Q_set[v] == false) {
						/* is v adjacent to Q? */
						bool vadjQ = true;
						for (int t = 0; t < C_size[k]; t++) {
							int q = C_set[k][t];
							if (Q_set[q] && adjacency[v][q] == 0) {
								vadjQ = false;
								break;
							}
						}
						if (vadjQ) {
							Q_set[v] = true;
							//cout << " + " << v << endl;
						}
					}
				}
				IloExpr restr(Xenv);
				for (int s = 0; s < C_size[k]; s++) {
					int v = C_set[k][s];
					if (Q_set[v]) {
						restr += Xvars[x_index[v][k]];
						C_covered[k][s] = true;
					}
				}
				restr -= Xvars[w_index[k]];
				Xmodel.add(restr <= 0);
				count++;
				/* remove the adjacencies covered by that clique */
				for (int s1 = 0; s1 < C_size[k] - 1; s1++) {
					int u = C_set[k][s1];
					if (Q_set[u]) {
						for (int s2 = s1 + 1; s2 < C_size[k]; s2++) {
							int v = C_set[k][s2];
							if (Q_set[v]) {
								adj_mat[u][v] = false;
								adj_mat[v][u] = false;
							}
						}
					}
				}
			}
		} while (flag);
	}
	delete[] Q_set;
	delete[] S_set;
	for (int v = 0; v < vertices; v++) delete[] adj_mat[v];
	delete[] adj_mat;
	cout << " clique constraints: " << count << endl;
#endif

	/* generate remaining constraints that cover colors */
	count = 0;
	for (int k = 0; k < colors; k++) {
		for (int s = 0; s < C_size[k]; s++) {
			if (C_covered[k][s] == false) {
				int v = C_set[k][s];
				Xmodel.add(Xvars[x_index[v][k]] - Xvars[w_index[k]] <= 0);
				count++;
			}
		}
	}
	cout << " isolated vertices constraints: " << count << endl;

	/* symmetry-breaking constraints: */
#ifdef SYMMETRYRESTR1
	/* w(k) >= w(k+1),     k in partition = {0, 1, 2, ...}
	   w(k) = 0,           for those k greater than |Vk| */
	count = 0;
	for (int i = 0; i < part_size; i++) {
		if (part_card[i] >= 2) {
			int k = part_set[i];
			int Vsize = C_size[k];
			int knext = next_part[k];
			int kcount = 0;
			do {
				if (kcount < Vsize-1) {
					Xmodel.add(Xvars[w_index[k]] - Xvars[w_index[knext]] >= 0);
					//cout << " w(" << k << ") - w(" << knext << ") >= 0" << endl;
				}
				else {
					Xmodel.add(Xvars[w_index[knext]] == 0);
					//cout << " w(" << knext << ") = 0" << endl;
				}
				count++;
				kcount++;
				k = knext;
				knext = next_part[k];
			} while (knext != -1);
		}
	}
	cout << " symmetry-breaking Type 1 constraints: " << count << endl;
#endif
#ifdef SYMMETRYRESTR2
	/* sum_{v in Vk} x(v,k) >= sum_{v in Vk} x(v,k+1),     k in partition = {0, 1, 2, ...} */
	count = 0;
	for (int i = 0; i < part_size; i++) {
		if (part_card[i] >= 2) {
			int first = part_set[i];
			int Vsize = C_size[first];
			int k = first;
			int knext = next_part[k];
			do {
				IloExpr restr(Xenv);
				for (int c = 0; c < Vsize; c++) {
					int v = C_set[first][c];
					restr += Xvars[x_index[v][k]];
					if (x_index[v][knext] == -1) bye("Internal error!");
					restr -= Xvars[x_index[v][knext]];
				}
				Xmodel.add(restr >= 0);
				count++;
				k = knext;
				knext = next_part[k];
			} while (knext != -1);
		}
	}
	cout << " symmetry-breaking Type 2 constraints: " << count << endl;
#endif
#ifdef SYMMETRYRESTR3
	/* x(v,k) = 0,    v in Vk,  v < k  (k in partition, v ordered in descending order from Gk) */
	count = 0;
	for (int i = 0; i < part_size; i++) {
		if (part_card[i] >= 2) {
			int first = part_set[i];
			int Vsize = C_size[first];
			/* compute the (modified) degree of each vertex of Gk */
			int *degree1 = new int[Vsize];
			for (int c = 0; c < Vsize; c++) {
				int v = C_set[first][c];
				int deg = 2;
				for (int d = 0; d < Vsize; d++) {
					if (c != d) {
						int w = C_set[first][d];
						if (adjacency[v][w] > 0) deg++;
					}
				}
				degree1[c] = deg;
			}
			/* compute the 2nd-order degree of each vertex */
			struct sortstr *degree2 = new struct sortstr[Vsize];
			for (int c = 0; c < Vsize; c++) {
				int v = C_set[first][c];
				int deg = 0;
				for (int d = 0; d < Vsize; d++) {
					if (c != d) {
						int w = C_set[first][d];
						if (adjacency[v][w] > 0) deg += degree1[d];
					}
				}
				degree2[c].i = c;
				degree2[c].val = deg;
			}
			/* sort vertices according to their 2nd-order degree */
			qsort(degree2, Vsize, sizeof(sortstr), cmpfunc);
			/*cout << "Order:";
			for (int j = 0; j < Vsize; j++) {
				int c = degree2[j].i;
				int val = degree2[j].val;
				int v = C_set[first][c];
				cout << " (" << v << "," << val << ")";
			}
			cout << endl;*/
			for (int j = 0; j < Vsize; j++) {
				int v = C_set[first][degree2[j].i];
				int k = first;
				for (int l = 0; l <= j; l++) {
					k = next_part[k];
					if (k == -1) goto skip_k;
				}
				do {
					Xmodel.add(Xvars[x_index[v][k]] == 0);
					count++;
					//cout << " x(" << v << "," << k << ") = 0" << endl;
					k = next_part[k];
				} while (k != -1);
skip_k:;
			}
			delete[] degree2;
			delete[] degree1;
		}
	}
	cout << " symmetry-breaking Type 3 constraints: " << count << endl;
#endif

	IloCplex cplex(Xmodel);
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
	cplex.setParam(IloCplex::BoolParam::MemoryEmphasis, CPX_ON); /* for memory saving */

#ifdef TUNEDPARAMS
	cplex.setParam(IloCplex::Param::MIP::Cuts::Cliques, 1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::ZeroHalfCut, 1);

	cplex.setParam(IloCplex::Param::MIP::Cuts::BQP, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::Covers, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::Disjunctive, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::FlowCovers, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::Gomory, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::GUBCovers, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::Implied, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::LiftProj, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::LocalImplied, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::MCFCut, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::MIRCut, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::PathCut, -1);
	cplex.setParam(IloCplex::Param::MIP::Cuts::RLT, -1);

	cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, 100);
	cplex.setParam(IloCplex::Param::MIP::Strategy::Backtrack, 0.1);
#endif

	cplex.extract(Xmodel);
#ifdef SAVELP
	cplex.exportModel(SAVELP);
	cout << "Integer formulation saved" << endl;
#endif

	/* solve it! */
	nodes_explored = -1;
	lowerbound = -BIGNUMBER;
	upperbound = BIGNUMBER;
	cplex.solve();
	IloCplex::CplexStatus status = cplex.getCplexStatus();

	int optim_flag = 1;
#ifdef ONLYRELAXATION
	/* LP treatment */
	if (status != IloCplex::Optimal) {
		switch (status) {
		case IloCplex::InfOrUnbd:
		case IloCplex::Infeasible: cout << "Infeasible :(" << endl; return 2;
		case IloCplex::AbortTimeLim: cout << "Time limit reached!" << endl; break;
		default: bye("Unexpected error :(");
		}
		optim_flag = 0;
	}
	else {
		/* optimality reached */
		set_color(10);
		cout << "LP relaxation solved! :)" << endl;
		lowerbound = cplex.getObjValue();
		cout << "  objective = " << lowerbound << endl;
		optim_flag = 3;
	}
#else
	/* MIP treatment */
	if (status != IloCplex::Optimal) {
		switch (status) {
		case IloCplex::InfOrUnbd:
		case IloCplex::Infeasible: cout << "Infeasible :(" << endl; return 2;
		case IloCplex::AbortTimeLim: cout << "Time limit reached!" << endl; break;
		default: bye("Unexpected error :(");
		}
		optim_flag = 0;
	}

	/* read bounds */
	lowerbound = cplex.getBestObjValue();
	upperbound = cplex.getObjValue();

	if (lowerbound <= 0.0) lowerbound = -BIGNUMBER;
	if (upperbound >= BIGNUMBER) upperbound = BIGNUMBER;
	else {
		/* save optimal solution */
		for (int k = 0; k < colors; k++) colors_used[k] = false;
		for (int v = 0; v < vertices; v++) {
			for (int k = 0; k < colors; k++) {
				int ivk = x_index[v][k];
				if (ivk != -1) {
					if (cplex.getValue(Xvars[ivk]) > 0.5) {
						optimal_coloring[v] = k;
						colors_used[k] = true;
					}
				}
			}
		}
	}

	if (lowerbound < upperbound) {
		set_color(12);
		cout << "Bounds: LB = " << lowerbound << ", UB = " << upperbound << "." << endl;
		if (0.0 < lowerbound && upperbound < BIGNUMBER) {
			/* Rel Gap = |best-integer - bestbound|/|bestinteger| */
			cout << "Relative gap = " << 100.0 * (upperbound - lowerbound) / upperbound << "." << endl;
		}
	}

	if (optim_flag == 1) {
		/* optimality reached */
		set_color(10);
		cout << "Optimality reached! :)" << endl;
		cout << "  objective = " << lowerbound << endl;
		nodes_explored = cplex.getNnodes();
		cout << "  nodes evaluated = " << nodes_explored << endl;
	}
#endif
	set_color(7);

	/* free memory */
	for (int k = 0; k < colors; k++) delete[] C_covered[k];
	delete[] C_covered;
	for (int v = 0; v < vertices; v++) delete[] x_index[v];
	delete[] x_index;
	delete[] w_index;
	return optim_flag;
}

/*
 * write_sol - Write output file
 */
void write_sol(char *filename, int optim_flag, double time_elapsed)
{
	/* open file */
	FILE *stream = fopen(filename, "wt");
	if (!stream) bye("Output file cannot be opened");
#ifdef ONLYRELAXATION
	fprintf(stream, "%d:%.2f:%.2f:-1:%.2f:0\n", optim_flag, lowerbound, upperbound, time_elapsed);
#else
	if (optim_flag == 2) fprintf(stream, "2:%.2f:%.2f:-1:%.2f:0\n", lowerbound, upperbound, time_elapsed);
	else {
		int col_used = 0;
		if (upperbound < BIGNUMBER) {
			for (int k = 0; k < colors; k++) if (colors_used[k]) col_used++;
		}
		if (optim_flag == 1) fprintf(stream, "1:%.2f:%.2f:%ld:%.2f:%d\n", lowerbound, upperbound, nodes_explored, time_elapsed, col_used);
		else fprintf(stream, "0:%.2f:%.2f:-1:%.2f:%d\n", lowerbound, upperbound, MAXTIME, col_used);
		if (upperbound < BIGNUMBER) {
			for (int k = 0; k < colors; k++) if (colors_used[k]) fprintf(stream, "%d\n", k);
			for (int v = 0; v < vertices; v++) fprintf(stream, "%d\n", optimal_coloring[v]);
		}
	}
#endif
	fclose(stream);
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
	cout << "LISTCOLA - Solves the Minimum Cost List Coloring Problem (w/Alt models)." << endl;
	set_color(7);

	if (argc < 2) {
		cout << "Usage: listcola file" << endl;
		cout << "  Read file.graph, file.cost and file.list and solve the instance" << endl;
		cout << "  It generates file.out with solution and statistics formatted as follows:" << endl;
		cout << "     optimality_flag:LB_total_weight:UB_total_weight:nodes_explored:time_elapsed:col_used" << endl;
		cout << "     first color used" << endl;
		cout << "     second color used" << endl;
		cout << "       ..." << endl;
		cout << "     last color used" << endl;
		cout << "     color of first vertex" << endl;
		cout << "     color of second vertex" << endl;
		cout << "       ..." << endl;
		cout << "     color of last vertex" << endl;
		bye("Bye!");
	}

	char *filename = argv[1];
	char filename_extension[210];

	/* read instance */
	double start_t = ECOclock();
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".graph");
	read_graph(filename_extension);
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".cost");
	read_cost(filename_extension);
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".list");
	read_list(filename_extension);
	double stop_t = ECOclock();
	set_color(15);
	cout << "Time of instance reading = " << stop_t - start_t << " sec." << endl;
	set_color(7);

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
		for (int s = 0; s < L_size[v]; s++) cout << " " << L_set[v][s];
		cout << " }" << endl;
	}
	for (int k = 0; k < colors; k++) {
		cout << "V(G" << k << ") = {";
		for (int s = 0; s < C_size[k]; s++) cout << " " << C_set[k][s];
		cout << " }" << endl;
	}
#endif

	/* Show some basic statistics */
	set_color(6);
	cout << "Statistics:" << endl;
	int clique_size = vertices * (vertices - 1) / 2;
	float density = 100.0 * (float)edges / (float)clique_size;
	cout << "  |V| = " << vertices << ", |E| = " << edges << " (density = " << density << "%), |C| = " << colors << "." << endl;
#ifdef SHOWSTATS
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
	/* Average density of V(Gk) */
	prom = 0.0;
	int count = 0;
	for (int k = 0; k < colors; k++) {
		int count_vertices = C_size[k];
		if (count_vertices > 1) {
			int count_edges = 0;
			for (int s1 = 0; s1 < count_vertices - 1; s1++) {
				int u = C_set[k][s1];
				for (int s2 = s1 + 1; s2 < count_vertices; s2++) {
					int v = C_set[k][s2];
					if (adjacency[u][v] > 0) count_edges++;
				}
			}
			prom += (float)count_edges / (float)(count_vertices * (count_vertices - 1) / 2);
			count++;
		}
	}
	prom /= (float)count;
	cout << "      average density of Gk = " << 100.0 * prom << "%" << endl;
	/* Superposition between Gk's: average of |V(Gk1) cap V(Gk2)| / |V(Gk1) cup V(Gk2)| */
	prom = 0.0;
	count = 0;
	for (int k1 = 0; k1 < colors - 1; k1++) {
		for (int k2 = k1 + 1; k2 < colors; k2++) {
			int count_intersection = 0;
			int count_union = 0;
			for (int v = 0; v < vertices; v++) {
				bool in_k1 = false;
				for (int s = 0; s < C_size[k1]; s++) {
					int u = C_set[k1][s];
					if (u == v) {
						in_k1 = true;
						break;
					}
				}
				bool in_k2 = false;
				for (int s = 0; s < C_size[k2]; s++) {
					int u = C_set[k2][s];
					if (u == v) {
						in_k2 = true;
						break;
					}
				}
				if (in_k1 || in_k2) count_union++;
				if (in_k1 && in_k2) count_intersection++;
			}
			prom += (float)count_intersection / (float)count_union;
			count++;
		}
	}
	prom /= (float)count;
	cout << " superposition between Gk's = " << 100.0 * prom << "%" << endl;
#endif

	/* optimize using an exhaustive enmeration of stable sets */
	set_color(7);
	header_part = new int[colors];
	next_part = new int[colors];
	part_set = new int[colors];
	part_card = new int[colors];
	part_size = 0;
	optimal_coloring = new int[vertices];
	colors_used = new bool[colors];
	start_t = ECOclock();
	find_color_part();
#ifdef  STABLEMODEL
	int optim_flag = optimize1();
#else
	int optim_flag = optimize2();
#endif
	stop_t = ECOclock();
	double time_elapsed = stop_t - start_t;
#ifndef ONLYRELAXATION
	if (optim_flag == 1) {
		/* show the answer on screen */
		set_color(3);
		cout << "Coloring:" << endl;
		for (int v = 0; v < vertices; v++) cout << "  f(" << v << ") = " << optimal_coloring[v] << endl;
		cout << "Colors used:";
		for (int k = 0; k < colors; k++) if (colors_used[k]) cout << " " << k;
		cout << endl;
	}
#endif
	set_color(15);
	cout << "Time of optimization = " << time_elapsed << " sec." << endl;
	set_color(7);

	/* save solution and stats on output file */
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".out");
	write_sol(filename_extension, optim_flag, time_elapsed);

	/* free memory */
	delete[] colors_used;
	delete[] optimal_coloring;
	delete[] part_card;
	delete[] part_set;
	delete[] next_part;
	delete[] header_part;
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
