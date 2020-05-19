/*
 * GENMUINST - Generate a random instance of the mu-coloring Problem from a graph
 * Made in 2018-2020 by Daniel Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>

#ifdef VISUALC
#include <ctime>
#endif

ILOSTLBEGIN

using namespace std;

#define INFDIST 9999999
#define EPSILON 0.00001
#define MAXTIME 600.0

/* FUNCTIONS */

/*
 * bye - finish executing and show a message
 */
void bye(char *string)
{
	cout << string << endl;
	exit(-1);
}

/*
 * urnd - generate numbers with uniform random distribution in the interval is [a, b)
 */
float urnd(float a, float b)
{
	return a + rand() * (b - a) / ((float)RAND_MAX + 1.0);
}

/*
 * main - Main program
 */
int main(int argc, char **argv)
{
	cout << "GENMUINST - Generate a random instance of mu-coloring Problem." << endl;

	if (argc < 5) {
		cout << "Usage: genmuinst file dis_colors max_repeat tprob" << endl;
		cout << "It takes file.graph and generates file.cost and file.list." << endl;
		cout << "Parameters:" << endl;
		cout << "  dis_colors = number of distinguishable colors (card. of K)" << endl;
		cout << "  max_repeat = max times the color is repeated (uniform distribution between 1 and max_repeat)" << endl;
		cout << "  tprob = prob. (1%-100%) of removing a vertex from a color to the next one" << endl;
		bye("Bye!");
	}

	char *filename = argv[1];
	int dis_colors = atoi(argv[2]);
	if (dis_colors < 3 || dis_colors > 10000) bye("Number of colors out range!");
	int max_repeat = atoi(argv[3]);
	if (max_repeat < 0 || max_repeat > 100) bye("max_repeat param. out of range!");
	int tprob = atoi(argv[4]);
	if (tprob < 1 || tprob > 100) bye("tprob param. out of range!");

	/* open graph file */
	char filename_extension[210];
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".graph");
	FILE *stream = fopen(filename_extension, "rt");
	if (!stream) bye("Graph file cannot be opened");
	int vertices, edges;
	fscanf(stream, "%d:%d\n", &vertices, &edges);

	/* do not accept graph of less than 4 vertices or stable sets */
	if (vertices < 4) bye("Number of vertices out range!");
	if (edges < 1 || edges > vertices*(vertices - 1) / 2) bye("Number of edges out of range!");

	/* ask for memory */
	bool **adjacency = new bool*[vertices];
	for (int u = 0; u < vertices; u++) {
		adjacency[u] = new bool[vertices];
		for (int v = 0; v < vertices; v++) adjacency[u][v] = false;
	}

	/* read edges */
	for (int e = 0; e < edges; e++) {
		int u, v;
		fscanf(stream, "%d,%d\n", &u, &v);
		if (u < 0 || u >= v || v >= vertices) {
			cout << "Error reading edge " << e + 1 << "!" << endl;
			bye("Bye!");
		}
		if (adjacency[u][v]) {
			cout << "A repeated edge was found: (" << u << ", " << v << ")" << endl;
			bye("Bye!");
		}
		else {
			adjacency[u][v] = true;
			adjacency[v][u] = true;
		}
	}
	fclose(stream);

	/* compute a clique of maximum size */
	IloEnv Xenv;
	IloNumVarArray Xvars(Xenv);
	for (int v = 0; v < vertices; v++) Xvars.add(IloNumVar(Xenv, 0.0, 1.0, ILOBOOL));

	/* generate objective function */
	IloModel Xmodel(Xenv);
	IloExpr fobj(Xenv, 0);
	for (int v = 0; v < vertices; v++) fobj += 1.0 * Xvars[v];
	Xmodel.add(IloMaximize(Xenv, fobj));

	/* generate anti-adjacency constraints */
	for (int u = 0; u < vertices - 1; u++) {
		for (int v = u + 1; v < vertices; v++) {
			if (adjacency[u][v] == false) Xmodel.add(Xvars[u] + Xvars[v] <= 1);
		}
	}

	IloCplex cplex(Xmodel);
	cplex.setDefaults();
	cplex.setOut(Xenv.getNullStream());
	cplex.setWarning(Xenv.getNullStream());
	cplex.setParam(IloCplex::NumParam::WorkMem, 2048);
	cplex.setParam(IloCplex::NumParam::TreLim, 2048);
	cplex.setParam(IloCplex::IntParam::NodeFileInd, 0);
	cplex.setParam(IloCplex::NumParam::TiLim, MAXTIME);
	cplex.setParam(IloCplex::NumParam::EpGap, 0.0);
	cplex.setParam(IloCplex::NumParam::EpAGap, 0.0);
	cplex.setParam(IloCplex::NumParam::EpInt, EPSILON);
	cplex.setParam(IloCplex::IntParam::Threads, 1);
	cplex.setParam(IloCplex::IntParam::RandomSeed, 1);

	cplex.extract(Xmodel);
	//cplex.exportModel("form.lp");

	/* solve it! */
	cplex.solve();
	IloCplex::CplexStatus status = cplex.getCplexStatus();
	if (status != IloCplex::Optimal) {
		switch (status) {
		case IloCplex::InfOrUnbd:
		case IloCplex::Infeasible: bye("Infeasible :(");
		case IloCplex::AbortTimeLim: bye("Time limit reached :(");
		default: bye("Unexpected error :(");
		}
	}

	/* show the clique found */
	int max_clique_size = 0;
	cout << "Max Clique = {";
	for (int v = 0; v < vertices; v++) {
		if (cplex.getValue(Xvars[v]) > 0.5) {
			cout << " " << v;
			max_clique_size++;
		}
	}
	cout << " }, size = " << max_clique_size << endl;

	/* generate multiplicity of colors with "max_repeat" parameter:
	     repeated[c] = false   if c is not a repeated color
	     repeated[c] = true    if c is color q, where repeated[q] = false, repeated[q+1] = ... = repeated[c-1] = true
	   Extreme case:
	     if max_repeat = 1 then repeated[k] = false for all k (all colors are different)
	Exceptions:
	     multiplicity of first color is clique_size to avoid trivial infeasible instances */
	srand(time(0));
	int colors = 0;
	bool *repeated = new bool[dis_colors*max_repeat];
	for (int k = 0; k < dis_colors; k++) {
		repeated[colors++] = false;
		int multiplicity = max_clique_size - 1;
		if (k > 0) multiplicity = (int)urnd(0.0, (float)max_repeat);
		for (int m = 0; m < multiplicity; m++) repeated[colors++] = true;
	}

	/* show some stats */
	cout << "Statistics:" << endl;
	int clique_size = vertices * (vertices - 1) / 2;
	float density = 100.0 * (float)edges / (float)clique_size;
	cout << "  |V| = " << vertices << ", |E| = " << edges << " (density = " << density << "%)" << endl;
	cout << "  |K| = " << dis_colors << ", |C| = " << colors << endl;
	cout << "  max_repeat = " << max_repeat << ", t prob. = " << tprob << "%" << endl;

	/* write cost file */
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".cost");
	stream = fopen(filename_extension, "wt");
	if (!stream) bye("Cost file cannot be created");
	fprintf(stream, "%d\n", colors);
	for (int k = 0; k < colors; k++) fprintf(stream, "1 ");
	fprintf(stream, "\n");
	fclose(stream);

	/* generate mu's */
	int *mu = new int[vertices];
	int max_colors_reached = -1;
	for (int v = 0; v < vertices; v++) {
		if (v == vertices - 1 && max_colors_reached < colors) {
			/* last vertex uses every color if there are some colors unused */
			mu[v] = colors;
			continue;
		}
		int c = colors;
		for (int k = 1; k < colors; k++) {
			if (repeated[k] == false && urnd(0.0, 100.0) <= (float)tprob) {
				c = k;
				break;
			}
		}
		if (c > max_colors_reached) max_colors_reached = c;
		mu[v] = c;
	}

	/* write list file */
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".list");
	stream = fopen(filename_extension, "wt");
	if (!stream) bye("List file cannot be created");
	fprintf(stream, "%d:%d\n", vertices, colors);
	for (int v = 0; v < vertices; v++) {
		fprintf(stream, "%d  ", mu[v]);
		for (int s = 0; s < mu[v]; s++) fprintf(stream, "%d ", s);
		fprintf(stream, "\n");
	}
	fclose(stream);

	delete[] mu;
	delete[] repeated;
	for (int v = 0; v < vertices; v++) delete[] adjacency[v];
	delete[] adjacency;
	return 0;
}
