/*
 * CHECKER - Checks output solution
 * Made in 2018-2019 by Daniel Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

using namespace std;

#define EPSILON 0.00001

/* GLOBAL VARIABLES */

int vertices, edges; /* number of vertices and edges */
int *edge_u, *edge_v; /* array of endpoints of edges */
int **adjacency; /* adjacency matrix: 0 means no adjacency; >0 gives the index to the edge + 1 */

int colors; /* number of colors */
int *cost; /* cost of each color */
int *L_size; /* size of each set L(v), per vertex */
int **L_set; /* elements of L(v), per vertex */

/* FUNCTIONS */

/*
 * bye - finish executing and show a message
 */
void bye(char *string)
{
	cout << string << endl;
	exit(1);
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
	adjacency = new int*[vertices];
	for (int u = 0; u < vertices; u++) {
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
			edge_u[e] = u;
			edge_v[e] = v;
			adjacency[u][v] = e + 1;
			adjacency[v][u] = e + 1;
		}
	}
	fclose(stream);
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
	L_set = new int*[vertices];

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
		}
	}
	fclose(stream);
}

/*
 * check_coloring - check if the given coloring is valid, uses the correct colors and has the correct cost
 */
bool check_coloring(int *f, bool *colors_used, double obj) {
	/* colors of each vertex belong to the list and are from colors used ??? */
	for (int v = 0; v < vertices; v++) {
		int k = f[v];
		bool happy = false;
		for (int s = 0; s < L_size[v]; s++) {
			if (L_set[v][s] == k) {
				happy = true;
				break;
			}
		}
		if (happy == false) {
			cout << "Ouch! Vertex " << v << " is colored with " << k << " which does not belong to its list :S" << endl;
			return false;
		}
		if (colors_used[k] == false) {
			cout << "Ouch! Vertex " << v << " is colored with " << k << " which is not used :S" << endl;
			return false;
		}
	}

	/* are there conflicting edges ??? */
	for (int e = 0; e < edges; e++) {
		int u = edge_u[e];
		int v = edge_v[e];
		if (f[u] == f[v]) {
			cout << "Ouch! Vertices " << u << " and " << v << " have the same color " << f[v] << " :S" << endl;
			return false;
		}
	}

	/* has the coloring the proper total cost ??? */
	double akku = 0.0;
	for (int k = 0; k < colors; k++) if (colors_used[k]) akku += cost[k];
	if (akku < obj - EPSILON || akku > obj + EPSILON) {
		cout << "Ouch! Cost of solution is " << akku << " but it is reported as " << obj << " :S" << endl;
		return false;
	}
	return true;
}

/*
 * main - Main program
 */
int main(int argc, char **argv)
{
	cout << "CHECKER - Checks output solution." << endl;

	if (argc < 2) {
		cout << "Usage: checker file" << endl;
		cout << "  Read file.graph, file.cost, file.list and file.out" << endl;
		bye("Bye!");
	}

	char *filename = argv[1];
	char filename_extension[210];

	/* read instance */
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".graph");
	read_graph(filename_extension);
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".cost");
	read_cost(filename_extension);
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".list");
	read_list(filename_extension);

	/* read output file */
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".out");
	FILE *stream = fopen(filename_extension, "rt");
	if (!stream) bye("Output file cannot be opened");

	int optim_flag, col_used;
	double lowerbound, upperbound, time_elapsed;
	long nodes_explored;
	fscanf(stream, "%d:%lf:%lf:%ld:%lf:%d\n", &optim_flag, &lowerbound, &upperbound, &nodes_explored, &time_elapsed, &col_used);
	if (optim_flag < 0 || optim_flag > 2) bye("Optimality flag out of range!");
	if (optim_flag == 1 && (lowerbound < upperbound - EPSILON || lowerbound > upperbound + EPSILON)) {
		bye("Reported as optimal but lowerbound != upperbound!");
	}
	if (optim_flag == 0 && lowerbound >= upperbound) bye("Reported as non-optimal but lowerbound >= upperbound!");
	if (col_used < 0 || col_used > colors) bye("Number of colors used out of range!");
	if (col_used > 0) {
		/* check coloring */
		bool *colors_used = new bool[colors];
		for (int k = 0; k < colors; k++) colors_used[k] = false;
		for (int i = 0; i < col_used; i++) {
			int k;
			fscanf(stream, "%d", &k);
			if (k < 0 || k >= colors) bye("Color out of range (during read of list of colors used)!");
			colors_used[k] = true;
		}

		int *f = new int[vertices];
		for (int v = 0; v < vertices; v++) {
			int k;
			fscanf(stream, "%d", &k);
			if (k < 0 || k >= colors) bye("Color out of range (during read of coloring)!");
			f[v] = k;
		}

		if (check_coloring(f, colors_used, upperbound) == false) bye("Error during checking!");
		delete[] f;
		delete[] colors_used;
	}
	fclose(stream);

	/* check the solution file */
	if (optim_flag == 1) {
		strncpy(filename_extension, filename, 200);
		strcat(filename_extension, ".sol");
		FILE *stream = fopen(filename_extension, "rt");
		if (!stream) {
			/* if the file does not exist, generate a new one */
			FILE *stream = fopen(filename_extension, "wt");
			if (!stream) bye("Solution file cannot be opened");
			fprintf(stream, "%lf\n", upperbound);
			fclose(stream);
		}
		else {
			/* if the file exists, check if the optimal objective value is correct */
			double value;
			fscanf(stream, "%lf\n", &value);
			if (upperbound < value - EPSILON || upperbound > value + EPSILON) {
				printf("Solution value is %lf but reported one is %lf\n", value, upperbound);
				bye("Bye!");
			}
			fclose(stream);
		}
	}
	cout << "Ok!" << endl;

	/* free memory */
	for (int v = 0; v < vertices; v++) delete[] L_set[v];
	delete[] L_set;
	delete[] L_size;
	delete[] cost;
	delete[] edge_v;
	delete[] edge_u;
	for (int v = 0; v < vertices; v++) delete[] adjacency[v];
	delete[] adjacency;
	return 0;
}
