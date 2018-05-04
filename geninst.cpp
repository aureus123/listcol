/*
 * GENINST - Generate instances of Minimum Cost List Coloring Problem
 * Made in 2018 by Daniel Severin and Mauro Lucci
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace std;

/* CONSTANTS */

#define PI 3.14159265358979323846
#define COEF_DENSITY 38.8/100.0
#define EDGES_PER_COLOR 3.14
#define LOG_COLORS_PER_AVERAGE 0.785
#define AVERAGES_PER_SIGMA 1.5
#define COLORS_PER_LIST_SIZE 8


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
* urnd - generate numbers with uniform random distribution
*  if flag=false, the interval is [a, b]
*  if flag=true, the interval is [a, b)
*/
float urnd(float a, float b, bool flag)
{
	return a + rand() * (b - a) / (float)(RAND_MAX + (flag ? 1 : 0));
}

/*
* grnd - generate numbers with gaussian random distribution
* (it uses Box-Muller algorithm)
*/
float grnd(float mu, float sigma)
{
	float z, u1, u2;

	/* u1 and u2 are given from an uniform distribution of (0, 1] */
	u1 = 1.0 - urnd(0, 1, true);
	u2 = 1.0 - urnd(0, 1, true);

	if (urnd(0, 1, 1) < 0.5) z = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
	else z = sqrt(-2.0 * log(u1)) * sin(2.0 * PI * u2);
	return z * sigma + mu;
}

/*
 * main - Main program
 */
int main(int argc, char **argv)
{
	cout << "GENINST - Generate an instance of Minimum Cost List Coloring Problem." << endl;

	if (argc < 3) {
		cout << "Usage: geninst file vertices" << endl;
		cout << "The files file.graph, file.cost and file.list are generated." << endl;
		bye("Bye!");
	}

	char *filename = argv[1];
	int vertices = atoi(argv[2]);
	if (vertices < 4 || vertices > 10000) bye("Number of vertices out range!");

	int clique_size = vertices * (vertices - 1) / 2;
	int *edge_u = new int[clique_size];
	int *edge_v = new int[clique_size];

	/* generate the random graph and the number of colors */
	int edges = 0;
	for (int u = 0; u < vertices - 1; u++) {
		for (int v = u + 1; v < vertices; v++) {
			if (urnd(0, 1, false) <= COEF_DENSITY) {
				edge_u[edges] = u;
				edge_v[edges] = v;
				edges++;
			}
		}
	}
	int colors = (int)((float)edges / EDGES_PER_COLOR);

	cout << "Statistics:" << endl;
	float density = 100.0 * (float)edges / (float)clique_size;
	cout << "  |V| = " << vertices << ", |E| = " << edges << " (density = " << density << "%), |C| = " << colors << "." << endl;

	/* write graph file */
	char filename_extension[210];
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".graph");
	FILE *stream = fopen(filename_extension, "wt");
	if (!stream) bye("Graph file cannot be created");
	fprintf(stream, "%d:%d\n", vertices, edges);
	for (int e = 0; e < edges; e++) fprintf(stream, "%d,%d\n", edge_u[e], edge_v[e]);
	fclose(stream);

	/* free memory */
	delete[] edge_v;
	delete[] edge_u;

	/* generate costs and write them to file */
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".cost");
	stream = fopen(filename_extension, "wt");
	if (!stream) bye("Cost file cannot be created");
	fprintf(stream, "%d\n", colors);
	for (int k = 0; k < colors; k++) {
		int cost = (int)urnd(400, 560, false);
		if (cost < 480) cost = 480;
		fprintf(stream, "%d ", cost);
	}
	fprintf(stream, "\n");
	fclose(stream);

	/* generate vertices of graph Gk */
	float prom = log((double)colors) / LOG_COLORS_PER_AVERAGE;
	float sigma = prom / AVERAGES_PER_SIGMA;
	cout << "Using size of V(Gk) with prom = " << prom << ", sigma = " << sigma << "." << endl;
	int *L_size = new int[vertices];
	int **L_set = new int*[vertices];
	int max_colors_available = colors / COLORS_PER_LIST_SIZE;
	for (int v = 0; v < vertices; v++) {
		L_size[v] = 0;
		L_set[v] = new int[max_colors_available];
	}
	for (int k = 0; k < colors; k++) {
		float C_size = grnd(prom, sigma);
		if (C_size < 3.0) C_size = 3.0;
		float prob = C_size / (float)vertices;
		int count = 0;
		for (int v = 0; v < vertices; v++) {
			if (urnd(0, 1, false) <= prob) {
				if (L_size[v] >= max_colors_available) bye("Out of available colors!");
				L_set[v][L_size[v]] = k;
				L_size[v]++;
				count++;
			}
		}
		if (count == 0) {
			/* force an addition of at least 1 vertex to the color */
			int v = (int)urnd(0, vertices, true);
			if (L_size[v] >= max_colors_available) bye("Out of available colors!");
			L_set[v][L_size[v]] = k;
			L_size[v]++;
		}
	}

	/* write them to list file */
	strncpy(filename_extension, filename, 200);
	strcat(filename_extension, ".list");
	stream = fopen(filename_extension, "wt");
	if (!stream) bye("List file cannot be created");
	fprintf(stream, "%d:%d\n", vertices, colors);
	for (int v = 0; v < vertices; v++) {
		int lv = L_size[v];
		fprintf(stream, "%d  ", lv);
		for (int s = 0; s < lv; s++) fprintf(stream, "%d ", L_set[v][s]);
		fprintf(stream, "\n");
	}
	fclose(stream);

	/* free memory */
	for (int v = 0; v < vertices; v++) delete[] L_set[v];
	delete[] L_set;
	delete[] L_size;
	return 0;
}
